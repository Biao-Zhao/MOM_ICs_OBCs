#!/usr/bin/env python3
"""Create a simple two-dimensional SIS2 initial-condition file.

The CMEMS daily sea-ice product supplies ice concentration and the mean ice
thickness over the ice-covered part of each source cell.  Concentration and
grid-cell-equivalent ice volume are conservatively remapped to the MOM6/SIS2
tracer grid.  The output contains the two fields used by SIS2's official
``file`` initialization path:

``aice``
    Total sea-ice concentration (0--1).
``hm``
    Ice mass per unit ice-covered area (kg m-2).
"""

import argparse
import os
from pathlib import Path

import numpy as np
import xarray as xr


def coordinate_bounds(values, *, latitude=False):
    """Return cell-edge coordinates for a monotonic 1-D center coordinate."""
    centers = np.asarray(values, dtype=np.float64)
    if centers.ndim != 1 or centers.size < 2:
        raise ValueError("Source coordinates must be one-dimensional with at least two points")
    delta = np.diff(centers)
    if not np.all(delta > 0.0):
        raise ValueError("Source coordinates must increase monotonically")

    bounds = np.empty(centers.size + 1, dtype=np.float64)
    bounds[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    bounds[0] = centers[0] - 0.5 * delta[0]
    bounds[-1] = centers[-1] + 0.5 * delta[-1]
    if latitude:
        bounds = np.clip(bounds, -90.0, 90.0)
    return bounds


def source_grid(dataset):
    """Build the rectilinear xESMF source grid, including cell bounds."""
    lon = dataset["longitude"].astype(np.float64)
    lat = dataset["latitude"].astype(np.float64)
    return xr.Dataset(
        coords={
            "lon": ("lon", lon.values),
            "lat": ("lat", lat.values),
            "lon_b": ("lon_b", coordinate_bounds(lon.values)),
            "lat_b": (
                "lat_b",
                coordinate_bounds(lat.values, latitude=True),
            ),
        }
    )


def target_tracer_grid(hgrid):
    """Build tracer centers and tracer-cell corners from a MOM6 supergrid."""
    required = {"x", "y"}
    missing = required.difference(hgrid.variables)
    if missing:
        raise ValueError(f"MOM6 hgrid is missing variables: {sorted(missing)}")

    centers = hgrid[["x", "y"]].isel(
        nxp=slice(1, None, 2),
        nyp=slice(1, None, 2),
    )
    corners = hgrid[["x", "y"]].isel(
        nxp=slice(0, None, 2),
        nyp=slice(0, None, 2),
    )
    return xr.Dataset(
        coords={
            "lon": (("yh", "xh"), centers["x"].values),
            "lat": (("yh", "xh"), centers["y"].values),
            "lon_b": (("y_b", "x_b"), corners["x"].values),
            "lat_b": (("y_b", "x_b"), corners["y"].values),
        }
    )


def logical_coordinates(hgrid):
    """Return the one-dimensional MOM6 tracer-grid logical coordinates."""
    xh = hgrid["x"].isel(nxp=slice(1, None, 2), nyp=0).values
    yh = hgrid["y"].isel(nxp=0, nyp=slice(1, None, 2)).values
    return np.asarray(xh, dtype=np.float64), np.asarray(yh, dtype=np.float64)


def make_regridder(source, target, weight_file, reuse_weights):
    """Create or reuse conservative CMEMS-to-SIS2 remapping weights."""
    try:
        import xesmf
    except ImportError as error:
        raise RuntimeError(
            "write_SIS2_IC.py requires xesmf and an ESMF-compatible environment"
        ) from error

    reuse = bool(reuse_weights and os.path.exists(weight_file))
    action = "Reusing" if reuse else "Generating"
    print(f"{action} conservative regridding weights: {weight_file}")
    return xesmf.Regridder(
        source,
        target,
        method="conservative",
        filename=weight_file,
        reuse_weights=reuse,
        periodic=False,
        ignore_degenerate=True,
        unmapped_to_nan=True,
    )


def write_sis2_initial(args):
    """Conservatively remap CMEMS ice fields and write a SIS2 IC file."""
    input_path = Path(args.input)
    grid_path = Path(args.grid)
    output_path = Path(args.output)
    weight_dir = Path(args.weight_dir)

    if not input_path.is_file():
        raise FileNotFoundError(f"CMEMS sea-ice file not found: {input_path}")
    if not grid_path.is_file():
        raise FileNotFoundError(f"MOM6 hgrid file not found: {grid_path}")
    weight_dir.mkdir(parents=True, exist_ok=True)

    region_values = (args.min_lon, args.max_lon, args.min_lat, args.max_lat)
    regional = all(value is not None for value in region_values)
    if regional:
        region_selection = {
            "longitude": slice(args.min_lon, args.max_lon),
            "latitude": slice(args.min_lat, args.max_lat),
        }
        region_label = (
            f"lon{args.min_lon:g}_{args.max_lon:g}_"
            f"lat{args.min_lat:g}_{args.max_lat:g}"
        ).replace("-", "m").replace(".", "p")
        print(
            f"Source region: lon=[{args.min_lon}, {args.max_lon}], "
            f"lat=[{args.min_lat}, {args.max_lat}]"
        )
    else:
        region_selection = {}
        region_label = "global"
        print("Source region: global (using the complete input grid)")

    with xr.open_dataset(input_path) as source_ds, xr.open_dataset(grid_path) as hgrid:
        required = {"siconc", "sithick"}
        missing = required.difference(source_ds.variables)
        if missing:
            raise ValueError(f"CMEMS sea-ice file is missing variables: {sorted(missing)}")
        if source_ds.sizes.get("time", 1) != 1:
            raise ValueError("Expected exactly one daily-mean sea-ice time record")

        selected_source = source_ds.sel(region_selection)
        if (
            selected_source.sizes.get("longitude", 0) < 2
            or selected_source.sizes.get("latitude", 0) < 2
        ):
            raise ValueError("The selected sea-ice source region is empty or too small")

        concentration = selected_source["siconc"]
        thickness = selected_source["sithick"]
        if "time" in concentration.dims:
            concentration = concentration.isel(time=0, drop=True)
            thickness = thickness.isel(time=0, drop=True)

        concentration = concentration.fillna(0.0).clip(0.0, 1.0)
        thickness = thickness.fillna(0.0).clip(min=0.0)
        ice_volume = concentration * thickness

        src_grid = source_grid(selected_source)
        dst_grid = target_tracer_grid(hgrid)
        if regional:
            weight_name = (
                f"regrid_glorys_{args.resolution}_ice_conservative_"
                f"{region_label}.nc"
            )
        else:
            weight_name = f"regrid_glorys_{args.resolution}_ice_conservative.nc"
        weight_file = weight_dir / weight_name
        regridder = make_regridder(
            src_grid,
            dst_grid,
            str(weight_file),
            args.reuse_weights,
        )

        remapped_concentration = regridder(concentration).fillna(0.0)
        remapped_volume = regridder(ice_volume).fillna(0.0)
        remapped_concentration = remapped_concentration.clip(0.0, 1.0)
        remapped_volume = remapped_volume.clip(min=0.0)

        wet_ice = remapped_concentration >= args.minimum_concentration
        mean_thickness = xr.where(
            wet_ice,
            remapped_volume / remapped_concentration,
            0.0,
        )
        aice = xr.where(wet_ice, remapped_concentration, 0.0).astype(np.float32)
        hm = (args.ice_density * mean_thickness).astype(np.float32)

        xh, yh = logical_coordinates(hgrid)
        aice = aice.assign_coords(xh=("xh", xh), yh=("yh", yh))
        hm = hm.assign_coords(xh=("xh", xh), yh=("yh", yh))
        aice.name = "aice"
        hm.name = "hm"
        aice.attrs = {
            "long_name": "total sea ice concentration",
            "units": "1",
            "valid_min": np.float32(0.0),
            "valid_max": np.float32(1.0),
        }
        hm.attrs = {
            "long_name": "sea ice mass per unit ice-covered area",
            "units": "kg m-2",
        }

        output = xr.Dataset({"aice": aice, "hm": hm})
        source_time = source_ds.get("time")
        if source_time is not None:
            output.attrs["source_time"] = str(source_time.values[0])
        output.attrs.update(
            {
                "title": "SIS2 two-dimensional initial condition",
                "source_file": input_path.name,
                "source_product": "CMEMS daily-mean sea ice",
                "ice_density_kg_m3": args.ice_density,
                "minimum_ice_concentration": args.minimum_concentration,
                "history": "Conservative remapping of concentration and ice volume",
            }
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        encoding = {
            "aice": {"dtype": "float32", "_FillValue": np.float32(1.0e20)},
            "hm": {"dtype": "float32", "_FillValue": np.float32(1.0e20)},
            "xh": {"_FillValue": None},
            "yh": {"_FillValue": None},
        }
        print(f"Writing SIS2 initial condition: {output_path}")
        output.to_netcdf(output_path, encoding=encoding)
        print(
            "Output ranges: "
            f"aice={float(aice.min()):.6g}..{float(aice.max()):.6g}, "
            f"hm={float(hm.min()):.6g}..{float(hm.max()):.6g} kg m-2"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a 2-D SIS2 initial condition from CMEMS sea ice"
    )
    parser.add_argument("--input", required=True, help="CMEMS daily sea-ice NetCDF")
    parser.add_argument("--grid", required=True, help="MOM6 ocean_hgrid.nc")
    parser.add_argument("--output", required=True, help="Output SIS2 initial NetCDF")
    parser.add_argument("--weight-dir", required=True, help="Regridding-weight directory")
    parser.add_argument("--resolution", required=True, help="Grid label used in weight filenames")
    parser.add_argument("--min-lon", type=float, help="Minimum source longitude")
    parser.add_argument("--max-lon", type=float, help="Maximum source longitude")
    parser.add_argument("--min-lat", type=float, help="Minimum source latitude")
    parser.add_argument("--max-lat", type=float, help="Maximum source latitude")
    parser.add_argument(
        "--ice-density",
        type=float,
        default=905.0,
        help="SIS2 nominal ice density in kg m-3 (default: 905)",
    )
    parser.add_argument(
        "--minimum-concentration",
        type=float,
        default=1.0e-6,
        help="Concentrations below this value are set to zero",
    )
    parser.add_argument(
        "--reuse-weights",
        action="store_true",
        help="Reuse an existing compatible conservative weight file",
    )
    args = parser.parse_args()
    if args.ice_density <= 0.0:
        parser.error("--ice-density must be positive")
    if not 0.0 <= args.minimum_concentration < 1.0:
        parser.error("--minimum-concentration must be in [0, 1)")
    region_values = (args.min_lon, args.max_lon, args.min_lat, args.max_lat)
    if any(value is not None for value in region_values) and not all(
        value is not None for value in region_values
    ):
        parser.error(
            "set all of --min-lon, --max-lon, --min-lat, and --max-lat, "
            "or omit all four for global input"
        )
    if all(value is not None for value in region_values):
        if args.min_lon >= args.max_lon:
            parser.error("--min-lon must be smaller than --max-lon")
        if args.min_lat >= args.max_lat:
            parser.error("--min-lat must be smaller than --max-lat")
    return args


def main():
    write_sis2_initial(parse_args())


if __name__ == "__main__":
    main()
