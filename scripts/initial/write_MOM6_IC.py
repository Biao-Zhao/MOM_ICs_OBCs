#!/usr/bin/env python3
"""
Test version of write_MOM6_IC.py.

Differences from the original:
1. Deep trailing NaNs are filled with a vectorized ffill along zl.
2. Earth-relative u and v are regridded directly to the MOM6 U and V
   locations instead of first being regridded to the full supergrid.
3. MOM6 state variables are stored as float32, matching the precision of
   the source GLORYS fields.

The tracer path is unchanged: temperature, salinity, and SSH are regridded
directly to the MOM6 tracer grid.
"""

import argparse
import gc
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor
from time import perf_counter

import numpy as np
import xarray
import xesmf
import yaml

from HCtFlood import kara as flood

from depths import vgrid_to_layers


_KARA_SOURCE_PLANES = None
TEMP_SALT_VALIDITY_VERSION = "surface_bfill_v1"


def _flood_kara_plane(index):
    """Flood one 2-D plane with the unchanged HCTFlood Kara kernel."""
    plane = _KARA_SOURCE_PLANES[index]
    return index, flood.flood_kara_ma(np.ma.masked_invalid(plane))


def flood_kara_parallel_levels(
    data,
    workers,
    xdim="lon",
    ydim="lat",
    zdim="z",
    tdim="time",
):
    """Run the unchanged Kara algorithm concurrently across vertical levels."""
    if workers <= 1:
        return flood.flood_kara(
            data,
            xdim=xdim,
            ydim=ydim,
            zdim=zdim,
            tdim=tdim,
        )

    if "fork" not in multiprocessing.get_all_start_methods():
        raise RuntimeError(
            "Parallel Kara flooding requires the multiprocessing 'fork' "
            "start method available on Linux."
        )

    if tdim not in data.dims:
        data = data.expand_dims(dim=tdim)
    if zdim not in data.dims:
        data = data.expand_dims(dim=zdim)

    ordered = data.transpose(tdim, zdim, ydim, xdim)
    values = np.asarray(ordered.data)
    nrec, nlev, ny, nx = values.shape
    planes = values.reshape((-1, ny, nx))
    output = np.empty_like(planes)

    # Compile the Numba kernel before forking so every worker inherits it.
    warmup = np.ma.masked_invalid(
        np.array([[1.0, 1.0], [1.0, np.nan]], dtype=values.dtype)
    )
    flood.flood_kara_ma(warmup)

    global _KARA_SOURCE_PLANES
    _KARA_SOURCE_PLANES = planes
    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
    ) as executor:
        for index, flooded_plane in executor.map(
            _flood_kara_plane,
            range(len(planes)),
            chunksize=1,
        ):
            output[index] = flooded_plane
    _KARA_SOURCE_PLANES = None

    output = output.reshape((nrec, nlev, ny, nx))
    return xarray.DataArray(
        output,
        name=str(data.name),
        coords={
            tdim: data[tdim],
            zdim: data[zdim],
            ydim: data[ydim],
            xdim: data[xdim],
        },
        dims=(tdim, zdim, ydim, xdim),
    )


def report_time(label, start_time):
    """Print elapsed wall-clock time for one processing phase."""
    elapsed = perf_counter() - start_time
    print(f"[TIMING] {label}: {elapsed:.2f} seconds", flush=True)
    return perf_counter()


def make_regridder(source, target, filename, reuse_weights, periodic):
    """
    Reuse a compatible weight file when requested and present.

    Direct U/V regridding uses new weight filenames, so the first run creates
    them even when reuse_weights is True. Subsequent runs reuse them.
    Bilinear interpolation is retained inside the source domain; destination
    points outside it use the nearest source point.
    """
    reuse_this_file = bool(reuse_weights and os.path.exists(filename))
    action = "Reusing" if reuse_this_file else "Generating"
    print(f"{action} regridding weights: {filename}")

    return xesmf.Regridder(
        source,
        target,
        method="bilinear",
        filename=filename,
        reuse_weights=reuse_this_file,
        periodic=periodic,
        extrap_method="nearest_s2d",
    )


def write_initial(config):
    total_start = perf_counter()

    temp_file = config["glorys_temperature"]
    sal_file = config["glorys_salinity"]
    ssh_file = config["glorys_sea_surface_height"]
    u_file = config["glorys_zonal_velocity"]
    v_file = config["glorys_meridional_velocity"]
    hour = int(str(config["ssh_time"]).zfill(2))

    vgrid_file = config["vgrid_file"]
    grid_file = config["grid_file"]
    output_file = config["output_file"]
    resolution = str(config["resolution"])
    weight_dir = config.get("weight_dir", ".")
    os.makedirs(weight_dir, exist_ok=True)
    reuse_weights = config.get("reuse_weights", False)
    kara_workers = int(config.get("kara_workers", 1))
    if kara_workers < 1:
        raise ValueError("kara_workers must be at least 1.")

    region_keys = ("min_lon", "max_lon", "min_lat", "max_lat")
    region_values = [config.get(key) for key in region_keys]
    if any(value is not None for value in region_values) and not all(
        value is not None for value in region_values
    ):
        raise ValueError(
            "Set all of min_lon, max_lon, min_lat, and max_lat, "
            "or leave all four unset for global input."
        )

    region_selection = {}
    if all(value is not None for value in region_values):
        lon_min, lon_max, lat_min, lat_max = map(float, region_values)
        region_selection = {
            "longitude": slice(lon_min, lon_max),
            "latitude": slice(lat_min, lat_max),
        }
        periodic_source = False
    else:
        periodic_source = True

    variable_names = config["variable_names"]
    temp_var = variable_names["temperature"]
    sal_var = variable_names["salinity"]
    ssh_var = variable_names["sea_surface_height"]
    u_var = variable_names["zonal_velocity"]
    v_var = variable_names["meridional_velocity"]

    print("Reading from the following GLORYS files:")
    print(f"  Temperature: {temp_file}")
    print(f"  Salinity:    {sal_file}")
    print(f"  SSH:         {ssh_file}")
    print(f"  U (zonal):   {u_file}")
    print(f"  V (merid.):  {v_file}")
    print(f"  SSH at hour: {hour}")
    if region_selection:
        print(
            f"  Region: lon=[{lon_min}, {lon_max}], "
            f"lat=[{lat_min}, {lat_max}]"
        )
    else:
        print("  Region: global (using the complete input grid)")

    ds_temp = (
        xarray.open_dataset(temp_file)[temp_var]
        .sel(region_selection)
        .rename({"longitude": "lon", "latitude": "lat"})
    )

    ds_sal = (
        xarray.open_dataset(sal_file)[sal_var]
        .sel(region_selection)
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    ds_ssh = (
        xarray.open_dataset(ssh_file)[ssh_var]
        .sel(region_selection)
        .isel(time=hour)
        .isel(depth=0, drop=True)
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    ds_u = (
        xarray.open_dataset(u_file)[u_var]
        .sel(region_selection)
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    ds_v = (
        xarray.open_dataset(v_file)[v_var]
        .sel(region_selection)
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    vgrid_ds = xarray.open_dataset(vgrid_file)
    vgrid = vgrid_ds["dz"]
    z = vgrid_to_layers(vgrid)
    ztarget = xarray.DataArray(
        z,
        name="zl",
        dims=["zl"],
        coords={"zl": z},
    )

    glorys = xarray.merge([ds_temp, ds_sal, ds_ssh, ds_u, ds_v])
    print("GLORYS dimensions:", glorys.dims)
    print(f"Kara vertical-level workers: {kara_workers}")

    # Keep the time treatment used by the original script.
    glorys["time"] = (("time",), ds_temp["time"].dt.floor("1d").data)

    target_grid = xarray.open_dataset(grid_file)
    angle_variable = target_grid["angle_dx"]
    angle_units = str(angle_variable.attrs.get("units", "")).lower()
    if "degree" in angle_units:
        print("Converting MOM6 angle_dx from degrees to radians")
        angle_to_radians = np.pi / 180.0
    elif "radian" in angle_units:
        angle_to_radians = 1.0
    else:
        maximum_angle = float(
            np.abs(angle_variable).max(skipna=True).values
        )
        if not np.isfinite(maximum_angle) or maximum_angle > 2.0 * np.pi:
            raise ValueError(
                "MOM6 angle_dx units are missing or unsupported, and its "
                "values are not consistent with radians"
            )
        print(
            "[WARNING] MOM6 angle_dx units are not specified; "
            "assuming radians"
        )
        angle_to_radians = 1.0

    # Tracer points: odd/odd locations on the MOM6 supergrid.
    target_t = (
        target_grid[["x", "y"]]
        .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
        .rename({"x": "lon", "y": "lat", "nxp": "xh", "nyp": "yh"})
    )

    # U points: even/odd locations on the MOM6 supergrid.
    target_u_native = target_grid[["x", "y", "angle_dx"]].isel(
        nxp=slice(0, None, 2),
        nyp=slice(1, None, 2),
    )
    target_u = target_u_native[["x", "y"]].rename(
        {"x": "lon", "y": "lat", "nxp": "xq", "nyp": "yh"}
    )
    angle_u = target_u_native["angle_dx"].rename(
        {"nxp": "xq", "nyp": "yh"}
    ) * angle_to_radians
    angle_u.attrs["units"] = "radians"

    # V points: odd/even locations on the MOM6 supergrid.
    target_v_native = target_grid[["x", "y", "angle_dx"]].isel(
        nxp=slice(1, None, 2),
        nyp=slice(0, None, 2),
    )
    target_v = target_v_native[["x", "y"]].rename(
        {"x": "lon", "y": "lat", "nxp": "xh", "nyp": "yq"}
    )
    angle_v = target_v_native["angle_dx"].rename(
        {"nxp": "xh", "nyp": "yq"}
    ) * angle_to_radians
    angle_v.attrs["units"] = "radians"

    print("Tracer target dimensions:", target_t.dims)
    print("U target dimensions:", target_u.dims)
    print("V target dimensions:", target_v.dims)

    phase_start = perf_counter()
    glorys_to_t = make_regridder(
        glorys,
        target_t,
        os.path.join(
        weight_dir,
        f"regrid_glorys_{resolution}_tracers.nc",
        ),
        reuse_weights,
        periodic_source,
    )
    glorys_to_u = make_regridder(
        glorys,
        target_u,
        os.path.join(
        weight_dir,
        f"regrid_glorys_{resolution}_u.nc",
        ),
        reuse_weights,
        periodic_source,
    )
    glorys_to_v = make_regridder(
        glorys,
        target_v,
        os.path.join(
        weight_dir,
        f"regrid_glorys_{resolution}_v.nc",
        ),
        reuse_weights,
        periodic_source,
    )
    report_time("create/reuse regridders", phase_start)

    output_time = ds_temp["time"].dt.floor("1d").data

    def flood_3d_source(source, label):
        """Vertically interpolate and flood one source field."""
        phase = perf_counter()
        source = source.assign_coords(time=output_time)
        # Extend the shallowest source value to MOM6 layers above the first
        # GLORYS level, and extend each water column's deepest valid value
        # downward before horizontal flooding.  Without the downward fill,
        # abyssal levels contain only a few trench points on a global grid;
        # HCTFlood then tries to propagate those values across all 4320
        # longitudes and can exceed its hard-coded 1000-iteration limit.
        reverted = (
            source.interp(depth=ztarget)
            .bfill("zl")
            .ffill("zl")
        )
        flooded_source = flood_kara_parallel_levels(
            reverted,
            workers=kara_workers,
            zdim="zl",
        )
        report_time(f"{label}: vertical interpolation + Kara fill", phase)
        return flooded_source

    def finish_3d_target(field, label):
        """Fill deep target NaNs and materialize one final field."""
        print(f"Filling deep NaNs for {label}...")
        original_dims = field.dims
        field = field.ffill("zl").transpose(*original_dims)
        field.load()
        return field

    # zl is shallow-to-deep, so ffill extends the deepest valid value
    # downward without one Python call per water column.
    zl_values = ztarget["zl"].values
    if len(zl_values) > 1 and not np.all(np.diff(zl_values) > 0):
        raise ValueError(
            "Expected zl to increase monotonically from shallow to deep."
        )
    validity_file = os.path.join(
        weight_dir,
        "temp_salt_valid_points.nc",
    )
    reuse_validity = False
    if os.path.isfile(validity_file):
        try:
            with xarray.open_dataset(validity_file) as existing:
                existing_mask = existing["temp_salt_valid"]
                expected_shape = (
                    len(zl_values),
                    target_t.sizes["yh"],
                    target_t.sizes["xh"],
                )
                reuse_validity = (
                    existing_mask.shape == expected_shape
                    and "zl" in existing.coords
                    and np.allclose(existing["zl"].values, zl_values)
                    and existing_mask.attrs.get("mask_version")
                    == TEMP_SALT_VALIDITY_VERSION
                )
        except (KeyError, OSError, ValueError):
            reuse_validity = False
    make_validity = not reuse_validity
    if make_validity:
        print("Creating density-gradient validity mask:", validity_file)
        phase = perf_counter()
        # Treat the thin target layers above the shallowest source level as
        # reliable: they inherit the first vertically valid source value and
        # introduce no new horizontal structure.  Do not ffill here, so
        # bottom extrapolation and horizontal flooding remain excluded from
        # density-gradient calculations.
        temp_valid = (
            ds_temp.assign_coords(time=output_time)
            .interp(depth=ztarget)
            .bfill("zl")
            .notnull()
        )
        salt_valid = (
            ds_sal.assign_coords(time=output_time)
            .interp(depth=ztarget)
            .bfill("zl")
            .notnull()
        )
        valid_fraction = glorys_to_t(
            (temp_valid & salt_valid).astype(np.float32)
        )
        temp_salt_valid = valid_fraction >= 0.999
        if "time" in temp_salt_valid.dims:
            temp_salt_valid = temp_salt_valid.all("time")
        temp_salt_valid = temp_salt_valid.astype(np.uint8)
        temp_salt_valid.name = "temp_salt_valid"
        temp_salt_valid.attrs = {
            "long_name": (
                "Temperature and salinity reliable for horizontal "
                "density gradients"
            ),
            "flag_values": np.array([0, 1], dtype=np.uint8),
            "flag_meanings": "excluded_from_gradient valid_for_gradient",
            "valid_fraction_threshold": np.float32(0.999),
            "mask_version": TEMP_SALT_VALIDITY_VERSION,
            "vertical_fill_policy": (
                "surface bfill allowed; bottom ffill excluded"
            ),
        }
        validity_partial = validity_file + ".partial"
        if os.path.exists(validity_partial):
            os.remove(validity_partial)
        temp_salt_valid.to_netcdf(
            validity_partial,
            format="NETCDF4",
            engine="netcdf4",
            encoding={
                "temp_salt_valid": {
                    "dtype": "uint8",
                    "_FillValue": None,
                    "zlib": True,
                    "complevel": 1,
                    "shuffle": True,
                }
            },
        )
        os.replace(validity_partial, validity_file)
        report_time("create pre-fill validity mask", phase)
        del temp_valid, salt_valid, valid_fraction, temp_salt_valid
        gc.collect()
    else:
        print("Reusing pre-fill validity mask:", validity_file)

    # Process and materialize each tracer independently so the global source
    # interpolation/flood graphs are released before the next field begins.
    flooded_temp = flood_3d_source(ds_temp, "temp")
    temp = finish_3d_target(
        glorys_to_t(flooded_temp).astype(np.float32), "temp"
    )
    temp.name = "temp"
    del flooded_temp
    gc.collect()

    flooded_sal = flood_3d_source(ds_sal, "salt")
    salt = finish_3d_target(
        glorys_to_t(flooded_sal).astype(np.float32), "salt"
    )
    salt.name = "salt"
    del flooded_sal
    gc.collect()

    phase_start = perf_counter()
    flooded_ssh = flood.flood_kara(ds_ssh)
    ssh = flooded_ssh.isel(z=0).drop_vars("z")
    if "time" in ssh.dims:
        ssh = ssh.assign_coords(time=output_time)
    else:
        ssh = ssh.expand_dims(time=output_time)
    ssh = glorys_to_t(ssh).astype(np.float32)
    ssh.load()
    ssh.name = "ssh"
    report_time("ssh: flood + regrid + compute", phase_start)
    del flooded_ssh
    gc.collect()

    # Velocity rotation needs both earth-relative components. Build both
    # target fields first, then materialize them in one shared Dask compute so
    # the source u/v interpolation and flooding graphs execute only once.
    flooded_u = flood_3d_source(ds_u, "earth u")
    flooded_v = flood_3d_source(ds_v, "earth v")
    velocity_source = xarray.merge((flooded_u, flooded_v))

    phase_start = perf_counter()
    earth_at_u = glorys_to_u(velocity_source)
    uo = (
        np.cos(angle_u) * earth_at_u[u_var]
        + np.sin(angle_u) * earth_at_u[v_var]
    ).astype(np.float32)
    uo.name = "u"
    earth_at_v = glorys_to_v(velocity_source)
    vo = (
        -np.sin(angle_v) * earth_at_v[u_var]
        + np.cos(angle_v) * earth_at_v[v_var]
    ).astype(np.float32)
    vo.name = "v"

    print("Filling deep NaNs for u and v...")
    uo = uo.ffill("zl").transpose(*uo.dims)
    vo = vo.ffill("zl").transpose(*vo.dims)
    velocity_target = xarray.merge((uo, vo))
    velocity_target.load()
    uo = velocity_target["u"]
    vo = velocity_target["v"]
    report_time("shared U/V regrid + rotation + compute", phase_start)
    del earth_at_u, earth_at_v, velocity_source, flooded_u, flooded_v
    gc.collect()

    interped = xarray.merge((temp, salt, ssh, uo, vo)).transpose(
        "time", "zl", "yh", "yq", "xh", "xq"
    )

    xh_1d = target_grid["x"].isel(
        nxp=slice(1, None, 2), nyp=0
    ).values
    yh_1d = target_grid["y"].isel(
        nxp=0, nyp=slice(1, None, 2)
    ).values
    xq_1d = target_grid["x"].isel(
        nxp=slice(0, None, 2), nyp=0
    ).values
    yq_1d = target_grid["y"].isel(
        nxp=0, nyp=slice(0, None, 2)
    ).values

    interped = interped.assign_coords(
        {
            "xh": ("xh", xh_1d),
            "yh": ("yh", yh_1d),
            "xq": ("xq", xq_1d),
            "yq": ("yq", yq_1d),
        }
    )

    all_vars = list(interped.data_vars) + list(interped.coords)
    encodings = {name: {"_FillValue": None} for name in all_vars}
    fill_value = np.float32(1.0e20)

    for var in ["temp", "salt", "u", "v", "ssh"]:
        encodings[var]["dtype"] = "float32"

    encodings["temp"]["_FillValue"] = fill_value
    encodings["salt"]["_FillValue"] = fill_value
    encodings["time"].update(
        {"dtype": "float64", "calendar": "gregorian"}
    )

    interped["zl"].attrs = {
        "long_name": "Layer pseudo-depth, -z*",
        "units": "meter",
        "cartesian_axis": "Z",
        "positive": "down",
    }

    output_folder = os.path.dirname(output_file)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    print("Variables in final dataset:", list(interped.data_vars))
    print("Writing:", output_file)

    phase_start = perf_counter()
    interped.to_netcdf(
        output_file,
        format="NETCDF4",
        engine="netcdf4",
        encoding=encodings,
    )
    report_time("NetCDF write", phase_start)
    report_time("TOTAL", total_start)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MOM6 ICs with direct U/V-grid regridding."
    )
    parser.add_argument(
        "--config_file",
        type=str,
        default="glorys_ic.yaml",
        help="Path to the YAML configuration file.",
    )
    args = parser.parse_args()

    with open(args.config_file, "r", encoding="utf-8") as yaml_file:
        config = yaml.safe_load(yaml_file)

    required = [
        "glorys_temperature",
        "glorys_salinity",
        "glorys_sea_surface_height",
        "glorys_zonal_velocity",
        "glorys_meridional_velocity",
        "ssh_time",
        "vgrid_file",
        "grid_file",
        "output_file",
        "resolution",
        "variable_names",
    ]
    missing = [key for key in required if key not in config]
    if missing:
        parser.error(
            "Missing required configuration keys: " + ", ".join(missing)
        )

    write_initial(config)


if __name__ == "__main__":
    main()

