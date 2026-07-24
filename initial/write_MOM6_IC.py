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
import os
from time import perf_counter

import numpy as np
import xarray
import xesmf
import yaml

from HCtFlood import kara as flood

from depths import vgrid_to_layers


def report_time(label, start_time):
    """Print elapsed wall-clock time for one processing phase."""
    elapsed = perf_counter() - start_time
    print(f"[TIMING] {label}: {elapsed:.2f} seconds", flush=True)
    return perf_counter()


def make_regridder(source, target, filename, reuse_weights):
    """
    Reuse a compatible weight file when requested and present.

    Direct U/V regridding uses new weight filenames, so the first run creates
    them even when reuse_weights is True. Subsequent runs reuse them.
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
        periodic=False,
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

    lon_min = float(config["min_lon"])
    lon_max = float(config["max_lon"])
    lat_min = float(config["min_lat"])
    lat_max = float(config["max_lat"])

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
    print(
        f"  Region: lon=[{lon_min}, {lon_max}], "
        f"lat=[{lat_min}, {lat_max}]"
    )

    ds_temp = (
        xarray.open_dataset(temp_file)[temp_var]
        .sel(
            longitude=slice(lon_min, lon_max),
            latitude=slice(lat_min, lat_max),
        )
        .rename({"longitude": "lon", "latitude": "lat"})
    )

    ds_sal = (
        xarray.open_dataset(sal_file)[sal_var]
        .sel(
            longitude=slice(lon_min, lon_max),
            latitude=slice(lat_min, lat_max),
        )
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    ds_ssh = (
        xarray.open_dataset(ssh_file)[ssh_var]
        .sel(
            longitude=slice(lon_min, lon_max),
            latitude=slice(lat_min, lat_max),
        )
        .isel(time=hour)
        .isel(depth=0, drop=True)
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    ds_u = (
        xarray.open_dataset(u_file)[u_var]
        .sel(
            longitude=slice(lon_min, lon_max),
            latitude=slice(lat_min, lat_max),
        )
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    ds_v = (
        xarray.open_dataset(v_file)[v_var]
        .sel(
            longitude=slice(lon_min, lon_max),
            latitude=slice(lat_min, lat_max),
        )
        .rename({"longitude": "lon", "latitude": "lat"})
        .assign_coords(lat=ds_temp.lat, lon=ds_temp.lon)
    )

    vgrid = xarray.open_dataarray(vgrid_file)
    z = vgrid_to_layers(vgrid)
    ztarget = xarray.DataArray(
        z,
        name="zl",
        dims=["zl"],
        coords={"zl": z},
    )

    glorys = xarray.merge([ds_temp, ds_sal, ds_ssh, ds_u, ds_v])
    print("GLORYS dimensions:", glorys.dims)

    # Keep the time treatment used by the original script.
    glorys["time"] = (("time",), ds_temp["time"].dt.floor("1d").data)

    # Interpolate vertically on the smaller GLORYS source grid.
    phase_start = perf_counter()
    revert = glorys.interp(depth=ztarget).bfill("zl")

    # Flood the four 3-D source fields over land.
    flooded = xarray.merge(
        flood.flood_kara(revert[var], zdim="zl")
        for var in [temp_var, sal_var, u_var, v_var]
    )

    # Flood SSH once, then select the surface added by flood_kara.
    flooded_ssh = flood.flood_kara(revert[ssh_var])
    surface_ssh = flooded_ssh.isel(z=0).drop_vars("z")
    surface_ssh["time"] = flooded.time
    flooded = xarray.merge(
        [flooded, surface_ssh.to_dataset(name=ssh_var)],
        compat="override",
    )
    report_time("vertical interpolation + flood", phase_start)

    target_grid = xarray.open_dataset(grid_file)

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
    )

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
    )

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
    )
    glorys_to_u = make_regridder(
        glorys,
        target_u,
        os.path.join(
        weight_dir,
        f"regrid_glorys_{resolution}_u.nc",
        ),
        reuse_weights,
    )
    glorys_to_v = make_regridder(
        glorys,
        target_v,
        os.path.join(
        weight_dir,
        f"regrid_glorys_{resolution}_v.nc",
        ),
        reuse_weights,
    )
    report_time("create/reuse regridders", phase_start)

    # Tracers already go directly to their final T locations.
    phase_start = perf_counter()
    interped_t = glorys_to_t(
        flooded[[temp_var, sal_var, ssh_var]]
    ).astype(np.float32)
    report_time("tracer regrid", phase_start)

    # Both earth-relative components are needed at U points because rotation
    # mixes u and v. Only the resulting model-relative u is retained.
    phase_start = perf_counter()
    earth_at_u = glorys_to_u(flooded[[u_var, v_var]])
    uo = (
        np.cos(angle_u) * earth_at_u[u_var]
        + np.sin(angle_u) * earth_at_u[v_var]
    ).astype(np.float32)
    uo.name = u_var
    report_time("U regrid + rotation", phase_start)

    # Both earth-relative components are also needed at V points. Only the
    # resulting model-relative v is retained.
    phase_start = perf_counter()
    earth_at_v = glorys_to_v(flooded[[u_var, v_var]])
    vo = (
        -np.sin(angle_v) * earth_at_v[u_var]
        + np.cos(angle_v) * earth_at_v[v_var]
    ).astype(np.float32)
    vo.name = v_var
    report_time("V regrid + rotation", phase_start)

    interped = xarray.merge((interped_t, uo, vo)).transpose(
        "time", "zl", "yh", "yq", "xh", "xq"
    )

    interped = interped.rename(
        {
            temp_var: "temp",
            sal_var: "salt",
            ssh_var: "ssh",
            u_var: "u",
            v_var: "v",
        }
    )

    # zl is shallow-to-deep, so ffill extends the deepest valid value
    # downward without one Python call per water column.
    zl_values = interped["zl"].values
    if len(zl_values) > 1 and not np.all(np.diff(zl_values) > 0):
        raise ValueError(
            "Expected zl to increase monotonically from shallow to deep."
        )

    phase_start = perf_counter()
    for var in ["temp", "salt", "u", "v"]:
        print(f"Filling deep NaNs for {var}...")
        original_dims = interped[var].dims
        interped[var] = interped[var].ffill("zl").transpose(*original_dims)
    report_time("vertical ffill", phase_start)

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
        "min_lon",
        "max_lon",
        "min_lat",
        "max_lat",
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

