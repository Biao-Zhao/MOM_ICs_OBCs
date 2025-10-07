#!/usr/bin/env python3

import xarray as xr
import numpy as np
from pathlib import Path
import sys
import yaml

# ============================================================================
# 1) Parse config path
if len(sys.argv) < 3 or sys.argv[1] != "--config":
    sys.exit("Usage: python3 merge_glorys_nc.py --config <yaml_file>")
cfg_file = sys.argv[2]

with open(cfg_file, "r") as f:
    cfg = yaml.safe_load(f) or {}

for k in ["thetao_fn", "so_fn", "uovo_fn", "ssh_fn", "merged_fn", "ssh_time"]:
    if k not in cfg:
        sys.exit(f"[ERROR] Missing '{k}' in YAML: {cfg_file}")

thetao_fn   = Path(cfg["thetao_fn"])
so_fn       = Path(cfg["so_fn"])
uovo_fn     = Path(cfg["uovo_fn"])
ssh_fn      = Path(cfg["ssh_fn"])
merged_fn = Path(cfg["merged_fn"])
hour    = int(str(cfg["ssh_time"]).zfill(2))  # e.g. "06" -> 6

bbox_keys = ["min_lon", "max_lon", "min_lat", "max_lat"]
have_bbox = all(k in cfg and cfg[k] is not None for k in bbox_keys)

# =============================================================================
# 2) Determine integer‐index slice from thetao’s native grid
# =============================================================================
ds_thetao = xr.open_dataset(thetao_fn, decode_times=True, mask_and_scale=True)
lat_vals  = ds_thetao.latitude.values
lon_vals  = ds_thetao.longitude.values

if have_bbox:
    lat_bounds = (float(cfg["min_lat"]), float(cfg["max_lat"]))
    lon_bounds = (float(cfg["min_lon"]), float(cfg["max_lon"]))
else:
    lat_bounds = (float(lat_vals[0]), float(lat_vals[-1]))
    lon_bounds = (float(lon_vals[0]), float(lon_vals[-1]))

ilat0 = np.searchsorted(lat_vals, lat_bounds[0], side="left")
ilat1 = np.searchsorted(lat_vals, lat_bounds[1], side="right")
ilon0 = np.searchsorted(lon_vals, lon_bounds[0], side="left")
ilon1 = np.searchsorted(lon_vals, lon_bounds[1], side="right")

#print(f"Latitude index slice:  {ilat0} … {ilat1-1}  → "
#      f"{lat_vals[ilat0]:.6f}° … {lat_vals[ilat1-1]:.6f}°")
#print(f"Longitude index slice: {ilon0} … {ilon1-1}  → "
#      f"{lon_vals[ilon0]:.6f}° … {lon_vals[ilon1-1]:.6f}°")

# =============================================================================
# 3) Subset thetao and capture its “true” coords
# =============================================================================
ds_thetao_sub = ds_thetao.isel(
    latitude = slice(ilat0, ilat1),
    longitude= slice(ilon0, ilon1)
)
lat_grid = ds_thetao_sub.latitude
lon_grid = ds_thetao_sub.longitude

# =============================================================================
# 4) Subset + re‐assign coords for so and uovo
# =============================================================================
ds_so_sub = (
    xr.open_dataset(so_fn, decode_times=True, mask_and_scale=True)
      .isel(latitude = slice(ilat0, ilat1),
            longitude= slice(ilon0, ilon1))
      .assign_coords(latitude=lat_grid, longitude=lon_grid)
)

ds_uovo_sub = (
    xr.open_dataset(uovo_fn, decode_times=True, mask_and_scale=True)
      .isel(latitude = slice(ilat0, ilat1),
            longitude= slice(ilon0, ilon1))
      .assign_coords(latitude=lat_grid, longitude=lon_grid)
)

# =============================================================================
# 5) Subset SSH: use ssh_time from YAML to pick the SSH at the specified hour, 
#                then assign coords & time
# =============================================================================
ds_ssh_raw = xr.open_dataset(ssh_fn, decode_times=False, mask_and_scale=True)

zos = (
    ds_ssh_raw["sea_surface_height"]
      .isel(time=hour, depth=0, drop=True)
      .rename("zos")
)

ds_ssh_sub = (
    zos.isel(latitude = slice(ilat0, ilat1),
             longitude= slice(ilon0, ilon1))
       .expand_dims(time=1)
       .assign_coords(time=ds_thetao_sub.time,
                      latitude=lat_grid,
                      longitude=lon_grid)
)    
# =============================================================================
# 6) Merge and write out final NetCDF
# =============================================================================
ds_combined = xr.merge([
    ds_thetao_sub,
    ds_so_sub,
    ds_uovo_sub,
    ds_ssh_sub
])

merged_fn.parent.mkdir(parents=True, exist_ok=True)
ds_combined.to_netcdf(merged_fn)
print(f"Wrote merged file → {merged_fn}")
