# author: 'Jing Chen'
# description: 'Combine GLORYS PHY fields into one NetCDF file'
# created: '2025-08-05'
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# 1) Define file paths and geographic bounds
# =============================================================================
date_str   = "20240922"
input_dir  = Path(f"/scratch/cimes/bz5265/MOM_ICs_OBCs/CMEMS_data/{date_str}")
output_file = Path(
    f"/scratch/cimes/bz5265/MOM_ICs_OBCs/scripts/temp/{date_str}/"
    f"GLOBAL_ANALYSISFORECAST_PHY_{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}.nc"
)

thetao_fn = input_dir / f"glo12_rg_6h-i_{date_str}-00h_3D-thetao_hcst.nc"
so_fn     = input_dir / f"glo12_rg_6h-i_{date_str}-00h_3D-so_hcst.nc"
uovo_fn   = input_dir / f"glo12_rg_6h-i_{date_str}-00h_3D-uovo_hcst.nc"
ssh_fn    = input_dir / f"MOL_{date_str}.nc"

# desired coordinate bounds
lat_bounds = (10.0, 55.0)    # degrees North
lon_bounds = (-105.0, -35.0)  # degrees East (negative = West)

# =============================================================================
# 2) Determine integer‐index slice from thetao’s native grid
# =============================================================================
ds_thetao = xr.open_dataset(thetao_fn, decode_times=True, mask_and_scale=True)
lat_vals  = ds_thetao.latitude.values
lon_vals  = ds_thetao.longitude.values

ilat0 = np.searchsorted(lat_vals, lat_bounds[0], side="left")
ilat1 = np.searchsorted(lat_vals, lat_bounds[1], side="right")
ilon0 = np.searchsorted(lon_vals, lon_bounds[0], side="left")
ilon1 = np.searchsorted(lon_vals, lon_bounds[1], side="right")

print(f"Latitude index slice:  {ilat0} … {ilat1-1}  → "
      f"{lat_vals[ilat0]:.6f}° … {lat_vals[ilat1-1]:.6f}°")
print(f"Longitude index slice: {ilon0} … {ilon1-1}  → "
      f"{lon_vals[ilon0]:.6f}° … {lon_vals[ilon1-1]:.6f}°")

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
# 5) Subset SSH: pick first time & surface depth, then assign coords & time
# =============================================================================
ds_ssh_raw = xr.open_dataset(ssh_fn, decode_times=False, mask_and_scale=True)
zos = (
    ds_ssh_raw["sea_surface_height"]
      .isel(time=0, depth=0, drop=True)
      .rename("zos")
)
ds_ssh_sub = (
    zos.isel(latitude = slice(ilat0, ilat1),
             longitude= slice(ilon0, ilon1))
       .expand_dims(time=1)
       .assign_coords(time=ds_thetao_sub.time,
                      latitude=lat_grid,
                      longitude=lon_grid)
       .to_dataset()
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
ds_combined.to_netcdf(output_file)
print(f"Wrote merged file → {output_file}")

# =============================================================================
# 7) Verify grid‐spacing is uniform by computing diffs and plotting
# =============================================================================
lat = ds_combined.latitude.values
lon = ds_combined.longitude.values
dlat = np.diff(lat)
dlon = np.diff(lon)

print("Unique Δlatitude:", np.unique(np.round(dlat,8)))
print("Unique Δlongitude:", np.unique(np.round(dlon,8)))

#fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,6))
#ax1.scatter(np.arange(dlat.size), dlat, s=10, marker='.', color='tab:blue')
#ax1.set_ylabel('ΔLatitude (°)')
#ax1.set_title('Latitude Grid Spacing')
#ax1.grid(True)

#ax2.scatter(np.arange(dlon.size), dlon, s=10, marker='.', color='tab:green')
#ax2.set_xlabel('Index (between consecutive points)')
#ax2.set_ylabel('ΔLongitude (°)')
#ax2.set_title('Longitude Grid Spacing')
#ax2.grid(True)

#plt.tight_layout()
#plt.show()
