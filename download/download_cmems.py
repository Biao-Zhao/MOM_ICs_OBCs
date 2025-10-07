import copernicusmarine as cm
import os

# ========================
# User settings
# ========================

# Base output directory
base_dir = "/scratch/cimes/bz5265/MOM_ICs_OBCs/CMEMS_data/"

# Download area (Florida box by default)
lon_min, lon_max = -90, -30
lat_min, lat_max =  0,  50

# Depth range
depth_min, depth_max = 0, 7000

# ========================
# Script starts
# ========================

# Ask for date input
date_str = input("Enter date (yyyymmdd): ").strip()

# Build datetime strings
hour_str = "18"
datetime_str = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}T{hour_str}:00:00"
start_day   = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}T00:00:00"
end_day     = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}T23:00:00"

# Make output folder for this date
output_dir = os.path.join(base_dir, date_str)
os.makedirs(output_dir, exist_ok=True)

# --- Temperature (thetao) ---
thetao_file = os.path.join(output_dir, f"glo12_rg_6h-i_{date_str}-{hour_str}h_3D-thetao_hcst.nc")
if os.path.exists(thetao_file): os.remove(thetao_file)
cm.subset(
    dataset_id="cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i",
    variables=["thetao"],
#    minimum_longitude=lon_min, maximum_longitude=lon_max,
#    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=datetime_str, end_datetime=datetime_str,
    minimum_depth=depth_min, maximum_depth=depth_max,
    output_filename=thetao_file,
    force_download=True
)

# --- Salinity (so) ---
so_file = os.path.join(output_dir, f"glo12_rg_6h-i_{date_str}-{hour_str}h_3D-so_hcst.nc")
if os.path.exists(so_file): os.remove(so_file)
cm.subset(
    dataset_id="cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i",
    variables=["so"],
#    minimum_longitude=lon_min, maximum_longitude=lon_max,
#    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=datetime_str, end_datetime=datetime_str,
    minimum_depth=depth_min, maximum_depth=depth_max,
    output_filename=so_file,
    force_download=True
)

# --- Currents (uo + vo together) ---
cur_file = os.path.join(output_dir, f"glo12_rg_6h-i_{date_str}-{hour_str}h_3D-uovo_hcst.nc")
if os.path.exists(cur_file): os.remove(cur_file)
cm.subset(
    dataset_id="cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i",
    variables=["uo", "vo"],
#    minimum_longitude=lon_min, maximum_longitude=lon_max,
#    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=datetime_str, end_datetime=datetime_str,
    minimum_depth=depth_min, maximum_depth=depth_max,
    output_filename=cur_file,
    force_download=True
)

# --- Sea level (zos, full day hourly) ---
zos_file = os.path.join(output_dir, f"MOL_{date_str}.nc")
if os.path.exists(zos_file): os.remove(zos_file)
cm.subset(
    dataset_id="cmems_mod_glo_phy_anfc_merged-sl_PT1H-i",
    variables=["sea_surface_height","total_sea_level"],
#    minimum_longitude=lon_min, maximum_longitude=lon_max,
#    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=start_day, end_datetime=end_day,
    output_filename=zos_file,
    force_download=True
)

print(f"\n✅ All downloads complete for {date_str}")
print(f"Saved to: {output_dir}")
