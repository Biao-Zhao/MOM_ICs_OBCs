#!/bin/bash
#
########################################################################################################
# prepare_MOM6_inputs.sh, written by Dr Biao Zhao                                                      #
#                                  --------- 2025.10.04 -------                                        #
#                                                                                                      #
# Purpose:                                                                                             #
#   Download CMEMS ocean reanalysis data for preparing initial and open boundary conditions            #
#   for MOM6 model.                                                                                    #
#                                                                                                      #
#                                                                                                      #
# Usage:                                                                                               #
#   Just edit the settings below, then run:                                                            #
#       chmod +x prepare_MOM6_inputs.sh                                                                #
#       ./prepare_MOM6_inputs.sh                                                                       #
#                                                                                                      #
# Note:                                                                                                #
#   - If region bounds are not set, global data will be downloaded by default.                         #
########################################################################################################


set -e

# ======== User defined parameters ========

# set 1: Perform download GLORYS data process. set 2: only make initial condition. set 3: only make bounday conditions.  set "all": run all steps
MODE="2"
echo "[INFO] Running mode: $MODE"

# Start and end date
START_DATE="2024-09-20"
END_DATE="2024-09-22"

# use 6-hourly data in UTC, could be changed according needs, for example ("00" "06" "12" "18")
TIME_SLOTS=("00")

# Spefify which boundaries will be processed, for example ("south" "north" "east" "west") 
SEGMENTS=("south" "north" "east")

#grid cases
res="C3200"

# Output directory for GLORYS and IC & OBC data
BASE_DIR="/scratch/cimes/bz5265/MOM_ICs_OBCs"
GLORYS_DIR="${BASE_DIR}/CMEMS"
IC_OUTPUT_DIR="${BASE_DIR}/ICs/${res}"
OBC_OUTPUT_DIR="${BASE_DIR}/OBCs/${res}"

# Paths of downloading, making initial $ open boundary condition scripts
DOWNLOAD_SCRIPT="${BASE_DIR}/scripts/download/download_cmems_glorys.py"
INITIAL_SCRIPT="${BASE_DIR}/scripts/initial/write_MOM6_IC_${res}.py"
BOUNDARY_MERGE_SCRIPT="${BASE_DIR}/scripts/boundary/merge_glorys_${res}.py"
BOUNDARY_MAKE_SCRIPT="${BASE_DIR}/scripts/boundary/write_MOM6_OBC_${res}.py"

# vertical grid and horizontal superrid file of MOM6
VGRID_FILE="${BASE_DIR}/grid/${res}/vgrid_75_2m.nc"
HGRID_FILE="${BASE_DIR}/grid/${res}/ocean_hgrid.nc"

# Optional region boundaries (comment out if you want global data)
MIN_LON=-105
MAX_LON=-35
MIN_LAT=10
MAX_LAT=55

# use python3, can be changed according the local enviroment 
EXE="python3"


# ===================================== Step 1: download glorys data  =================================
if [[ "$MODE" == "1" || "$MODE" == "all" ]]; then
# Calculate the day after END_DATE (for while loop)
END_NEXT=$(date -I -d "$END_DATE + 1 day")

CURRENT_DATE="$START_DATE"
while [[ "$CURRENT_DATE" != "$END_NEXT" ]]; do
  for HOUR in "${TIME_SLOTS[@]}"; do
    echo "[INFO] Downloading data for ${CURRENT_DATE} ${HOUR} UTC..."

    if [[ -n "${MIN_LON:-}" && -n "${MAX_LON:-}" && -n "${MIN_LAT:-}" && -n "${MAX_LAT:-}" ]]; then
      echo "Regional download"
      ${EXE} "$DOWNLOAD_SCRIPT"  --outdir "$GLORYS_DIR" --date "$CURRENT_DATE" --hour "$HOUR" --min-lon "$MIN_LON" --max-lon "$MAX_LON" --min-lat "$MIN_LAT" --max-lat "$MAX_LAT"
    else
      echo "Global download"
      ${EXE} "$DOWNLOAD_SCRIPT" --outdir "$GLORYS_DIR" --date "$CURRENT_DATE" --hour "$HOUR"
    fi
  done

  # Move to next day
  CURRENT_DATE=$(date -I -d "$CURRENT_DATE + 1 day")
done

echo "[INFO] Step 1: All downloads completed successfully"

else
  echo "[INFO] Skipping Step 1 (download)."
fi



# ===================================== Step 2: Making initial condition  =================================
if [[ "$MODE" == "2" || "$MODE" == "all" ]]; then

echo "[INFO] Writing YAML configs for START_DATE (${START_DATE}) at hours: ${TIME_SLOTS[*]}"
DATE_COMPACT="${START_DATE//-/}"   # YYYYMMDD

for HOUR in "${TIME_SLOTS[@]}"; do
  # Paths for this hour (3D fields include the hour; SSH is daily file)
  THETAO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-thetao_hcst.nc"
  SO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-so_hcst.nc"
  UOVO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-uovo_hcst.nc"
  SSH_PATH="${GLORYS_DIR}/${DATE_COMPACT}/MOL_${DATE_COMPACT}.nc"

  # IC output file
  IC_File="${IC_OUTPUT_DIR}/MOM6_IC_${DATE_COMPACT}${HOUR}_${res}.nc"

  YAML="${BASE_DIR}/scripts/initial/glorys_IC_${res}.yaml"

  echo "[INFO] Writing ${YAML}"
  cat > "${YAML}" <<EOF
glorys_temperature: ${THETAO_PATH}
glorys_salinity: ${SO_PATH}
glorys_sea_surface_height: ${SSH_PATH}
glorys_zonal_velocity: ${UOVO_PATH}
glorys_meridional_velocity: ${UOVO_PATH}
ssh_time: ${HOUR}
# Paths to model grid files
vgrid_file: ${VGRID_FILE}
grid_file: ${HGRID_FILE}
# define the area to cut out
min_lon: ${MIN_LON}
max_lon: ${MAX_LON}
min_lat: ${MIN_LAT}
max_lat: ${MAX_LAT}

# Output NetCDF file
output_file: ${IC_File}

# Whether to reuse existing regridding weights (if applicable)
reuse_weights: False

# Variable names inside the NetCDF files
variable_names:
  temperature: thetao
  salinity: so
  sea_surface_height: sea_surface_height
  zonal_velocity: uo
  meridional_velocity: vo
EOF

 ${EXE} ${INITIAL_SCRIPT} --config_file  ${YAML}

done

echo "[INFO] Step 2 finished successfully."

else
  echo "[INFO] Skipping pahse 2 (Making initial condition)."
fi



# ===================================== Step 3: Making open boundary condition  =================================
if [[ "$MODE" == "3" || "$MODE" == "all" ]]; then

for i in "${!SEGMENTS[@]}"; do SEGS+=($((i+1))); done

END_NEXT=$(date -I -d "$END_DATE + 1 day")

CURRENT_DATE="$START_DATE"

while [[ "$CURRENT_DATE" != "$END_NEXT" ]]; do
  for HOUR in "${TIME_SLOTS[@]}"; do

    echo "[INFO] Merging booundary data for ${CURRENT_DATE} ${HOUR} UTC..."

    DATE_COMPACT="${CURRENT_DATE//-/}"
    THETAO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-thetao_hcst.nc"
    SO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-so_hcst.nc"
    UOVO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-uovo_hcst.nc"
    SSH_PATH="${GLORYS_DIR}/${DATE_COMPACT}/MOL_${DATE_COMPACT}.nc"

    YAML_MERGE="${BASE_DIR}/scripts/boundary/config_merge_glorys.yaml"
  
    mkdir -p ${BASE_DIR}/scripts/temp/${DATE_COMPACT}/
    echo "[INFO] Writing ${YAML_MERGE}"

    cat > "${YAML_MERGE}" <<EOF
thetao_fn: ${THETAO_PATH}
so_fn: ${SO_PATH}
uovo_fn: ${UOVO_PATH}
ssh_fn: ${SSH_PATH}
merged_fn: ${BASE_DIR}/scripts/temp/${DATE_COMPACT}/merged_${DATE_COMPACT}-${HOUR}.nc
ssh_time: ${HOUR}
min_lon: ${MIN_LON}
max_lon: ${MAX_LON}
min_lat: ${MIN_LAT}
max_lat: ${MAX_LAT}
EOF

    # 1. merging the salinity, temperature, currents and sea surface elevation files into one single file on model domain
    ${EXE} $BOUNDARY_MERGE_SCRIPT  --config ${YAML_MERGE}

    YAML_REGRID="${BASE_DIR}/scripts/boundary/config_regrid_glorys.yaml"

    echo "[INFO] Writing ${YAML_REGRID}"
{
cat <<EOF
glorys_dir: ${BASE_DIR}/scripts/temp/${DATE_COMPACT}
output_dir: ${BASE_DIR}/scripts/temp/
hgrid: ${HGRID_FILE} 
ncrcat_years: true  # Set to false if you want to skip ncrcat_years
ncrcat_names:
  - 'thetao'
  - 'so'
  - 'zos'
  - 'uv'
variables:
  - 'thetao'
  - 'so'
  - 'zos'
  - 'uv'
segments:
EOF
  i=1
  for b in "${SEGMENTS[@]}"; do
    echo "  - id: ${i}"
    echo "    border: '${b}'"
    i=$((i+1))
  done
}> "${YAML_REGRID}"

  # 2. extract information at the boundaries specifieed in "config_regrid_glorys.yaml"
  ${EXE} ${BOUNDARY_MAKE_SCRIPT} --config ${YAML_REGRID} --year ${DATE_COMPACT:0:4} --month ${DATE_COMPACT:4:2} --day ${DATE_COMPACT:6:2} --hour ${HOUR}

  done
  # Move to next day
  CURRENT_DATE=$(date -I -d "$CURRENT_DATE + 1 day")
done

  # 3. Merging booundary data
  mkdir -p ${OBC_OUTPUT_DIR}/${START_DATE:0:4}${START_DATE:5:2}${START_DATE:8:2} 
  VARS=(thetao so zos uv)
  for var in "${VARS[@]}"; do
    for seg in "${SEGS[@]}"; do
      segnum=$(printf "%03d" "$seg")
      mapfile -t files < <(find "${BASE_DIR}/scripts/temp/" -maxdepth 1 -type f -name "${var}_${segnum}_*.nc" | sort -V)
      if (( ${#files[@]} == 0 )); then
        echo "[ERROR] No input files for ${var}_${segnum} in ${BASE_DIR}/scripts/temp/."
        echo "[ERROR] Expected pattern: ${var}_${segnum}_YYYYMMDD-HH.nc"
        exit 1
      fi
      out="${OBC_OUTPUT_DIR}/${START_DATE:0:4}${START_DATE:5:2}${START_DATE:8:2}/${var}_${segnum}.nc"
      printf '    %s\n' "${files[@]}"
      echo "[INFO] ncrcat ${#files[@]} files → ${out}"
      ncrcat "${files[@]}" "${out}"
    done
  done

  rm -r ${BASE_DIR}/scripts/temp/ 
echo "[INFO] Step 3 finished successfully."

fi

echo "[INFO] All done"
