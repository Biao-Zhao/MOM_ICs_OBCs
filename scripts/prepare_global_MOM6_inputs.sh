#!/bin/bash
########################################################################################################
# prepare_MOM6_inputs.sh, written by Dr Biao Zhao                                                      #
#                                  --------- 2025.10.04 -------                                        #
#                                                                                                      #
# Purpose:                                                                                             #
#   Download CMEMS ocean reanalysis data for preparing initial condition for MOM6 model.               #
#                                                                                                      #
#                                                                                                      #
# Usage:                                                                                               #
#   Just edit the settings below, then run:                                                            #
#       chmod +x prepare_global_MOM6_inputs.sh                                                         #
#   Dowload CMEMS data:                                                                                #
#       ./prepare_global_MOM6_inputs.sh 2022-11-28 12 2022-11-30 1                                     #
#   Generate one initial condition                                                                     #
#       ./prepare_global_MOM6_inputs.sh 2022-11-28 12 2022-11-28 2                                     #
#   Reconstruct geostrophic currents for an existing initial condition                                 #
#       ./prepare_global_MOM6_inputs.sh 2022-11-28 12 2022-11-28 4                                     #
#                                                                                                      #
########################################################################################################
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate regrid
export PYTHONUNBUFFERED=1
set -e

# ======== User defined parameters ========

# set 1: Perform download GLORYS data process. set 2: only make initial condition. 
# set 3: only make bounday conditions.         set 4: only reconstruct geostrophic currents.
# set "all": run all steps
MODE="3"

# Start and end date
START_DATE="2022-11-26"
END_DATE="2022-12-01"

# Specify START_HOUR will only generate single initial condition, otherwise comment this line out
START_HOUR="12"

# use 6-hourly data in UTC, could be changed according needs, for example ("00" "06" "12" "18") or ("00" "01" "02" "03" "04" ........) 
TIME_SLOTS=("00" "06" "12" "18")

#grid cases
res="C768"

#vertical levels
NK="75"

# Number of processes used to run Kara flooding across vertical levels
KARA_WORKERS=32

# Apply the depth-integrated Poisson (barotropic) correction after rebuilding
# the geostrophic currents. Set to "false" to keep the uncorrected currents.
APPLY_BAROTROPIC_CORRECTION="false"

# runscripts and work directory
BASE_DIR="/ncrc/home1/Biao.Zhao/grid_prep/MOM_ICs_OBCs"
WORK_DIR="/gpfs/f6/bil-coastal-gfdl/scratch/Biao.Zhao/MOM_ICs_OBCs"

# Paths of downloading and making initial condition scripts
DOWNLOAD_SCRIPT="${BASE_DIR}/scripts/download/download_cmems_glorys.py"
INITIAL_SCRIPT="${BASE_DIR}/scripts/initial/write_MOM6_IC.py"
SIS2_INITIAL_SCRIPT="${BASE_DIR}/scripts/initial/write_SIS2_IC.py"
GEO_RECONSTRUCTION_SCRIPT="${BASE_DIR}/scripts/geostrophic_adj/python/reconstruction_current.py"

# Re/Analysis data directory
GLORYS_DIR="${WORK_DIR}/CMEMS"
# vertical grid and horizontal superrid file of MOM6
VGRID_FILE="${WORK_DIR}/grid/${res}/vgrid_${NK}.nc"
HGRID_FILE="${WORK_DIR}/grid/${res}/ocean_hgrid.nc"

# Path of store generateed files 
IC_OUTPUT_DIR="${WORK_DIR}/ICs/${res}/NK${NK}"
OBC_OUTPUT_DIR="${WORK_DIR}/OBCs/${res}/NK${NK}"
REGRID_WEIGHT_DIR="${WORK_DIR}/regrid_weights/${res}"

# use python3, can be changed according the local enviroment 
EXE="python3"

if [ -n "$1" ]; then
    START_DATE="$1"
    echo "START_DATE = ${START_DATE}"
fi

if [ -n "$2" ]; then
    START_HOUR="$2"
    echo "START_HOUR = ${START_HOUR}"
fi

if [ -n "$3" ]; then
    END_DATE="$3"
    echo "END_DATE   = ${END_DATE}"
fi

if [ -n "$4" ]; then
    MODE="$4"
    echo "MODE   = ${MODE}"
fi


# ===================================== Step 1: download glorys data  =================================
if [[ "$MODE" == "1" || "$MODE" == "all" ]]; then
# Calculate the day after END_DATE (for while loop)
END_NEXT=$(date -I -d "$END_DATE + 1 day")

CURRENT_DATE="$START_DATE"
while [[ "$CURRENT_DATE" != "$END_NEXT" ]]; do
  for HOUR in "${TIME_SLOTS[@]}"; do
    echo "[INFO] Downloading data for ${CURRENT_DATE} ${HOUR} UTC..."

    DOWNLOAD_EXTRA_ARGS=()
    # Sea ice is a daily mean, so download it once with the first time slot
    # on the initialization date. Its source time does not use START_HOUR.
    if [[ "$CURRENT_DATE" == "$START_DATE" && "$HOUR" == "${TIME_SLOTS[0]}" ]]; then
      DOWNLOAD_EXTRA_ARGS+=(--download-sea-ice)
    fi

    if [[ -n "${MIN_LON:-}" && -n "${MAX_LON:-}" && -n "${MIN_LAT:-}" && -n "${MAX_LAT:-}" ]]; then
      echo "Regional download"
      ${EXE} "$DOWNLOAD_SCRIPT" --outdir "$GLORYS_DIR" --date "$CURRENT_DATE" --hour "$HOUR" "${DOWNLOAD_EXTRA_ARGS[@]}" --min-lon "$MIN_LON" --max-lon "$MAX_LON" --min-lat "$MIN_LAT" --max-lat "$MAX_LAT"
    else
      echo "Global download"
      ${EXE} "$DOWNLOAD_SCRIPT" --outdir "$GLORYS_DIR" --date "$CURRENT_DATE" --hour "$HOUR" "${DOWNLOAD_EXTRA_ARGS[@]}"
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

echo "[INFO] Writing YAML configs for START_DATE (${START_DATE})"
DATE_COMPACT="${START_DATE//-/}"   # YYYYMMDD

if [[ -n "$START_HOUR" ]]; then
   echo "START_HOUR is set to ${START_HOUR}, only this hour will be processed."
   HOURS_TO_RUN=("$START_HOUR")
else
   echo "Generating initial conditions for TIME_SLOTS: ${TIME_SLOTS[*]}"
   HOURS_TO_RUN=("${TIME_SLOTS[@]}")
fi

for HOUR in "${HOURS_TO_RUN[@]}"; do
  # Paths for this hour (3D fields include the hour; SSH is daily file)
  THETAO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-thetao_hcst.nc"
  SO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-so_hcst.nc"
  UOVO_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_6h-i_${DATE_COMPACT}-${HOUR}h_3D-uovo_hcst.nc"
  SSH_PATH="${GLORYS_DIR}/${DATE_COMPACT}/MOL_${DATE_COMPACT}.nc"

  # IC output file
  IC_File="${IC_OUTPUT_DIR}/MOM6_IC_${DATE_COMPACT}${HOUR}_${res}.nc"
  mkdir -p ${WORK_DIR}/scripts/initial
  YAML="${WORK_DIR}/scripts/initial/glorys_IC_${res}.yaml"

  echo "[INFO] Writing ${YAML}"
  cat > "${YAML}" <<EOF
glorys_temperature: ${THETAO_PATH}
glorys_salinity: ${SO_PATH}
glorys_sea_surface_height: ${SSH_PATH}
glorys_zonal_velocity: ${UOVO_PATH}
glorys_meridional_velocity: ${UOVO_PATH}
resolution: ${res}
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
weight_dir: ${REGRID_WEIGHT_DIR}
reuse_weights: True
kara_workers: ${KARA_WORKERS}

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

# SIS2 reads a two-dimensional concentration/mass file for a new run. The
# daily-mean sea-ice source does not depend on START_HOUR, so create it once.
SEAICE_PATH="${GLORYS_DIR}/${DATE_COMPACT}/glo12_rg_1d-m_${DATE_COMPACT}-${DATE_COMPACT}_2D-ice.nc"
SIS2_IC_FILE="${IC_OUTPUT_DIR}/SIS2_IC_${DATE_COMPACT}_${res}.nc"
SIS2_REGION_ARGS=()

if [[ -n "${MIN_LON:-}" && -n "${MAX_LON:-}" && -n "${MIN_LAT:-}" && -n "${MAX_LAT:-}" ]]; then
  SIS2_REGION_ARGS+=(--min-lon "${MIN_LON}" --max-lon "${MAX_LON}")
  SIS2_REGION_ARGS+=(--min-lat "${MIN_LAT}" --max-lat "${MAX_LAT}")
elif [[ -n "${MIN_LON:-}" || -n "${MAX_LON:-}" || -n "${MIN_LAT:-}" || -n "${MAX_LAT:-}" ]]; then
  echo "ERROR: Set all four regional bounds or leave all four unset for global input." >&2
  exit 1
fi

${EXE} "${SIS2_INITIAL_SCRIPT}" \
  --input "${SEAICE_PATH}" \
  --grid "${HGRID_FILE}" \
  --output "${SIS2_IC_FILE}" \
  --weight-dir "${REGRID_WEIGHT_DIR}" \
  --resolution "${res}" \
  "${SIS2_REGION_ARGS[@]}" \
  --reuse-weights

echo "[INFO] Step 2 finished successfully."

else
  echo "[INFO] Skipping step 2 (Making initial condition)."
fi



# ===================================== Step 3: Making open boundary condition  =================================
if [[ "$MODE" == "3" || "$MODE" == "all" ]]; then

echo "[INFO] Skipping Step 3, global configuration doesn't need open boundary condition"

fi


# ================= Step 4: Reconstruct geostrophic currents =================

if [[ "${MODE}" == "4" ]]; then

    DATE_COMPACT="${START_DATE//-/}"

    if [[ -z "${START_HOUR}" ]]; then
        echo "[ERROR] START_HOUR must be specified for mode=4."
        exit 1
    fi

    IC_FILE="${IC_OUTPUT_DIR}/MOM6_IC_${DATE_COMPACT}${START_HOUR}_${res}.nc"

    IC_GEO_FILE="${IC_OUTPUT_DIR}/MOM6_IC_${DATE_COMPACT}${START_HOUR}_${res}_geocurrents.nc"

    GEO_PLOT_DIR="${IC_OUTPUT_DIR}/diagnostics_${DATE_COMPACT}${START_HOUR}"

    echo "[INFO] Running mode 4: geostrophic-current reconstruction"
    echo "[INFO] Input:  ${IC_FILE}"
    echo "[INFO] Output: ${IC_GEO_FILE}"
    echo "[INFO] Plots:  ${GEO_PLOT_DIR}"
    echo "[INFO] Apply barotropic correction: ${APPLY_BAROTROPIC_CORRECTION}"

    if [[ ! -f "${IC_FILE}" ]]; then
        echo "[ERROR] Initial-condition file does not exist:"
        echo "        ${IC_FILE}"
        exit 1
    fi

    if [[ -f "${IC_GEO_FILE}" ]]; then
        echo "[INFO] Geostrophic IC already exists, skipping:"
        echo "        ${IC_GEO_FILE}"
    else
    GEO_CORRECTION_ARG=""
    if [[ "${APPLY_BAROTROPIC_CORRECTION}" == "false" ]]; then
        GEO_CORRECTION_ARG="--skip-barotropic-correction"
    fi
    ${EXE} -u "${GEO_RECONSTRUCTION_SCRIPT}" \
            --input "${IC_FILE}" \
            --output "${IC_GEO_FILE}" \
            --grid-dir "${WORK_DIR}/grid/${res}" \
            --plot-dir "${GEO_PLOT_DIR}" \
            --plot-format svg \
            ${GEO_CORRECTION_ARG}

    echo "[INFO] Step 4 finished successfully."
    fi

else
    echo "[INFO] Skipping Step 4 (geostrophic reconstruction)."
fi



echo "[INFO] All done"
