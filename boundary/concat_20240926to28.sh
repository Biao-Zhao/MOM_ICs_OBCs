#!/bin/bash
#
# Author: Jing Chen
# Date: July 7, 2025
# Description: Concatenate NetCDF files with the same variable prefix using ncrcat.
# group NetCDF files with the same variable prefix (e.g., so_001, thetao_002, etc.) together into a single concatenated file.
#  chmod +x concat_by_20240926to28.sh 
# ./concat_by_20240926to28.sh 
# mkdir -p ./no2024
# ls *.nc | grep -v 2024 | xargs -I {} cp {} ./no2024/
# Load NCO module (adjust version as needed)

# List of variable prefixes to process
prefixes=(
  so_001
  so_002
  so_003
  thetao_001
  thetao_002
  thetao_003
  uv_001
  uv_002
  uv_003
  zos_001
  zos_002
  zos_003
)

# Loop to concatenate files
for prefix in "${prefixes[@]}"; do
  echo "Processing ${prefix}..."
  ncrcat ${prefix}_*.nc ./${prefix}.nc
done

echo "All done."

