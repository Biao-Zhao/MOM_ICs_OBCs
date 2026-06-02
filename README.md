# prepare_MOM6_inputs.sh
### Automated CMEMS GLORYS Data Processing for MOM6

**Author:** Biao Zhao  
**Date:** 2025-10-04  

---

## Overview

This script integrates the Python scripts developed by **Dr. Jing Chen** for processing MOM6 initial and boundary conditions. It automates the downloading and processing of **CMEMS GLORYS analysis data**, generating **initial** and **open boundary condition** files for regional MOM6.

Main functions:
1. Download CMEMS GLORYS data  
2. Create MOM6 Initial Condition (IC) files  
3. Create MOM6 Open Boundary Condition (OBC) files  

Each step can be executed separately or together.

---

## Quick Start
## Quick Start

Make the script executable:

```bash
chmod +x prepare_MOM6_inputs.sh
```

Run the script with the default settings defined inside the file:

```bash
./prepare_MOM6_inputs.sh
```

The script also supports optional runtime arguments:

```bash
./prepare_MOM6_inputs.sh START_DATE START_HOUR END_DATE MODE
```

For example:

```bash
./prepare_MOM6_inputs.sh 2022-11-24 00 2022-12-03 1
```

This is equivalent to setting:

```bash
START_DATE="2022-11-24"
START_HOUR="00"
END_DATE="2022-12-03"
MODE="1"
```

If an argument is not provided, the corresponding default value defined inside `prepare_MOM6_inputs.sh` will be used.

---

## Available Modes

```bash
MODE="1"    # download CMEMS GLORYS data only
MODE="2"    # generate MOM6 initial condition only
MODE="3"    # generate MOM6 open boundary conditions only
MODE="all"  # run all steps
```
---
Set parameters (time range, region, resolution, etc.) in the user-defined section of the script before running.

## Directory Structure
Below is the expected directory organization under `BASE_DIR`.  
This structure ensures that downloaded data, grid files, and processing scripts are properly located.

```
BASE_DIR/
├── CMEMS/ → Downloaded GLORYS data
│
├── ICs/<res>/ → MOM6 initial condition outputs
│
├── OBCs/<res>/ → MOM6 boundary condition outputs
│
├── grid/<res>/ → MOM6 grid and vertical coordinate files
│
└── scripts/
├── download/ → download_cmems_glorys.py
├── initial/ → write_MOM6_IC_<res>.py
└── boundary/ → merge_glorys_<res>.py, write_MOM6_OBC_<res>.py

```   
---

## Workflow

The process consists of three main steps that can be executed separately or together.
```
Step 1 → Download GLORYS data
Step 2 → Generate MOM6 initial condition
Step 3 → Generate MOM6 open boundary condition
```
---

## Example Outputs
```
ICs/C3200/MOM6_IC_2024092000_C3200.nc
OBCs/C3200/20240920/thetao_001.nc
OBCs/C3200/20240920/so_002.nc
OBCs/C3200/20240920/uv_003.nc
```
