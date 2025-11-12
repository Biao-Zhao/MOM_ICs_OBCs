# prepare_MOM6_inputs.sh
### Automated CMEMS GLORYS Data Processing for MOM6

**Author:** Dr. Biao Zhao  
**Date:** 2025-10-04  

---

## Overview

This script automates the preparation of **initial and open boundary conditions** for the MOM6 ocean model.  
It downloads and processes **CMEMS GLORYS reanalysis data**, generating MOM6-ready NetCDF input files.

Main functions:
1. Download CMEMS GLORYS data  
2. Create MOM6 Initial Condition (IC) files  
3. Create MOM6 Open Boundary Condition (OBC) files  

Each step can be executed separately or together.

---
## Directory Structure
BASE_DIR/
├── CMEMS/                  →  Downloaded GLORYS data
├── ICs/<res>/              →  MOM6 initial condition outputs
├── OBCs/<res>/             →  MOM6 boundary condition outputs
├── grid/<res>/             →  MOM6 grid and vgrid files
└── scripts/
    ├── download/
    ├── initial/
    └── boundary/
---
## Quick Start

```bash
# Make the script executable
chmod +x prepare_MOM6_inputs.sh

# Run the script
./prepare_MOM6_inputs.sh

