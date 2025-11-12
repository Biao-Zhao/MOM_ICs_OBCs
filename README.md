<div align="center">

# 🌊 **prepare_MOM6_inputs.sh**
### Automated CMEMS GLORYS Data Processing for MOM6
**Author:** Dr. Biao Zhao · **Date:** 2025.10.04  
---

</div>

> ⚙️ *Automatically downloads and processes CMEMS GLORYS reanalysis data to generate initial and open boundary conditions for the MOM6 ocean model.*

---

## 🧭 **Overview**

`prepare_MOM6_inputs.sh` provides a **one-stop workflow** to prepare MOM6 input data:

1. ⬇️ Download CMEMS GLORYS data  
2. 🧊 Generate MOM6 **Initial Conditions (IC)**  
3. 🌐 Generate MOM6 **Open Boundary Conditions (OBC)**  

Each step can be run separately or all at once.

---

## 🚀 **Quick Start**

```bash
# Grant permission
chmod +x prepare_MOM6_inputs.sh

# Run the script
./prepare_MOM6_inputs.sh
