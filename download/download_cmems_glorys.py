#!/usr/bin/env python3
"""
download_cmems_glorys.py, written by Drs. Biao Zhao and Jing Chen
                                     -------- 2025.10.04 --------
Purpose:
  Download CMEMS ocean reanalysis data for a given date/hour.
  All parameters (output directory, date, hour, region) must be passed by the controller script.

Inputs (required):
  --outdir  : Base output directory
  --date    : Date (YYYY-MM-DD)
  --hour    : Hour (00, 06, 12, 18)

Optional:
  --min-lon, --max-lon, --min-lat, --max-lat  : If all given → regional subset; else → global

Behavior:
  - Downloads 3D fields (thetao, so, uo, vo) at given 6-hour timestamp
  - Downloads full-day hourly sea level only at 00 UTC
  - Output: <outdir>/YYYYMMDD/
"""

import argparse
import os
import sys
from datetime import datetime

import copernicusmarine as cm

# Dataset IDs
DATASET_THETAO = "cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i"
DATASET_SO     = "cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i"
DATASET_CUR    = "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i"
DATASET_ZOS    = "cmems_mod_glo_phy_anfc_merged-sl_PT1H-i"

# Depth range
DEPTH_MIN, DEPTH_MAX = 0, 7000


def parse_args():
    p = argparse.ArgumentParser(description="Download CMEMS data for given date/hour (regional or global).")
    p.add_argument("--outdir", required=True, help="Base output directory")
    p.add_argument("--date", required=True, help="Date (YYYY-MM-DD)")
    p.add_argument("--hour", required=True, help="Hour (00, 06, 12, 18)")
    p.add_argument("--min-lon", type=float)
    p.add_argument("--max-lon", type=float)
    p.add_argument("--min-lat", type=float)
    p.add_argument("--max-lat", type=float)
    args = p.parse_args()

    # Validate date
    try:
        datetime.strptime(args.date, "%Y-%m-%d")
    except ValueError:
        sys.exit("ERROR: --date must be in YYYY-MM-DD format.")

    # Validate hour
    if args.hour not in {"00", "06", "12", "18"}:
        sys.exit('ERROR: --hour must be one of "00", "06", "12", "18".')

    # Build region arguments (only if all provided)
    region_kwargs = {}
    if all(v is not None for v in (args.min_lon, args.max_lon, args.min_lat, args.max_lat)):
        region_kwargs = dict(
            minimum_longitude=args.min_lon,
            maximum_longitude=args.max_lon,
            minimum_latitude=args.min_lat,
            maximum_latitude=args.max_lat,
        )

    return argparse.Namespace(
        outdir=args.outdir,
        date=args.date,
        hour=args.hour,
        region_kwargs=region_kwargs,
    )


def yyyymmdd(date_str):
    return date_str.replace("-", "")


def main():
    args = parse_args()

    date_compact = yyyymmdd(args.date)
    timestamp = f"{args.date}T{args.hour}:00:00"
    day_start, day_end = f"{args.date}T00:00:00", f"{args.date}T23:00:00"

    out_dir = os.path.join(args.outdir, date_compact)
    os.makedirs(out_dir, exist_ok=True)

    # --- Temperature ---
    f_thetao = os.path.join(out_dir, f"glo12_rg_6h-i_{date_compact}-{args.hour}h_3D-thetao_hcst.nc")
    if os.path.exists(f_thetao):
        os.remove(f_thetao)
    cm.subset(
        dataset_id=DATASET_THETAO,
        variables=["thetao"],
        start_datetime=timestamp,
        end_datetime=timestamp,
        minimum_depth=DEPTH_MIN,
        maximum_depth=DEPTH_MAX,
        output_filename=f_thetao,
        force_download=True,
        **args.region_kwargs,
    )

    # --- Salinity ---
    f_so = os.path.join(out_dir, f"glo12_rg_6h-i_{date_compact}-{args.hour}h_3D-so_hcst.nc")
    if os.path.exists(f_so):
        os.remove(f_so)
    cm.subset(
        dataset_id=DATASET_SO,
        variables=["so"],
        start_datetime=timestamp,
        end_datetime=timestamp,
        minimum_depth=DEPTH_MIN,
        maximum_depth=DEPTH_MAX,
        output_filename=f_so,
        force_download=True,
        **args.region_kwargs,
    )

    # --- Currents ---
    f_cur = os.path.join(out_dir, f"glo12_rg_6h-i_{date_compact}-{args.hour}h_3D-uovo_hcst.nc")
    if os.path.exists(f_cur):
        os.remove(f_cur)
    cm.subset(
        dataset_id=DATASET_CUR,
        variables=["uo", "vo"],
        start_datetime=timestamp,
        end_datetime=timestamp,
        minimum_depth=DEPTH_MIN,
        maximum_depth=DEPTH_MAX,
        output_filename=f_cur,
        force_download=True,
        **args.region_kwargs,
    )

    # --- Sea level (only once per day) ---
    if args.hour == "00":
        f_zos = os.path.join(out_dir, f"MOL_{date_compact}.nc")
        if os.path.exists(f_zos):
            os.remove(f_zos)
        print(f"[INFO] Downloading daily sea level for {args.date} ...")
        cm.subset(
            dataset_id=DATASET_ZOS,
            variables=["sea_surface_height", "total_sea_level"],
            start_datetime=day_start,
            end_datetime=day_end,
            output_filename=f_zos,
            force_download=True,
            **args.region_kwargs,
        )
    else:
        print(f"[INFO] Skipping sea level for {args.date} (already handled at 00 UTC).")

    # --- Summary ---
    print(f"\n✅ Completed: {args.date} {args.hour} UTC")
    print(f"   Output: {out_dir}")
    if args.region_kwargs:
        print(f"   Region: lon[{args.region_kwargs['minimum_longitude']}, {args.region_kwargs['maximum_longitude']}], "
              f"lat[{args.region_kwargs['minimum_latitude']}, {args.region_kwargs['maximum_latitude']}]")
    else:
        print("   Region: global")


if __name__ == "__main__":
    main()

