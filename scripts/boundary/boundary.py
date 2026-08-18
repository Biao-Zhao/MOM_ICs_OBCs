#!/usr/bin/env python3
"""
Standalone MOM6 boundary utilities with persistent, resolution-specific
XESMF weights.

This file can be copied directly over boundary.py on the HPC. It contains
the functions and Segment methods used by write_MOM6_OBC.py and does not
import another boundary module.
"""

from os import path
import warnings

import numpy as np
import xarray
import xesmf


def check_angle_range(angle):
    """Verify that the model-grid angle is expressed in radians."""
    amax = float(angle.max())
    amin = float(angle.min())
    if amax > (2 * np.pi) or amin < (-2 * np.pi):
        raise ValueError(
            f"Grid angle ranges from [{amin}, {amax}]. "
            "Expected values in [-2pi, 2pi]. Are the units correct?"
        )


def rotate_uv(uearth, vearth, angle_earth_to_model_rad):
    """Rotate earth-relative velocity to model-relative velocity."""
    urot = (
        np.cos(angle_earth_to_model_rad) * uearth
        + np.sin(angle_earth_to_model_rad) * vearth
    )
    vrot = (
        -np.sin(angle_earth_to_model_rad) * uearth
        + np.cos(angle_earth_to_model_rad) * vearth
    )
    return urot, vrot


def fill_missing(arr, xdim="locations", zdim="z", fill="b"):
    """Fill missing boundary values horizontally, then vertically."""
    if fill == "f":
        filled = arr.ffill(dim=xdim, limit=None)
    elif fill == "b":
        filled = arr.bfill(dim=xdim, limit=None)
    else:
        raise ValueError("fill must be either 'f' or 'b'")

    if zdim is not None:
        filled = filled.ffill(dim=zdim, limit=None).fillna(0)
    return filled


def flood_missing(arr, **kwargs):
    """Flood source-grid land values with HCtFlood."""
    from HCtFlood import kara as hct

    flooded = hct.flood_kara(arr, **kwargs)
    if arr.ndim <= 3 and "zdim" not in kwargs:
        if "z" in arr.dims and len(arr.z) > 1:
            warnings.warn(
                "flood_kara used the default z dimension; "
                "not dropping it."
            )
        else:
            flooded = flooded.isel(z=0).drop_vars("z")
    return flooded


def reuse_regrid(*args, **kwargs):
    """Create a missing Regridder weight file or reuse an existing one."""
    filename = kwargs.pop("filename", None)
    reuse_weights = kwargs.pop("reuse_weights", False)

    if reuse_weights and path.isfile(filename):
        return xesmf.Regridder(
            *args,
            reuse_weights=True,
            filename=filename,
            **kwargs,
        )

    regrid = xesmf.Regridder(*args, **kwargs)
    if reuse_weights:
        regrid.to_netcdf(filename)
    return regrid


def z_to_dz(ds, max_depth=6500.0):
    """Convert source layer-center depths to layer thicknesses."""
    zi = 0.5 * (np.roll(ds["z"], shift=-1) + ds["z"])
    zi[-1] = max_depth
    dz = zi - np.roll(zi, shift=1)
    dz[0] = zi[0]

    nt = len(ds["time"])
    nz = len(ds["z"])
    nx = len(ds["locations"])
    dz = np.tile(dz.data[np.newaxis, :, np.newaxis], (nt, 1, nx))

    da_dz = xarray.DataArray(
        dz,
        coords=[
            ("time", ds["time"].data),
            ("z", ds["z"].data),
            ("locations", ds["locations"].data),
        ],
    )
    for var in ["time", "z", "locations"]:
        da_dz[var].attrs = ds[var].attrs
    return da_dz


def cached_regrid(*args, filename, **kwargs):
    """Create a missing weight file or reuse an existing one."""
    exists = path.isfile(filename)
    action = "Reusing" if exists else "Generating"
    print(f"{action} regridding weights: {filename}", flush=True)

    return reuse_regrid(
        *args,
        filename=filename,
        reuse_weights=True,
        **kwargs,
    )


class Segment:
    """A MOM6 boundary segment with persistent weight-file support."""

    def __init__(
        self,
        num,
        border,
        hgrid,
        in_degrees=False,
        output_dir=".",
        regrid_dir=None,
        resolution="unknown",
    ):
        self.num = num
        self.border = border
        self.resolution = str(resolution)
        self.hgrid = hgrid.copy(deep=True)

        angle_units = hgrid["angle_dx"].attrs.get("units")
        if angle_units == "degrees" or in_degrees:
            print("Converting grid angle from degrees to radians")
            self.hgrid["angle_dx"] = np.radians(
                self.hgrid["angle_dx"]
            )
        check_angle_range(self.hgrid["angle_dx"])

        self.segstr = f"segment_{self.num:03d}"
        self.output_dir = output_dir
        self.regrid_dir = (
            output_dir if regrid_dir is None else regrid_dir
        )

    @property
    def coords(self):
        """Return longitude, latitude, and angle along this boundary."""
        if self.border == "south":
            return xarray.Dataset(
                {
                    "lon": self.hgrid["x"].isel(nyp=0),
                    "lat": self.hgrid["y"].isel(nyp=0),
                    "angle": self.hgrid["angle_dx"].isel(nyp=0),
                }
            )
        if self.border == "north":
            return xarray.Dataset(
                {
                    "lon": self.hgrid["x"].isel(nyp=-1),
                    "lat": self.hgrid["y"].isel(nyp=-1),
                    "angle": self.hgrid["angle_dx"].isel(nyp=-1),
                }
            )
        if self.border == "west":
            return xarray.Dataset(
                {
                    "lon": self.hgrid["x"].isel(nxp=0),
                    "lat": self.hgrid["y"].isel(nxp=0),
                    "angle": self.hgrid["angle_dx"].isel(nxp=0),
                }
            )
        if self.border == "east":
            return xarray.Dataset(
                {
                    "lon": self.hgrid["x"].isel(nxp=-1),
                    "lat": self.hgrid["y"].isel(nxp=-1),
                    "angle": self.hgrid["angle_dx"].isel(nxp=-1),
                }
            )
        raise ValueError(
            "border must be south, north, west, or east"
        )

    @property
    def nx(self):
        if self.border in ["south", "north"]:
            return len(self.coords["lon"])
        return 1

    @property
    def ny(self):
        if self.border in ["west", "east"]:
            return len(self.coords["lat"])
        return 1

    def to_netcdf(
        self,
        ds,
        varnames,
        suffix=None,
        additional_encoding=None,
    ):
        """Write one segment variable group to NetCDF4."""
        for var in ds:
            ds[var].encoding["_FillValue"] = 1.0e20

        if suffix is None:
            fname = f"{varnames}_{self.num:03d}.nc"
        else:
            fname = f"{varnames}_{self.num:03d}_{suffix}.nc"

        ds[f"lon_{self.segstr}"].encoding["dtype"] = "float64"
        ds[f"lat_{self.segstr}"].encoding["dtype"] = "float64"

        if (
            "calendar" not in ds["time"].attrs
            and "modulo" not in ds["time"].attrs
        ):
            ds["time"].encoding["calendar"] = "gregorian"
            ds["time"].encoding["dtype"] = "float64"
            ds["time"].encoding["_FillValue"] = 1.0e20

        allowed_keys = {
            "_FillValue",
            "dtype",
            "zlib",
            "complevel",
            "chunksizes",
            "scale_factor",
            "add_offset",
        }
        enc = {}
        for name in list(ds.data_vars) + ["time"]:
            raw = ds[name].encoding or {}
            enc[name] = {
                key: raw[key]
                for key in raw
                if key in allowed_keys
            }

        if additional_encoding is not None:
            enc.update(additional_encoding)

        ds.to_netcdf(
            path.join(self.output_dir, fname),
            mode="w",
            format="NETCDF4",
            engine="netcdf4",
            encoding=enc,
            unlimited_dims="time",
        )

    def expand_dims(self, ds):
        """Add the length-one dimension normal to the boundary."""
        offset = 0 if ("z" in ds.coords or "constituent" in ds.dims) else 1
        if self.border in ["south", "north"]:
            return ds.expand_dims(
                f"ny_{self.segstr}",
                2 - offset,
            )
        return ds.expand_dims(
            f"nx_{self.segstr}",
            3 - offset,
        )

    def rename_dims(self, ds):
        """Rename coordinates and dimensions to MOM6 segment names."""
        ds = ds.rename(
            {
                "lon": f"lon_{self.segstr}",
                "lat": f"lat_{self.segstr}",
            }
        )
        if "z" in ds.coords:
            ds = ds.rename({"z": f"nz_{self.segstr}"})
        if self.border in ["south", "north"]:
            return ds.rename(
                {"locations": f"nx_{self.segstr}"}
            )
        return ds.rename(
            {"locations": f"ny_{self.segstr}"}
        )

    def weight_filename(self, suffix):
        """Return a resolution-specific weight filename."""
        return path.join(
            self.regrid_dir,
            f"regrid_{self.resolution}_{self.segstr}_{suffix}.nc",
        )

    def regrid_velocity(
        self,
        usource,
        vsource,
        method="nearest_s2d",
        periodic=False,
        write=True,
        flood=False,
        fill="b",
        xdim="lon",
        ydim="lat",
        zdim="z",
        rotate=True,
        time_attrs=None,
        time_encoding=None,
        **kwargs,
    ):
        """Regrid earth-relative velocity using persistent U/V weights."""
        if flood:
            usource = flood_missing(
                usource,
                xdim=xdim,
                ydim=ydim,
                zdim=zdim,
            ).load()
            vsource = flood_missing(
                vsource,
                xdim=xdim,
                ydim=ydim,
                zdim=zdim,
            ).load()

        uregrid = cached_regrid(
            usource,
            self.coords,
            method=method,
            locstream_out=True,
            periodic=periodic,
            filename=self.weight_filename("u"),
        )
        vregrid = cached_regrid(
            vsource,
            self.coords,
            method=method,
            locstream_out=True,
            periodic=periodic,
            filename=self.weight_filename("v"),
        )

        udest = uregrid(usource)
        vdest = vregrid(vsource)

        if isinstance(udest, xarray.Dataset):
            udest = udest.to_array().squeeze()
        if isinstance(vdest, xarray.Dataset):
            vdest = vdest.to_array().squeeze()

        if rotate:
            if self.border in ["south", "north"]:
                udest = udest.rename({"nxp": "locations"})
                vdest = vdest.rename({"nxp": "locations"})
                angle = self.coords["angle"].rename(
                    {"nxp": "locations"}
                )
            elif self.border in ["west", "east"]:
                udest = udest.rename({"nyp": "locations"})
                vdest = vdest.rename({"nyp": "locations"})
                angle = self.coords["angle"].rename(
                    {"nyp": "locations"}
                )
            udest, vdest = rotate_uv(udest, vdest, angle)

        ds_uv = xarray.Dataset(
            {
                f"u_{self.segstr}": udest,
                f"v_{self.segstr}": vdest,
            }
        )

        ds_uv = fill_missing(ds_uv, fill=fill)
        ds_uv = ds_uv.transpose("time", "z", "locations")

        dz = z_to_dz(ds_uv)
        ds_uv[f"dz_u_{self.segstr}"] = dz
        ds_uv[f"dz_v_{self.segstr}"] = dz
        ds_uv["z"] = np.arange(len(ds_uv["z"]))

        ds_uv = self.expand_dims(ds_uv)

        if "lon" not in ds_uv.variables:
            ds_uv.update(
                {"lon": ("locations", self.coords["lon"].values)}
            )
        if "lat" not in ds_uv.variables:
            ds_uv.update(
                {"lat": ("locations", self.coords["lat"].values)}
            )

        ds_uv = self.rename_dims(ds_uv)

        if time_attrs:
            ds_uv["time"].attrs = time_attrs
        if time_encoding:
            ds_uv["time"].encoding = time_encoding

        if write:
            self.to_netcdf(ds_uv, "uv", **kwargs)

        return ds_uv

    def regrid_tracer(
        self,
        tsource,
        method="nearest_s2d",
        periodic=False,
        write=True,
        flood=False,
        fill="b",
        xdim="lon",
        ydim="lat",
        zdim="z",
        regrid_suffix="t",
        source_var=None,
        time_attrs=None,
        time_encoding=None,
        **kwargs,
    ):
        """Regrid a tracer using a persistent segment weight file."""
        if source_var is None:
            name = tsource.name
            if flood:
                tsource = flood_missing(
                    tsource,
                    xdim=xdim,
                    ydim=ydim,
                    zdim=zdim,
                ).load()
        else:
            name = source_var
            if flood:
                tsource[name] = flood_missing(
                    tsource[name],
                    xdim=xdim,
                    ydim=ydim,
                    zdim=zdim,
                ).load()

        regrid = cached_regrid(
            tsource,
            self.coords,
            method=method,
            locstream_out=True,
            periodic=periodic,
            filename=self.weight_filename(regrid_suffix),
        )
        tdest = regrid(tsource)

        if not isinstance(tdest, xarray.Dataset):
            tdest.name = name
            tdest = tdest.to_dataset()

        xname = list(tdest.dims)[-1]
        tdest = tdest.rename({xname: "locations"})

        if "z" in tsource.coords:
            tdest = fill_missing(tdest, fill=fill)
            tdest = tdest.transpose("time", "z", "locations")
            dz = z_to_dz(tdest)
            tdest[f"dz_{name}_{self.segstr}"] = dz
            tdest["z"] = np.arange(len(tdest["z"]))
        else:
            tdest = fill_missing(tdest, zdim=None, fill=fill)
            tdest = tdest.transpose("time", "locations")

        tdest = self.expand_dims(tdest)
        tdest["lon"] = (("locations",), self.coords["lon"].data)
        tdest["lat"] = (("locations",), self.coords["lat"].data)

        tdest = self.rename_dims(tdest)
        tdest = tdest.rename({name: f"{name}_{self.segstr}"})

        if time_attrs:
            tdest["time"].attrs = time_attrs
        if time_encoding:
            tdest["time"].encoding = time_encoding

        if write:
            self.to_netcdf(tdest, name, **kwargs)

        return tdest

