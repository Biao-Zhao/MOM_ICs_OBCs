#!/usr/bin/env python3
"""Reconstruct geostrophic currents in a MOM6 initial-condition file.

Requirements
------------
* Python 3.10 or newer.
* NumPy, SciPy, and netCDF4.
* Matplotlib, unless ``--no-plots`` is used.

For example, install the dependencies in an active Conda environment with::

    conda install -c conda-forge numpy scipy netcdf4 matplotlib

Required input files
--------------------
``--input`` must point to a MOM6 initial-condition NetCDF file containing
``temp``, ``salt``, ``ssh``, ``u``, ``v``, ``zl``, ``xh``, ``yh``, ``xq``,
and ``yq``.

``--grid-dir`` must point to a directory containing:

* ``ocean_hgrid.nc``
* ``ocean_mask.nc``
* ``topog.nc``

Usage example
-------------
::

    python -u reconstruction_current.py \\
      --input /path/to/MOM6_IC_2022112712_C3200.nc \\
      --output /path/to/MOM6_IC_2022112712_C3200_geocurrents_python.nc \\
      --grid-dir /path/to/grid/C3200 \\
      --work-dir "$TMPDIR" \\
      --plot-dir /path/to/diagnostics \\
      --reference /path/to/MOM6_IC_2022112712_C3200_geocurrents.nc \\
      --overwrite

Important options
-----------------
``--work-dir DIR``
    Perform the repeated writes to a temporary NetCDF file in ``DIR``, then
    safely publish the completed file to ``--output``.  On an HPC system,
    use a node-local SSD or ``$TMPDIR`` when it has enough free space.  Omit
    this option if local space is insufficient.

``--plot-dir DIR``
    Save diagnostic PNG files in ``DIR``.  If omitted, the default directory
    is ``<output_stem>_diagnostics`` beside the output file.

``--plot-format {png,svg}``
    Select the diagnostic figure format.  Use ``svg`` on HPC systems where
    the Matplotlib/FreeType PNG renderer reports ``FT_Render_Glyph`` or
    ``raster overflow``.  SVG text bypasses the failing glyph rasterizer.

``--no-plots``
    Disable diagnostic plotting.  Matplotlib is not required in this mode.

``--skip-barotropic-correction``
    Write the reconstructed geostrophic currents without solving the
    depth-integrated Poisson problem or applying its barotropic correction.

``--reference FILE``
    Compare the generated U/V fields with an existing result, one vertical
    level at a time, and include the statistics in the summary JSON file.

``--overwrite``
    Permit replacement of an existing ``--output`` file after the new file
    has been completed successfully.  Without this option, the program stops
    if the output already exists.  Never set ``--input`` and ``--output`` to
    the same path.

Run ``python reconstruction_current.py --help`` for all available options.

Method
------

This is a Python implementation of the MATLAB workflow in
``reconstruction_current.m``.  The numerical logic is intentionally kept the
same:

1. Compute in-situ density with the Jackett and McDougall (1995) polynomial
   used by ``rho_eos.m``.
2. Compute SSH-referenced geostrophic currents and thermal-wind shear on MOM6
   tracer cells.
3. Interpolate the currents to MOM6 U/V faces.
4. Solve the same depth-integrated Poisson problem with conjugate gradients.
5. Add the depth-independent velocity correction and overwrite ``u`` and
   ``v`` in a copy of the original IC file.

Unlike the MATLAB implementation, three-dimensional fields are processed one
vertical level at a time.  This avoids keeping several multi-gigabyte C3200
or C9600 arrays in memory at once.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
from netCDF4 import Dataset
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import cg


RHO0 = 1025.0
GRAVITY = 9.8
OMEGA = 7.292115e-5


@dataclass
class GridMetrics:
    """MOM6 native-grid distances in Python ``(y, x)`` order."""

    dx_t: np.ndarray
    dy_t: np.ndarray
    dx_u: np.ndarray
    dy_u: np.ndarray
    dx_v: np.ndarray
    dy_v: np.ndarray


@dataclass
class SolverDiagnostics:
    """Diagnostics from the barotropic correction solve."""

    wet_points: int
    nonzeros: int
    iterations: int
    info: int
    relative_residual: float
    divergence_std_before: float
    divergence_std_theory: float
    divergence_std_after: float
    divergence_ratio: float


@dataclass
class PlotDiagnostics:
    """Small fields retained for diagnostic plotting."""

    lon: np.ndarray
    lat: np.ndarray
    z: np.ndarray
    plot_level: int
    mask_t: np.ndarray
    original_u_surface: np.ndarray
    original_v_surface: np.ndarray
    geostrophic_u_surface: np.ndarray
    geostrophic_v_surface: np.ndarray
    adjusted_u_surface: np.ndarray
    adjusted_v_surface: np.ndarray
    correction_u_t: np.ndarray
    correction_v_t: np.ndarray
    original_u_section: np.ndarray
    geostrophic_u_section: np.ndarray
    adjusted_u_section: np.ndarray
    section_wet: np.ndarray
    section_x_index: int
    rhs_t: np.ndarray
    divergence_after: np.ndarray


def log(message: str) -> None:
    """Print an immediately visible progress message."""

    print(message, flush=True)


def elapsed(label: str, start: float) -> None:
    """Print elapsed wall time for one phase."""

    log(f"[TIMING] {label}: {perf_counter() - start:.2f} seconds")


def publish_output(partial_path: Path, output_path: Path) -> None:
    """Safely publish a completed file copied from another filesystem.

    Copying directly to ``output_path`` can leave a truncated file with the
    final name when a batch job is interrupted.  Copy to a process-specific
    name in the destination directory first, verify its size, and only then
    atomically replace the final path.
    """

    incoming_path = output_path.with_name(
        f"{output_path.name}.incoming.{os.getpid()}"
    )
    try:
        shutil.copy2(partial_path, incoming_path)
        source_size = partial_path.stat().st_size
        copied_size = incoming_path.stat().st_size
        if copied_size != source_size:
            raise OSError(
                f"Incomplete output copy: {copied_size} of "
                f"{source_size} bytes"
            )
        os.replace(incoming_path, output_path)
    finally:
        if incoming_path.exists():
            incoming_path.unlink()

    partial_path.unlink()


def read_array(data: Any) -> np.ndarray:
    """Return NetCDF data as float64, replacing masked values with NaN."""

    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    return np.asarray(data, dtype=np.float64)


def rho_eos(temperature: np.ndarray, salinity: np.ndarray, z: float) -> np.ndarray:
    """Port of ``rho_eos.m`` for one vertical level.

    Parameters
    ----------
    temperature
        Potential temperature in degrees Celsius.
    salinity
        Practical salinity in PSU.
    z
        Negative layer-center depth in meters.
    """

    a00 = +19092.56
    a01 = +209.8925
    a02 = -3.041638
    a03 = -1.852732e-3
    a04 = -1.361629e-5
    a10 = 104.4077
    a11 = -6.500517
    a12 = +0.1553190
    a13 = 2.326469e-4
    as0 = -5.587545
    as1 = +0.7390729
    as2 = -1.909078e-2
    b00 = +4.721788e-1
    b01 = +1.028859e-2
    b02 = -2.512549e-4
    b03 = -5.939910e-7
    b10 = -1.571896e-2
    b11 = -2.598241e-4
    b12 = +7.267926e-6
    bs1 = +2.042967e-3
    e00 = +1.045941e-5
    e01 = -5.782165e-10
    e02 = +1.296821e-7
    e10 = -2.595994e-7
    e11 = -1.248266e-9
    e12 = -3.508914e-9

    qr = +999.842594
    q01 = +6.793952e-2
    q02 = -9.095290e-3
    q03 = +1.001685e-4
    q04 = -1.120083e-6
    q05 = +6.536332e-9
    q10 = +0.824493
    q11 = -4.08990e-3
    q12 = +7.64380e-5
    q13 = -8.24670e-7
    q14 = +5.38750e-9
    qs0 = -5.72466e-3
    qs1 = +1.02270e-4
    qs2 = -1.65460e-6
    q20 = +4.8314e-4

    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        sqrt_s = np.sqrt(salinity)

        k0 = (
            a00
            + temperature
            * (
                a01
                + temperature
                * (a02 + temperature * (a03 + temperature * a04))
            )
            + salinity
            * (
                a10
                + temperature
                * (a11 + temperature * (a12 + temperature * a13))
                + sqrt_s * (as0 + temperature * (as1 + temperature * as2))
            )
        )

        k1 = (
            b00
            + temperature
            * (b01 + temperature * (b02 + temperature * b03))
            + salinity
            * (
                b10
                + temperature * (b11 + temperature * b12)
                + sqrt_s * bs1
            )
        )

        k2 = (
            e00
            + temperature * (e01 + temperature * e02)
            + salinity
            * (
                e10
                + temperature * (e11 + temperature * e12)
            )
        )

        rho1 = (
            qr
            + temperature
            * (
                q01
                + temperature
                * (
                    q02
                    + temperature
                    * (
                        q03
                        + temperature * (q04 + temperature * q05)
                    )
                )
            )
            + salinity
            * (
                q10
                + temperature
                * (
                    q11
                    + temperature
                    * (
                        q12
                        + temperature
                        * (q13 + temperature * q14)
                    )
                )
                + sqrt_s
                * (qs0 + temperature * (qs1 + temperature * qs2))
                + salinity * q20
            )
        )

        rho = rho1 / (1.0 + 0.1 * z / (k0 - z * (k1 - z * k2)))

    return rho


def build_grid_metrics(
    hgrid_path: Path,
    ny: int,
    nx: int,
) -> GridMetrics:
    """Combine MOM6 supergrid distances into native T/U/V metrics."""

    start = perf_counter()
    with Dataset(hgrid_path, "r") as hgrid:
        dx_sg = read_array(hgrid.variables["dx"][:])
        dy_sg = read_array(hgrid.variables["dy"][:])

    expected_dx = (2 * ny + 1, 2 * nx)
    expected_dy = (2 * ny, 2 * nx + 1)
    if dx_sg.shape != expected_dx:
        raise ValueError(
            f"Unexpected dx shape {dx_sg.shape}; expected {expected_dx}"
        )
    if dy_sg.shape != expected_dy:
        raise ValueError(
            f"Unexpected dy shape {dy_sg.shape}; expected {expected_dy}"
        )

    # U-face length in the model-y direction: (ny, nx+1).
    dy_u = dy_sg[0::2, 0::2] + dy_sg[1::2, 0::2]

    # V-face length in the model-x direction: (ny+1, nx).
    dx_v = dx_sg[0::2, 0::2] + dx_sg[0::2, 1::2]

    # Distance between adjacent T points across U faces: (ny, nx+1).
    dx_u = np.zeros((ny, nx + 1), dtype=np.float64)
    dx_u[:, 1:nx] = (
        dx_sg[1::2, 1:-1:2] + dx_sg[1::2, 2::2]
    )
    dx_u[:, 0] = dx_sg[1::2, 0]
    dx_u[:, nx] = dx_sg[1::2, -1]

    # Distance between adjacent T points across V faces: (ny+1, nx).
    dy_v = np.zeros((ny + 1, nx), dtype=np.float64)
    dy_v[1:ny, :] = (
        dy_sg[1:-1:2, 1::2] + dy_sg[2::2, 1::2]
    )
    dy_v[0, :] = dy_sg[0, 1::2]
    dy_v[ny, :] = dy_sg[-1, 1::2]

    # Arrays expected by the original finite-difference formulas.
    dx_t = np.zeros((ny, nx), dtype=np.float64)
    dx_t[:, : nx - 1] = dx_u[:, 1:nx]
    dx_t[:, nx - 1] = dx_u[:, nx - 1]

    dy_t = np.zeros((ny, nx), dtype=np.float64)
    dy_t[: ny - 1, :] = dy_v[1:ny, :]
    dy_t[ny - 1, :] = dy_v[ny - 1, :]

    log(
        f"dx_T range: {np.nanmin(dx_t):.2f} - "
        f"{np.nanmax(dx_t):.2f} m"
    )
    log(
        f"dy_T range: {np.nanmin(dy_t):.2f} - "
        f"{np.nanmax(dy_t):.2f} m"
    )
    elapsed("read/combine horizontal metrics", start)

    return GridMetrics(
        dx_t=dx_t,
        dy_t=dy_t,
        dx_u=dx_u,
        dy_u=dy_u,
        dx_v=dx_v,
        dy_v=dy_v,
    )


def make_face_masks(mask_t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Create the same two-dimensional U/V masks as the MATLAB code."""

    ny, nx = mask_t.shape
    umask = np.zeros((ny, nx + 1), dtype=np.float64)
    vmask = np.zeros((ny + 1, nx), dtype=np.float64)

    umask[:, 1:nx] = mask_t[:, : nx - 1] * mask_t[:, 1:nx]
    umask[:, 0] = mask_t[:, 0]
    umask[:, nx] = mask_t[:, nx - 1]

    vmask[1:ny, :] = mask_t[: ny - 1, :] * mask_t[1:ny, :]
    vmask[0, :] = mask_t[0, :]
    vmask[ny, :] = mask_t[ny - 1, :]

    return umask, vmask


def surface_geostrophic_current(
    ssh: np.ndarray,
    mask_t: np.ndarray,
    metrics: GridMetrics,
    f2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the SSH-referenced current on tracer cells."""

    ny, nx = ssh.shape
    wet = mask_t != 0.0
    dssh_dx = np.zeros((ny, nx), dtype=np.float64)
    dssh_dy = np.zeros((ny, nx), dtype=np.float64)

    valid_x = wet[:, : nx - 2] & wet[:, 2:nx]
    numerator_x = ssh[:, 2:nx] - ssh[:, : nx - 2]
    denominator_x = (
        metrics.dx_t[:, 1 : nx - 1]
        + metrics.dx_t[:, : nx - 2]
    )
    dssh_dx[:, 1 : nx - 1] = np.where(
        valid_x,
        numerator_x / denominator_x,
        0.0,
    )

    valid_y = wet[: ny - 2, :] & wet[2:ny, :]
    numerator_y = ssh[2:ny, :] - ssh[: ny - 2, :]
    denominator_y = (
        metrics.dy_t[1 : ny - 1, :]
        + metrics.dy_t[: ny - 2, :]
    )
    dssh_dy[1 : ny - 1, :] = np.where(
        valid_y,
        numerator_y / denominator_y,
        0.0,
    )

    valid = wet[:, 0] & wet[:, 1]
    dssh_dx[:, 0] = np.where(
        valid,
        (ssh[:, 1] - ssh[:, 0]) / metrics.dx_t[:, 0],
        0.0,
    )
    valid = wet[:, nx - 2] & wet[:, nx - 1]
    dssh_dx[:, nx - 1] = np.where(
        valid,
        (ssh[:, nx - 1] - ssh[:, nx - 2])
        / metrics.dx_t[:, nx - 1],
        0.0,
    )

    valid = wet[0, :] & wet[1, :]
    dssh_dy[0, :] = np.where(
        valid,
        (ssh[1, :] - ssh[0, :]) / metrics.dy_t[0, :],
        0.0,
    )
    valid = wet[ny - 2, :] & wet[ny - 1, :]
    dssh_dy[ny - 1, :] = np.where(
        valid,
        (ssh[ny - 1, :] - ssh[ny - 2, :])
        / metrics.dy_t[ny - 1, :],
        0.0,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        u_ssh = -(GRAVITY / f2d) * dssh_dy * mask_t
        v_ssh = +(GRAVITY / f2d) * dssh_dx * mask_t

    u_ssh[~np.isfinite(u_ssh)] = 0.0
    v_ssh[~np.isfinite(v_ssh)] = 0.0
    return u_ssh, v_ssh


def thermal_wind_shear(
    rho: np.ndarray,
    wet: np.ndarray,
    metrics: GridMetrics,
    g_over_rho0_f: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute ``du/dz`` and ``dv/dz`` on one tracer level."""

    ny, nx = rho.shape
    du_dz = np.zeros((ny, nx), dtype=np.float64)
    dv_dz = np.zeros((ny, nx), dtype=np.float64)

    valid_x = wet[:, : nx - 2] & wet[:, 2:nx]
    drho_dx = (
        (rho[:, 2:nx] - rho[:, : nx - 2])
        / (
            metrics.dx_t[:, 1 : nx - 1]
            + metrics.dx_t[:, : nx - 2]
        )
    )
    dv_dz[:, 1 : nx - 1] = np.where(
        valid_x,
        g_over_rho0_f[:, 1 : nx - 1] * drho_dx,
        0.0,
    )

    valid_y = wet[: ny - 2, :] & wet[2:ny, :]
    drho_dy = (
        (rho[2:ny, :] - rho[: ny - 2, :])
        / (
            metrics.dy_t[1 : ny - 1, :]
            + metrics.dy_t[: ny - 2, :]
        )
    )
    du_dz[1 : ny - 1, :] = np.where(
        valid_y,
        -g_over_rho0_f[1 : ny - 1, :] * drho_dy,
        0.0,
    )

    du_dz[~np.isfinite(du_dz)] = 0.0
    dv_dz[~np.isfinite(dv_dz)] = 0.0
    return du_dz, dv_dz


def tracer_to_faces(
    u_t: np.ndarray,
    v_t: np.ndarray,
    wet: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate tracer currents to U/V faces using wet-neighbor masks."""

    ny, nx = wet.shape
    u_face = np.zeros((ny, nx + 1), dtype=np.float64)
    v_face = np.zeros((ny + 1, nx), dtype=np.float64)

    wet_u = wet[:, : nx - 1] & wet[:, 1:nx]
    u_face[:, 1:nx] = np.where(
        wet_u,
        0.5 * (u_t[:, : nx - 1] + u_t[:, 1:nx]),
        0.0,
    )

    wet_v = wet[: ny - 1, :] & wet[1:ny, :]
    v_face[1:ny, :] = np.where(
        wet_v,
        0.5 * (v_t[: ny - 1, :] + v_t[1:ny, :]),
        0.0,
    )

    u_face[~np.isfinite(u_face)] = 0.0
    v_face[~np.isfinite(v_face)] = 0.0
    return u_face, v_face


def thickness_to_faces(dz_t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate layer thickness from tracer cells to U/V faces."""

    ny, nx = dz_t.shape
    dz_u = np.empty((ny, nx + 1), dtype=np.float64)
    dz_v = np.empty((ny + 1, nx), dtype=np.float64)

    dz_u[:, 1:nx] = 0.5 * (dz_t[:, : nx - 1] + dz_t[:, 1:nx])
    dz_u[:, 0] = dz_t[:, 0]
    dz_u[:, nx] = dz_t[:, nx - 1]

    dz_v[1:ny, :] = 0.5 * (dz_t[: ny - 1, :] + dz_t[1:ny, :])
    dz_v[0, :] = dz_t[0, :]
    dz_v[ny, :] = dz_t[ny - 1, :]

    return dz_u, dz_v


def face_depths(h_t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate total water-column thickness to U/V faces."""

    ny, nx = h_t.shape
    he = np.empty((ny, nx + 1), dtype=np.float64)
    hn = np.empty((ny + 1, nx), dtype=np.float64)

    he[:, 1:nx] = 0.5 * (h_t[:, : nx - 1] + h_t[:, 1:nx])
    he[:, 0] = h_t[:, 0]
    he[:, nx] = h_t[:, nx - 1]

    hn[1:ny, :] = 0.5 * (h_t[: ny - 1, :] + h_t[1:ny, :])
    hn[0, :] = h_t[0, :]
    hn[ny, :] = h_t[ny - 1, :]
    return he, hn


def build_poisson_matrix(
    h_t: np.ndarray,
    mask_eff: np.ndarray,
    area: np.ndarray,
    metrics: GridMetrics,
) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
    """Vectorized construction of the MATLAB five-point sparse operator."""

    start = perf_counter()
    ny, nx = h_t.shape
    he, hn = face_depths(h_t)

    wet_flat = np.flatnonzero(mask_eff.ravel(order="C"))
    nw = wet_flat.size
    if nw == 0:
        raise ValueError("No wet cells are available for the Poisson solve")
    if nw >= np.iinfo(np.int32).max:
        raise ValueError("Poisson system exceeds int32 sparse-index capacity")

    index_map = np.full(mask_eff.size, -1, dtype=np.int32)
    index_map[wet_flat] = np.arange(nw, dtype=np.int32)
    index_map = index_map.reshape(mask_eff.shape)

    rows = np.empty(5 * nw, dtype=np.int32)
    cols = np.empty(5 * nw, dtype=np.int32)
    vals = np.empty(5 * nw, dtype=np.float64)
    diagonal = np.zeros(nw, dtype=np.float64)
    pointer = 0

    def append_direction(
        source: np.ndarray,
        neighbor: np.ndarray,
        coefficient: np.ndarray,
    ) -> None:
        nonlocal pointer
        count = source.size
        rows[pointer : pointer + count] = source
        cols[pointer : pointer + count] = neighbor
        vals[pointer : pointer + count] = -coefficient
        diagonal[source] += coefficient
        pointer += count

    # East: cell (y, x) to (y, x+1).
    valid = mask_eff[:, : nx - 1] & mask_eff[:, 1:nx]
    coefficient = (
        he[:, 1:nx]
        * metrics.dy_u[:, 1:nx]
        / (metrics.dx_u[:, 1:nx] * area[:, : nx - 1])
    )
    append_direction(
        index_map[:, : nx - 1][valid],
        index_map[:, 1:nx][valid],
        coefficient[valid],
    )

    # West: cell (y, x) to (y, x-1).
    coefficient = (
        he[:, 1:nx]
        * metrics.dy_u[:, 1:nx]
        / (metrics.dx_u[:, 1:nx] * area[:, 1:nx])
    )
    append_direction(
        index_map[:, 1:nx][valid],
        index_map[:, : nx - 1][valid],
        coefficient[valid],
    )

    # North: cell (y, x) to (y+1, x).
    valid = mask_eff[: ny - 1, :] & mask_eff[1:ny, :]
    coefficient = (
        hn[1:ny, :]
        * metrics.dx_v[1:ny, :]
        / (metrics.dy_v[1:ny, :] * area[: ny - 1, :])
    )
    append_direction(
        index_map[: ny - 1, :][valid],
        index_map[1:ny, :][valid],
        coefficient[valid],
    )

    # South: cell (y, x) to (y-1, x).
    coefficient = (
        hn[1:ny, :]
        * metrics.dx_v[1:ny, :]
        / (metrics.dy_v[1:ny, :] * area[1:ny, :])
    )
    append_direction(
        index_map[1:ny, :][valid],
        index_map[: ny - 1, :][valid],
        coefficient[valid],
    )

    diagonal_indices = np.arange(nw, dtype=np.int32)
    rows[pointer : pointer + nw] = diagonal_indices
    cols[pointer : pointer + nw] = diagonal_indices
    vals[pointer : pointer + nw] = diagonal
    pointer += nw

    matrix = coo_matrix(
        (vals[:pointer], (rows[:pointer], cols[:pointer])),
        shape=(nw, nw),
    ).tocsr()
    matrix.sum_duplicates()

    elapsed("build Poisson matrix", start)
    log(
        f"Poisson matrix: wet_points={nw}, nonzeros={matrix.nnz}"
    )
    return matrix, wet_flat, diagonal


def solve_poisson(
    matrix: csr_matrix,
    rhs: np.ndarray,
    rtol: float,
    maxiter: int,
) -> tuple[np.ndarray, int, int, float]:
    """Solve with the same unpreconditioned conjugate-gradient method."""

    start = perf_counter()
    iterations = 0

    def count_iteration(_: np.ndarray) -> None:
        nonlocal iterations
        iterations += 1

    solution, info = cg(
        matrix,
        rhs,
        x0=None,
        rtol=rtol,
        atol=0.0,
        maxiter=maxiter,
        callback=count_iteration,
    )
    residual = np.linalg.norm(matrix @ solution - rhs)
    relative_residual = residual / max(np.linalg.norm(rhs), np.finfo(float).eps)

    log(
        "cg: "
        f"info={info}, relative_residual={relative_residual:.6e}, "
        f"iterations={iterations}, ||x||={np.linalg.norm(solution):.6e}"
    )
    elapsed("Poisson CG solve", start)
    return solution, info, iterations, float(relative_residual)


def moving_nanmean(data: np.ndarray, window: int = 5) -> np.ndarray:
    """Centered moving mean along the last axis, ignoring NaNs."""

    if window <= 1:
        return data.copy()
    kernel = np.ones(window, dtype=np.float64)
    result = np.empty_like(data, dtype=np.float64)
    for index in range(data.shape[0]):
        values = data[index]
        valid = np.isfinite(values)
        numerator = np.convolve(
            np.where(valid, values, 0.0),
            kernel,
            mode="same",
        )
        denominator = np.convolve(
            valid.astype(np.float64),
            kernel,
            mode="same",
        )
        result[index] = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan),
            where=denominator > 0,
        )
    return result


def make_diagnostic_plots(
    data: PlotDiagnostics,
    plot_dir: Path,
    plot_format: str = "png",
) -> None:
    """Create diagnostic plots corresponding to the MATLAB figures."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    # Ignore user/site matplotlibrc settings that can request extremely large
    # fonts on an HPC node and trigger FreeType "raster overflow" errors.
    plt.rcdefaults()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 14,
            "text.usetex": False,
            "svg.fonttype": "none",
        }
    )

    plot_dir.mkdir(parents=True, exist_ok=True)
    start = perf_counter()

    def save_figure(fig: Any, stem: str) -> Path:
        """Save a figure and fall back from broken HPC PNG fonts to SVG."""

        requested_path = plot_dir / f"{stem}.{plot_format}"
        try:
            if plot_format == "png":
                fig.savefig(requested_path, dpi=120)
            else:
                fig.savefig(requested_path, format="svg")
            return requested_path
        except RuntimeError as error:
            message = str(error)
            freetype_error = (
                plot_format == "png"
                and (
                    "FT_Render_Glyph" in message
                    or "raster overflow" in message
                )
            )
            if not freetype_error:
                raise

            fallback_path = plot_dir / f"{stem}.svg"
            log(
                "[WARNING] PNG font rendering failed; "
                f"saving SVG instead: {fallback_path}"
            )
            fig.savefig(fallback_path, format="svg")
            return fallback_path

    def speed(u: np.ndarray, v: np.ndarray) -> np.ndarray:
        return np.sqrt(u * u + v * v)

    mask = data.mask_t != 0

    # Full C3200 maps contain more than two million cells per panel, and
    # C9600 contains about nine times more.  Downsample only the diagnostic
    # maps to keep Matplotlib memory and rasterization bounded.  All numerical
    # calculations and NetCDF output continue to use the full grid.
    maximum_map_cells = 250_000
    map_stride = max(
        1,
        int(math.ceil(math.sqrt(mask.size / maximum_map_cells))),
    )
    maximum_arrows = 4_000
    quiver_stride = max(
        map_stride,
        int(math.ceil(math.sqrt(mask.size / maximum_arrows))),
    )
    map_y = slice(None, None, map_stride)
    map_x = slice(None, None, map_stride)
    quiver_y = slice(None, None, quiver_stride)
    quiver_x = slice(None, None, quiver_stride)
    map_lon = data.lon[map_x]
    map_lat = data.lat[map_y]
    map_mask = mask[map_y, map_x]

    # Calculate plotting quantities only at sampled points.  This avoids
    # allocating several additional full-grid C9600 arrays just for figures.
    original_speed = np.where(
        map_mask,
        speed(
            data.original_u_surface[map_y, map_x],
            data.original_v_surface[map_y, map_x],
        ),
        np.nan,
    )
    geo_speed = np.where(
        map_mask,
        speed(
            data.geostrophic_u_surface[map_y, map_x],
            data.geostrophic_v_surface[map_y, map_x],
        ),
        np.nan,
    )
    adjusted_speed = np.where(
        map_mask,
        speed(
            data.adjusted_u_surface[map_y, map_x],
            data.adjusted_v_surface[map_y, map_x],
        ),
        np.nan,
    )
    speed_difference = original_speed - geo_speed

    log(
        "Diagnostic map sampling: "
        f"every {map_stride} grid point(s); "
        f"quiver every {quiver_stride} grid point(s)"
    )

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    surface_panels = (
        ("Original speed", original_speed, 0.0, 2.0),
        ("Geostrophic speed", geo_speed, 0.0, 2.0),
        ("Original - geostrophic", speed_difference, -1.0, 1.0),
        ("Adjusted speed", adjusted_speed, 0.0, 2.0),
    )
    for axis, (title, field, vmin, vmax) in zip(
        axes.flat,
        surface_panels,
        strict=True,
    ):
        image = axis.pcolormesh(
            map_lon,
            map_lat,
            field,
            shading="auto",
            rasterized=True,
            cmap="RdYlBu_r" if vmin >= 0 else "RdBu_r",
            vmin=vmin,
            vmax=vmax,
        )
        axis.set_title(title)
        axis.set_xlabel("Longitude")
        axis.set_ylabel("Latitude")
        fig.colorbar(image, ax=axis, label="m s$^{-1}$")

    axes[0, 0].quiver(
        data.lon[quiver_x],
        data.lat[quiver_y],
        data.original_u_surface[quiver_y, quiver_x],
        data.original_v_surface[quiver_y, quiver_x],
        color="k",
        scale=35,
    )
    axes[0, 1].quiver(
        data.lon[quiver_x],
        data.lat[quiver_y],
        data.geostrophic_u_surface[quiver_y, quiver_x],
        data.geostrophic_v_surface[quiver_y, quiver_x],
        color="k",
        scale=35,
    )
    axes[1, 1].quiver(
        data.lon[quiver_x],
        data.lat[quiver_y],
        data.adjusted_u_surface[quiver_y, quiver_x],
        data.adjusted_v_surface[quiver_y, quiver_x],
        color="k",
        scale=35,
    )
    depth = -data.z[data.plot_level]
    fig.suptitle(f"Surface-current diagnostics at {depth:.2f} m")
    save_figure(fig, "surface_currents")
    plt.close(fig)

    section_fields = (
        ("Original zonal current", data.original_u_section),
        ("Geostrophic zonal current", data.geostrophic_u_section),
        ("Adjusted zonal current", data.adjusted_u_section),
    )
    fig, axes = plt.subplots(3, 1, figsize=(15, 12), constrained_layout=True)
    levels = np.arange(-0.8, 1.21, 0.1)
    for axis, (title, section) in zip(
        axes,
        section_fields,
        strict=True,
    ):
        section = moving_nanmean(
            np.where(data.section_wet, section, np.nan),
            window=5,
        )
        contour = axis.contourf(
            data.lat,
            data.z,
            section,
            levels=levels,
            cmap="Spectral_r",
            extend="both",
        )
        axis.contour(
            data.lat,
            data.z,
            section,
            levels=np.arange(-0.8, 1.21, 0.2),
            colors="0.35",
            linewidths=0.5,
        )
        axis.set_ylim(-2000.0, 0.0)
        axis.set_title(title)
        axis.set_xlabel("Latitude")
        axis.set_ylabel("Depth (m)")
        fig.colorbar(contour, ax=axis, label="m s$^{-1}$")
    fig.suptitle(
        f"Zonal-current section at xh index {data.section_x_index}"
    )
    save_figure(fig, "zonal_current_sections")
    plt.close(fig)

    correction_speed = speed(
        data.correction_u_t[map_y, map_x],
        data.correction_v_t[map_y, map_x],
    )
    correction_speed = np.where(map_mask, correction_speed, np.nan)
    fig, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
    image = axis.pcolormesh(
        map_lon,
        map_lat,
        correction_speed,
        shading="auto",
        rasterized=True,
        cmap="RdYlBu_r",
    )
    axis.quiver(
        data.lon[quiver_x],
        data.lat[quiver_y],
        data.correction_u_t[quiver_y, quiver_x],
        data.correction_v_t[quiver_y, quiver_x],
        color="k",
        scale=5,
    )
    axis.set_title("Depth-independent velocity correction")
    axis.set_xlabel("Longitude")
    axis.set_ylabel("Latitude")
    fig.colorbar(image, ax=axis, label="m s$^{-1}$")
    save_figure(fig, "barotropic_correction")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(13, 10), constrained_layout=True)
    for axis, title, field in (
        (axes[0], "Divergence before", data.rhs_t),
        (axes[1], "Divergence after", data.divergence_after),
    ):
        image = axis.pcolormesh(
            map_lon,
            map_lat,
            np.where(map_mask, field[map_y, map_x], np.nan),
            shading="auto",
            rasterized=True,
            cmap="RdBu_r",
            vmin=-0.01,
            vmax=0.01,
        )
        axis.set_title(title)
        axis.set_xlabel("Longitude")
        axis.set_ylabel("Latitude")
        fig.colorbar(image, ax=axis)
    save_figure(fig, "divergence_before_after")
    plt.close(fig)
    elapsed("diagnostic plots", start)


def compare_uv_files(
    candidate_path: Path,
    reference_path: Path,
    time_index: int = 0,
) -> dict[str, dict[str, float | int]]:
    """Compare U/V fields level-by-level without loading full files."""

    results: dict[str, dict[str, float | int]] = {}
    with Dataset(candidate_path, "r") as candidate, Dataset(
        reference_path,
        "r",
    ) as reference:
        for name in ("u", "v"):
            candidate_var = candidate.variables[name]
            reference_var = reference.variables[name]
            if candidate_var.shape != reference_var.shape:
                raise ValueError(
                    f"{name} shape mismatch: {candidate_var.shape} vs "
                    f"{reference_var.shape}"
                )

            count = 0
            sum_squared = 0.0
            reference_squared = 0.0
            maximum = 0.0
            for level in range(candidate_var.shape[1]):
                lhs = read_array(candidate_var[time_index, level, ...])
                rhs = read_array(reference_var[time_index, level, ...])
                valid = np.isfinite(lhs) & np.isfinite(rhs)
                difference = lhs[valid] - rhs[valid]
                count += difference.size
                if difference.size:
                    sum_squared += float(np.dot(difference, difference))
                    reference_values = rhs[valid]
                    reference_squared += float(
                        np.dot(reference_values, reference_values)
                    )
                    maximum = max(
                        maximum,
                        float(np.max(np.abs(difference))),
                    )

            rms = math.sqrt(sum_squared / max(count, 1))
            relative_rms = math.sqrt(
                sum_squared / max(reference_squared, np.finfo(float).eps)
            )
            results[name] = {
                "count": count,
                "max_abs": maximum,
                "rms": rms,
                "relative_rms": relative_rms,
            }
    return results


def reconstruct(
    input_path: Path,
    output_path: Path,
    grid_dir: Path,
    work_dir: Path | None,
    time_index: int,
    plot_level: int,
    section_x_index: int | None,
    cg_rtol: float,
    cg_maxiter: int,
    require_convergence: bool,
    apply_barotropic_correction: bool,
    overwrite: bool,
) -> tuple[SolverDiagnostics, PlotDiagnostics]:
    """Run the complete reconstruction and write the adjusted IC file."""

    total_start = perf_counter()
    hgrid_path = grid_dir / "ocean_hgrid.nc"
    mask_path = grid_dir / "ocean_mask.nc"
    topog_path = grid_dir / "topog.nc"
    for path in (input_path, hgrid_path, mask_path, topog_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"{output_path} already exists; use --overwrite to replace it"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if work_dir is None:
        partial_path = output_path.with_name(output_path.name + ".partial")
    else:
        work_dir.mkdir(parents=True, exist_ok=True)
        partial_path = work_dir / (output_path.name + ".partial")
    if partial_path.exists():
        partial_path.unlink()

    with Dataset(input_path, "r") as source:
        required = {
            "temp",
            "salt",
            "ssh",
            "u",
            "v",
            "zl",
            "xh",
            "yh",
            "xq",
            "yq",
        }
        missing = sorted(required.difference(source.variables))
        if missing:
            raise KeyError(f"Missing IC variables: {missing}")

        nt, nk, ny, nx = source.variables["temp"].shape
        if source.variables["u"].shape != (nt, nk, ny, nx + 1):
            raise ValueError("Unexpected MOM6 u dimensions")
        if source.variables["v"].shape != (nt, nk, ny + 1, nx):
            raise ValueError("Unexpected MOM6 v dimensions")
        if not 0 <= time_index < nt:
            raise IndexError(f"time-index {time_index} outside [0, {nt})")
        if not 0 <= plot_level < nk:
            raise IndexError(f"plot-level {plot_level} outside [0, {nk})")

        if section_x_index is None:
            section_x_index = min(1294, nx - 1)
        if not 0 <= section_x_index < nx:
            raise IndexError(
                f"section-x-index {section_x_index} outside [0, {nx})"
            )

        lon = read_array(source.variables["xh"][:])
        lat = read_array(source.variables["yh"][:])
        z = -read_array(source.variables["zl"][:])

    log(
        f"IC dimensions: time={nt}, zl={nk}, yh={ny}, xh={nx}, "
        f"xq={nx + 1}, yq={ny + 1}"
    )

    with Dataset(mask_path, "r") as mask_ds:
        mask_t = read_array(mask_ds.variables["mask"][:])
        area = read_array(mask_ds.variables["areaO"][:])
    with Dataset(topog_path, "r") as topog_ds:
        depth = read_array(topog_ds.variables["depth"][:])

    expected_shape = (ny, nx)
    for name, array in (
        ("mask", mask_t),
        ("areaO", area),
        ("depth", depth),
    ):
        if array.shape != expected_shape:
            raise ValueError(
                f"{name} shape {array.shape}; expected {expected_shape}"
            )

    metrics = build_grid_metrics(hgrid_path, ny=ny, nx=nx)
    umask, vmask = make_face_masks(mask_t)

    f2d = 2.0 * OMEGA * np.sin(np.deg2rad(lat))[:, np.newaxis]
    f2d = np.broadcast_to(f2d, (ny, nx))
    g_over_rho0_f = GRAVITY / (RHO0 * f2d)

    interfaces = np.empty(nk + 1, dtype=np.float64)
    interfaces[0] = 0.0
    interfaces[1:nk] = 0.5 * (z[: nk - 1] + z[1:nk])
    interfaces[nk] = z[nk - 1] + (z[nk - 1] - z[nk - 2])
    dz_geom = interfaces[:nk] - interfaces[1 : nk + 1]
    if np.any(dz_geom <= 0.0):
        raise ValueError("Computed non-positive geometric layer thickness")

    log(f"Copying input IC to temporary output: {partial_path}")
    copy_start = perf_counter()
    shutil.copy2(input_path, partial_path)
    elapsed("copy IC file", copy_start)

    h_t = np.zeros((ny, nx), dtype=np.float64)
    hu = np.zeros((ny, nx + 1), dtype=np.float64)
    hv = np.zeros((ny + 1, nx), dtype=np.float64)

    u_rel = np.zeros((ny, nx), dtype=np.float64)
    v_rel = np.zeros((ny, nx), dtype=np.float64)
    previous_du_dz: np.ndarray | None = None
    previous_dv_dz: np.ndarray | None = None
    previous_wet: np.ndarray | None = None

    original_u_surface = np.zeros((ny, nx), dtype=np.float64)
    original_v_surface = np.zeros((ny, nx), dtype=np.float64)
    geostrophic_u_surface = np.zeros((ny, nx), dtype=np.float64)
    geostrophic_v_surface = np.zeros((ny, nx), dtype=np.float64)
    original_u_section = np.zeros((nk, ny), dtype=np.float64)
    geostrophic_u_section = np.zeros((nk, ny), dtype=np.float64)
    section_wet = np.zeros((nk, ny), dtype=bool)

    try:
        first_pass_start = perf_counter()
        with Dataset(input_path, "r") as source, Dataset(
            partial_path,
            "r+",
        ) as destination:
            temp_var = source.variables["temp"]
            salt_var = source.variables["salt"]
            ssh_var = source.variables["ssh"]
            source_u = source.variables["u"]
            source_v = source.variables["v"]
            destination_u = destination.variables["u"]
            destination_v = destination.variables["v"]

            ssh = read_array(ssh_var[time_index, :, :])
            u_ssh, v_ssh = surface_geostrophic_current(
                ssh,
                mask_t,
                metrics,
                f2d,
            )
            log(
                "[INFO] Preserving original currents next to coasts "
                "and in the bottom wet layer."
            )

            for level in range(nk):
                wet = (mask_t == 1.0) & (z[level] > -depth)
                safe_t = np.zeros_like(wet)
                safe_t[1:-1, 1:-1] = (
                    wet[1:-1, 1:-1]
                    & wet[1:-1, :-2]
                    & wet[1:-1, 2:]
                    & wet[:-2, 1:-1]
                    & wet[2:, 1:-1]
                )
                if level + 1 < nk:
                    safe_t &= z[level + 1] > -depth

                dz_t = dz_geom[level] * wet
                h_t += dz_t

                if level > 0:
                    if (
                        previous_du_dz is None
                        or previous_dv_dz is None
                        or previous_wet is None
                    ):
                        raise RuntimeError("Missing previous-level shear")
                    previous_dz = dz_geom[level - 1] * previous_wet
                    u_rel += previous_du_dz * previous_dz
                    v_rel += previous_dv_dz * previous_dz

                # Equivalent to the final MATLAB dry-cell mask because the
                # water column is vertically contiguous from the surface.
                u_rel[~wet] = 0.0
                v_rel[~wet] = 0.0

                temperature = read_array(
                    temp_var[time_index, level, :, :]
                )
                salinity = read_array(
                    salt_var[time_index, level, :, :]
                )
                rho = rho_eos(temperature, salinity, float(z[level]))
                du_dz, dv_dz = thermal_wind_shear(
                    rho,
                    safe_t,
                    metrics,
                    g_over_rho0_f,
                )

                u_t = (u_rel + u_ssh) * wet
                v_t = (v_rel + v_ssh) * wet
                u_geo, v_geo = tracer_to_faces(u_t, v_t, wet)

                original_u = read_array(
                    source_u[time_index, level, :, :]
                )
                original_v = read_array(
                    source_v[time_index, level, :, :]
                )
                safe_u = np.zeros((ny, nx + 1), dtype=bool)
                safe_v = np.zeros((ny + 1, nx), dtype=bool)
                safe_u[:, 1:nx] = safe_t[:, :-1] & safe_t[:, 1:]
                safe_v[1:ny, :] = safe_t[:-1, :] & safe_t[1:, :]

                u_geo = np.where(safe_u, u_geo, original_u)
                v_geo = np.where(safe_v, v_geo, original_v)

                destination_u[time_index, level, :, :] = u_geo
                destination_v[time_index, level, :, :] = v_geo

                dz_u, dz_v = thickness_to_faces(dz_t)
                hu += u_geo * dz_u
                hv += v_geo * dz_v

                # Retain only the fields needed by the diagnostic figures.
                u_section_faces = read_array(
                    source_u[
                        time_index,
                        level,
                        :,
                        section_x_index : section_x_index + 2,
                    ]
                )
                original_u_section[level, :] = np.mean(
                    u_section_faces,
                    axis=1,
                )
                geostrophic_u_section[level, :] = 0.5 * (
                    u_geo[:, section_x_index]
                    + u_geo[:, section_x_index + 1]
                )
                section_wet[level, :] = wet[:, section_x_index]

                if level == plot_level:
                    original_u = read_array(
                        source_u[time_index, level, :, :]
                    )
                    original_v = read_array(
                        source_v[time_index, level, :, :]
                    )
                    original_u_surface = 0.5 * (
                        original_u[:, :nx] + original_u[:, 1 : nx + 1]
                    )
                    original_v_surface = 0.5 * (
                        original_v[:ny, :] + original_v[1 : ny + 1, :]
                    )
                    geostrophic_u_surface = 0.5 * (
                        u_geo[:, :nx] + u_geo[:, 1 : nx + 1]
                    )
                    geostrophic_v_surface = 0.5 * (
                        v_geo[:ny, :] + v_geo[1 : ny + 1, :]
                    )

                previous_du_dz = du_dz
                previous_dv_dz = dv_dz
                previous_wet = wet

                if level == 0 or (level + 1) % 5 == 0 or level + 1 == nk:
                    log(f"Processed vertical level {level + 1}/{nk}")

            destination.sync()

        elapsed("density/geostrophic first pass", first_pass_start)

        hu *= umask
        hv *= vmask
        h_t[h_t < 0.0] = 0.0

        fx_w = hu[:, :nx] * metrics.dy_u[:, :nx]
        fx_e = hu[:, 1 : nx + 1] * metrics.dy_u[:, 1 : nx + 1]
        fy_s = hv[:ny, :] * metrics.dx_v[:ny, :]
        fy_n = hv[1 : ny + 1, :] * metrics.dx_v[1 : ny + 1, :]
        rhs_t = (fx_e - fx_w + fy_n - fy_s) / area
        rhs_t *= mask_t

        mask_eff = (mask_t == 1.0) & (h_t > 0.0)
        matrix, wet_flat, _ = build_poisson_matrix(
            h_t,
            mask_eff,
            area,
            metrics,
        )
        rhs = rhs_t.ravel(order="C")[wet_flat].copy()
        rhs -= np.mean(rhs)
        if apply_barotropic_correction:
            chi_vector, info, iterations, relative_residual = solve_poisson(
                matrix,
                rhs,
                rtol=cg_rtol,
                maxiter=cg_maxiter,
            )
            if info != 0 and require_convergence:
                raise RuntimeError(
                    "Conjugate-gradient solver did not converge: "
                    f"info={info}, "
                    f"relative_residual={relative_residual:.6e}"
                )
            if info != 0:
                log(
                    "[WARNING] CG did not converge within the requested "
                    "iterations. Continuing to match the MATLAB workflow, "
                    "which also writes the current iterate when pcg flag != 0."
                )
        else:
            log("[INFO] Skipping barotropic correction.")
            chi_vector = np.zeros(wet_flat.size, dtype=np.float64)
            info, iterations, relative_residual = 0, 0, 1.0

        chi_t = np.zeros((ny, nx), dtype=np.float64)
        chi_t.ravel(order="C")[wet_flat] = chi_vector

        dchi_u = np.zeros((ny, nx + 1), dtype=np.float64)
        dchi_v = np.zeros((ny + 1, nx), dtype=np.float64)
        dchi_u[:, 1:nx] = (
            chi_t[:, 1:nx] - chi_t[:, : nx - 1]
        ) * umask[:, 1:nx]
        dchi_v[1:ny, :] = (
            chi_t[1:ny, :] - chi_t[: ny - 1, :]
        ) * vmask[1:ny, :]

        uc = (dchi_u / metrics.dx_u) * umask
        vc = (dchi_v / metrics.dy_v) * vmask
        uc[~np.isfinite(uc)] = 0.0
        vc[~np.isfinite(vc)] = 0.0

        second_pass_start = perf_counter()
        with Dataset(partial_path, "r+") as destination:
            destination_u = destination.variables["u"]
            destination_v = destination.variables["v"]
            levels_to_correct = range(nk) if apply_barotropic_correction else ()
            for level in levels_to_correct:
                u_level = read_array(
                    destination_u[time_index, level, :, :]
                )
                v_level = read_array(
                    destination_v[time_index, level, :, :]
                )
                destination_u[time_index, level, :, :] = u_level + uc
                destination_v[time_index, level, :, :] = v_level + vc
                if level == 0 or (level + 1) % 10 == 0 or level + 1 == nk:
                    if apply_barotropic_correction:
                        log(f"Applied correction to level {level + 1}/{nk}")
            destination.sync()
        if apply_barotropic_correction:
            elapsed("apply/write barotropic correction", second_pass_start)

        he, hn = face_depths(h_t)
        hu_new = hu + he * dchi_u / metrics.dx_u
        hv_new = hv + hn * dchi_v / metrics.dy_v

        fx_w_new = hu_new[:, :nx] * metrics.dy_u[:, :nx]
        fx_e_new = hu_new[:, 1 : nx + 1] * metrics.dy_u[:, 1 : nx + 1]
        fy_s_new = hv_new[:ny, :] * metrics.dx_v[:ny, :]
        fy_n_new = hv_new[1 : ny + 1, :] * metrics.dx_v[1 : ny + 1, :]
        divergence_after = (
            fx_e_new - fx_w_new + fy_n_new - fy_s_new
        ) / area
        divergence_after *= mask_eff

        achi = matrix @ chi_vector
        divergence_theory = np.zeros((ny, nx), dtype=np.float64)
        divergence_theory.ravel(order="C")[wet_flat] = rhs - achi

        before_std = float(np.std(rhs_t[mask_eff], ddof=0))
        theory_std = float(np.std(divergence_theory[mask_eff], ddof=0))
        after_std = float(np.std(divergence_after[mask_eff], ddof=0))
        ratio = after_std / before_std
        log(
            "RMS: "
            f"before={before_std:.6e} "
            f"theory(after)={theory_std:.6e} "
            f"reconstructed(after)={after_std:.6e} "
            f"ratio={ratio:.6e}"
        )

        correction_u_t = 0.5 * (uc[:, :nx] + uc[:, 1 : nx + 1])
        correction_v_t = 0.5 * (vc[:ny, :] + vc[1 : ny + 1, :])
        adjusted_u_surface = geostrophic_u_surface + correction_u_t
        adjusted_v_surface = geostrophic_v_surface + correction_v_t
        adjusted_u_section = (
            geostrophic_u_section
            + correction_u_t[:, section_x_index][np.newaxis, :]
        )

        solver_diagnostics = SolverDiagnostics(
            wet_points=int(wet_flat.size),
            nonzeros=int(matrix.nnz),
            iterations=iterations,
            info=int(info),
            relative_residual=relative_residual,
            divergence_std_before=before_std,
            divergence_std_theory=theory_std,
            divergence_std_after=after_std,
            divergence_ratio=ratio,
        )

        plot_diagnostics = PlotDiagnostics(
            lon=lon,
            lat=lat,
            z=z,
            plot_level=plot_level,
            mask_t=mask_t,
            original_u_surface=original_u_surface,
            original_v_surface=original_v_surface,
            geostrophic_u_surface=geostrophic_u_surface,
            geostrophic_v_surface=geostrophic_v_surface,
            adjusted_u_surface=adjusted_u_surface,
            adjusted_v_surface=adjusted_v_surface,
            correction_u_t=correction_u_t,
            correction_v_t=correction_v_t,
            original_u_section=original_u_section,
            geostrophic_u_section=geostrophic_u_section,
            adjusted_u_section=adjusted_u_section,
            section_wet=section_wet,
            section_x_index=section_x_index,
            rhs_t=rhs_t,
            divergence_after=divergence_after,
        )

        if work_dir is None:
            os.replace(partial_path, output_path)
        else:
            move_start = perf_counter()
            publish_output(partial_path, output_path)
            elapsed("publish IC from work directory", move_start)
        log(f"Wrote adjusted IC: {output_path}")
        elapsed("TOTAL reconstruction", total_start)
        return solver_diagnostics, plot_diagnostics

    except BaseException:
        if partial_path.exists():
            partial_path.unlink()
        raise


def default_output_path(input_path: Path) -> Path:
    """Return ``*_geocurrents_python.nc`` beside the input file."""

    return input_path.with_name(
        input_path.stem + "_geocurrents_python.nc"
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct MOM6 geostrophic currents using the same numerical "
            "logic as geostrophic_adj/reconstruction_current.m"
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--grid-dir",
        required=True,
        type=Path,
        help="Directory containing ocean_hgrid.nc, ocean_mask.nc and topog.nc",
    )
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--time-index", type=int, default=0)
    parser.add_argument("--plot-level", type=int, default=0)
    parser.add_argument("--section-x-index", type=int)
    parser.add_argument(
        "--work-dir",
        type=Path,
        help=(
            "Optional fast local directory for the temporary 6+ GB IC copy; "
            "the completed file is moved to --output"
        ),
    )
    parser.add_argument("--plot-dir", type=Path)
    parser.add_argument(
        "--plot-format",
        choices=("png", "svg"),
        default="png",
        help=(
            "Diagnostic figure format; use svg to bypass HPC FreeType "
            "PNG glyph-rendering errors"
        ),
    )
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--cg-rtol", type=float, default=1.0e-8)
    parser.add_argument("--cg-maxiter", type=int, default=5000)
    parser.add_argument(
        "--require-convergence",
        action="store_true",
        help="Fail instead of matching MATLAB when CG reaches maxiter",
    )
    parser.add_argument(
        "--skip-barotropic-correction",
        action="store_true",
        help=(
            "Write reconstructed geostrophic currents without the "
            "depth-integrated Poisson correction"
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    """CLI entry point."""

    args = parse_arguments()
    input_path = args.input.resolve()
    output_path = (
        args.output.resolve()
        if args.output is not None
        else default_output_path(input_path)
    )
    grid_dir = args.grid_dir.resolve()
    plot_dir = (
        args.plot_dir.resolve()
        if args.plot_dir is not None
        else output_path.with_name(output_path.stem + "_diagnostics")
    )

    log(f"Input:     {input_path}")
    log(f"Output:    {output_path}")
    log(f"Grid dir:  {grid_dir}")
    if args.reference is not None:
        log(f"Reference: {args.reference.resolve()}")

    solver, plot_data = reconstruct(
        input_path=input_path,
        output_path=output_path,
        grid_dir=grid_dir,
        work_dir=(
            args.work_dir.resolve()
            if args.work_dir is not None
            else None
        ),
        time_index=args.time_index,
        plot_level=args.plot_level,
        section_x_index=args.section_x_index,
        cg_rtol=args.cg_rtol,
        cg_maxiter=args.cg_maxiter,
        require_convergence=args.require_convergence,
        apply_barotropic_correction=not args.skip_barotropic_correction,
        overwrite=args.overwrite,
    )

    summary: dict[str, Any] = {"solver": asdict(solver)}
    if not args.no_plots:
        make_diagnostic_plots(
            plot_data,
            plot_dir,
            plot_format=args.plot_format,
        )
        summary["plot_dir"] = str(plot_dir)
        summary["plot_format"] = args.plot_format

    if args.reference is not None:
        comparison_start = perf_counter()
        comparison = compare_uv_files(
            output_path,
            args.reference.resolve(),
            time_index=args.time_index,
        )
        summary["reference_comparison"] = comparison
        log("Reference comparison:")
        log(json.dumps(comparison, indent=2))
        elapsed("reference comparison", comparison_start)

#    summary_path = output_path.with_suffix(".summary.json")
#    summary_path.write_text(
#        json.dumps(summary, indent=2),
#        encoding="utf-8",
#    )
#    log(f"Wrote summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        log("Interrupted")
        raise SystemExit(130)
