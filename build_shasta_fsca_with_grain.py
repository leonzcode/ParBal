"""Merge grain_size + dust from the WUS SPIRES mosaic into the Shasta fSCA cube.

Why: UCSB's published Shasta_2019.h5 contains ONLY snow_fraction (2.5 MB);
ParBal's dailyEnergy needs grain_size (um) and dust (ppm) for the snow albedo
model. Their WUS-wide mosaic2019.h5 carries all variables on the "Western USA
Albers" 500 m grid, so we resample its grain/dust onto the Shasta grid
(nearest neighbor - these are per-pixel snow properties, no smearing wanted).

Keeps the OFFICIAL Shasta snow_fraction untouched (same vintage as the
official reconstruction we validate against); only grain_size and dust come
from the newer WUS mosaic (the 2020-vintage grain/dust were never published).

Output: a new fsca_Shasta_2019.h5 with all three variables; the original
fsca-only cube is kept as fsca_Shasta_2019_fscaonly_orig.h5.

Orientation reminder: MATLAB h5 arrays appear dimension-reversed in h5py;
MATLAB (row,col,day) == h5py (day,col,row).
"""
import os
import shutil

import h5py
import numpy as np
from pyproj import CRS, Transformer
from scipy import ndimage

FILLV = 65535  # uint16 fill in the WUS mosaic (no snow / no retrieval)

INPUTS = r"D:\code\ucsb\data\ParBal\Shasta\inputs"
MOSAIC = INPUTS + r"\wus_mosaic2019.h5"
CUBE   = INPUTS + r"\fsca_Shasta_2019.h5"
ORIG   = INPUTS + r"\fsca_Shasta_2019_fscaonly_orig.h5"
GROUP  = "Grid/MODIS_GRID_500m"


def albers_from_attrs(attrs):
    """Build a pyproj CRS from ParBal's eqaconicstd /Grid attributes."""
    assert attrs["mapprojection"] in (b"eqaconicstd", "eqaconicstd")
    a, ecc = attrs["geoid"].ravel()
    b = a * np.sqrt(1 - ecc**2)
    sp = attrs["mapparallels"].ravel()
    lat0, lon0 = attrs["origin"].ravel()[:2]
    return CRS.from_proj4(
        f"+proj=aea +lat_1={sp[0]} +lat_2={sp[-1]} +lat_0={lat0} +lon_0={lon0} "
        f"+x_0={attrs['falseeasting'].ravel()[0]} +y_0={attrs['falsenorthing'].ravel()[0]} "
        f"+a={a} +b={b} +units=m +no_defs")


# --- the official Shasta cube: grid + dates + snow_fraction stay as-is -----
SRC = ORIG if os.path.exists(ORIG) else CUBE           # rerun-safe
with h5py.File(SRC, "r") as f:
    sh_rm = f[GROUP].attrs["ReferencingMatrix"]        # h5py (2,3)
    sh_dates = f.attrs["MATLABdates"]
    sh_crs = albers_from_attrs(dict(f["Grid"].attrs))
    sf_all = f[GROUP + "/snow_fraction"][()].transpose(0, 2, 1)  # (365,row,col)
    nrows, ncols = sf_all.shape[1:]                    # MATLAB order

# Shasta pixel centers in Teale Albers (MATLAB orientation (row,col))
cols = np.arange(1, ncols + 1)
rows = np.arange(1, nrows + 1)
X, Y = np.meshgrid(sh_rm[0, 2] + sh_rm[0, 1] * cols,
                   sh_rm[1, 2] + sh_rm[1, 0] * rows)

# --- the WUS mosaic: locate those centers on its grid ----------------------
with h5py.File(MOSAIC, "r") as f:
    wus_rm = f[GROUP].attrs["ReferencingMatrix"]
    wus_dates = f.attrs["MATLABdates"]
    wus_crs = albers_from_attrs(dict(f["Grid"].attrs))
    assert np.array_equal(np.asarray(sh_dates), np.asarray(wus_dates)), \
        "mosaic dates differ from Shasta cube dates"

    xw, yw = Transformer.from_crs(sh_crs, wus_crs, always_xy=True).transform(X, Y)
    colw = np.rint((xw - wus_rm[0, 2]) / wus_rm[0, 1]).astype(int)   # MATLAB 1-based
    roww = np.rint((yw - wus_rm[1, 2]) / wus_rm[1, 0]).astype(int)
    print(f"Shasta window in mosaic: rows {roww.min()}..{roww.max()}, "
          f"cols {colw.min()}..{colw.max()}")

    out = {}
    r0, r1 = roww.min(), roww.max()
    c0, c1 = colw.min(), colw.max()
    ri = roww - r0      # 0-based into the window
    ci = colw - c0
    for name in ("grain_size", "dust", "snow_fraction"):
        ds = f[f"{GROUP}/{name}"]
        # h5py [day, col-1, row-1]; read the bounding window once
        W = ds[:, c0 - 1:c1, r0 - 1:r1]                # (365, wc, wr)
        cube = W[:, ci, ri]                            # (365, nrows, ncols), MATLAB (row,col) per day
        if name == "snow_fraction":
            wus_sf = cube                              # only for the alignment check
        else:
            out[name] = (cube, ds.dtype, dict(ds.attrs))

# --- alignment sanity: official fsca vs mosaic fsca on the same grid -------
off_sf = sf_all[196].astype(float)                     # Apr 15
m = wus_sf[196].astype(float)
m[m == 255] = 0                                        # uint8 fill
both = (off_sf > 0) | (m > 0)
corr = np.corrcoef(off_sf[both], m[both])[0, 1]
print(f"Apr 15 fsca: official snow px {(off_sf>0).sum()}, mosaic snow px {(m>0).sum()}, "
      f"corr {corr:.3f} (geometry is wrong if this is low)")

# --- fills + gap-fill -------------------------------------------------------
# Mosaic grain/dust are FILLV where the 2025-vintage SPIRES saw no snow. The
# official (2020-vintage) fsca is more inclusive, so some official-snow pixels
# lack grain/dust. Fill them per day from the nearest valid snow pixel
# (snow properties vary smoothly); for days where the mosaic has no valid
# pixels at all, carry the previous day's fields (snow evolves slowly),
# with a backward pass for any leading days.
grain = out["grain_size"][0]
dust = out["dust"][0]
grain[grain == FILLV] = 0
dust[dust == FILLV] = 0
pending, prev = [], None
filled_px = 0
for dd in range(grain.shape[0]):
    g, du = grain[dd], dust[dd]
    need = (sf_all[dd] > 0) & (g == 0)
    valid = g > 0
    if valid.any():
        idx = tuple(ndimage.distance_transform_edt(
            ~valid, return_distances=False, return_indices=True))
        g[need] = g[idx][need]
        du[need] = du[idx][need]
        prev = dd
        for pdd in pending:                            # backward-fill leading days
            pneed = (sf_all[pdd] > 0) & (grain[pdd] == 0)
            grain[pdd][pneed] = g[pneed]
            dust[pdd][pneed] = du[pneed]
        pending = []
    elif need.any():
        if prev is not None:                           # carry forward
            g[need] = grain[prev][need]
            du[need] = dust[prev][need]
        else:
            pending.append(dd)
    filled_px += int(need.sum())
print(f"gap-filled grain/dust on {filled_px} px-days "
      f"(official snow without 2025-vintage retrieval)")

g = grain[196] / out["grain_size"][2]["divisor"][0]
d = dust[196] / out["dust"][2]["divisor"][0]
print(f"Apr 15 grain: {g[g>0].min():.0f}..{g.max():.0f} um on {(g>0).sum()} px; "
      f"dust: {d.max():.1f} ppm max")

# --- write the merged cube --------------------------------------------------
if not os.path.exists(ORIG):
    shutil.move(CUBE, ORIG)
shutil.copy(ORIG, CUBE)                                # start from the original
with h5py.File(CUBE, "a") as f:
    for name, (cube, dtype, attrs) in out.items():
        ds = f.create_dataset(f"{GROUP}/{name}",
                              data=cube.transpose(0, 2, 1),  # -> MATLAB (row,col,day)
                              dtype=dtype, compression="gzip", compression_opts=1)
        for k, v in attrs.items():
            ds.attrs[k] = v
print(f"wrote {CUBE} (original kept as {os.path.basename(ORIG)})")
