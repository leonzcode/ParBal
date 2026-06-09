"""Build the two static support files ParBal Stage 1 needs for the Shasta run.

1. GLDAS_topo.h5  (coarse model terrain -> elevation correction)
   - /Grid/Z, /Grid/slope, /Grid/aspect on EXACTLY the subset weather grid
     (read_LDAS.m georeferences the weather stack with THIS file's
     ReferencingMatrix via load_coarse_topo -> makeInterp -> subsetGLDAS,
     so the grids must be identical: 4x6 cells, 0.25 deg).
   - Z comes from NASA's static GLDASp4_elevation_025d.nc4 (the terrain the
     GLDAS-2.1 model assumed); slope/aspect computed from it (near-zero at
     25 km cells; loaded by load_coarse_topo but unused in the GLDAS path).
   - Values stored as plain float32: load_coarse_topo h5-branch applies NO
     divisor scaling (unlike the fine topo files).

2. Shasta_landcover.h5  (forest map -> canopy corrections in makeFOREST)
   - /Grid/cc: canopy cover fraction 0-1 on the fine 169x166 Shasta grid,
     resampled from UCSB's MODIS-tile cc_h08v04.mat (sinusoidal, 463.3 m).
   - /Grid/Z:  veg-type fraction layers (deciduous, coniferous). Shasta is
     conifer country: deciduous=0, coniferous=cc, which reproduces what
     makeFOREST.m:44 assumes for unlabeled forest anyway.
   - include_vars_melt.m reads these datasets directly (no reprojection),
     so they must already match the fine topo grid exactly.

MATLAB/h5py orientation note: MATLAB h5read returns dimensions reversed
relative to h5py, so every array is written transposed from its MATLAB
orientation, and the (3,2) ReferencingMatrix is written as h5py (2,3).
"""
import glob

import h5py
import numpy as np
from pyproj import CRS, Transformer

DATA   = r"E:\ucsb\data\ParBal\Shasta"
TOPO   = DATA + r"\inputs\ShastaTopography.h5"
CCMAT  = DATA + r"\inputs\cc_h08v04.mat"
ELEV   = DATA + r"\inputs\GLDASp4_elevation_025d.nc4"
OUT_LT = DATA + r"\GLDAS\GLDAS_topo.h5"    # path in run_shasta_stage1.m; the dir
# must contain 'GLDAS' in its path: include_vars_melt.m:105-107 picks the LDAS
# variable list by matching NLDAS/GLDAS in the coarse-topo file's directory
OUT_LC = DATA + r"\inputs\Shasta_landcover.h5"

CELL_MODIS = 463.31271657  # m, MODIS sinusoidal 500 m grid
CELL_GLDAS = 0.25          # deg


# --------------------------------------------------------------------------
# 1. GLDAS coarse topo
# --------------------------------------------------------------------------
print("=== GLDAS_topo.h5 ===")

# the authoritative grid is whatever the downloaded weather files use
wx = sorted(glob.glob(DATA + r"\GLDAS\2019\091\*.nc4"))[0]
with h5py.File(wx, "r") as f:
    wlat = f["lat"][:].astype(np.float64)   # ascending (S->N)
    wlon = f["lon"][:].astype(np.float64)   # ascending (W->E)
print(f"weather grid: {wlat.size} lat x {wlon.size} lon, "
      f"lat {wlat[0]}..{wlat[-1]}, lon {wlon[0]}..{wlon[-1]}")

with h5py.File(ELEV, "r") as f:
    glat = f["lat"][:].astype(np.float64)   # global, ascending
    glon = f["lon"][:].astype(np.float64)
    gz   = f["GLDAS_elevation"][0, :, :]    # (lat, lon)

ilat = np.searchsorted(glat, wlat)
ilon = np.searchsorted(glon, wlon)
assert np.allclose(glat[ilat], wlat) and np.allclose(glon[ilon], wlon), \
    "weather grid does not sit on the GLDAS 0.25-deg lattice"

Z = gz[np.ix_(ilat, ilon)]                  # (4,6), row 0 = south
assert not np.any(Z < -999), "fill values in elevation subset"
Z = np.flipud(Z)                            # row 0 = north (read_LDAS rot90 convention)
lats_desc = wlat[::-1]

# slope/aspect of the coarse surface (loaded by load_coarse_topo, unused in
# the GLDAS path; computed anyway for completeness)
dy_m = CELL_GLDAS * 111_320.0                                  # N-S cell, m
dx_m = CELL_GLDAS * 111_320.0 * np.cos(np.deg2rad(lats_desc))  # E-W cell, m
dZdr, dZdc = np.gradient(Z.astype(np.float64))                 # per row/col
dzdy = -dZdr / dy_m                 # row index increases southward -> negate
dzdx = dZdc / dx_m[:, None]
slope = np.rad2deg(np.arctan(np.hypot(dzdx, dzdy)))
az_north_cw = np.rad2deg(np.arctan2(-dzdx, -dzdy)) % 360.0     # downslope, cw from N
aspect = 180.0 - az_north_cw                                   # ParBal: deg ccw from south

# ReferencingMatrix exactly as read_LDAS.m builds it for the weather stack:
# makerefmat(lon(1), lat(end), dlon, -dlat) -> MATLAB (3,2) -> h5py (2,3)
refmat = np.array([[0.0,        CELL_GLDAS, wlon[0] - CELL_GLDAS],
                   [-CELL_GLDAS, 0.0,       wlat[-1] + CELL_GLDAS]])

with h5py.File(OUT_LT, "w") as f:
    g = f.create_group("Grid")
    for name, arr in (("Z", Z), ("slope", slope), ("aspect", aspect)):
        g.create_dataset(name, data=arr.astype(np.float32).T)  # transpose for MATLAB
    g.attrs["ReferencingMatrix"] = refmat
    g.attrs["mapprojection"] = "geographic"
print(f"Z (row 0 = north):\n{np.round(Z).astype(int)}")
print(f"slope max {slope.max():.2f} deg")
print(f"wrote {OUT_LT}")


# --------------------------------------------------------------------------
# 2. Shasta land cover
# --------------------------------------------------------------------------
print("\n=== Shasta_landcover.h5 ===")

with h5py.File(TOPO, "r") as f:
    rm = f["Grid"].attrs["ReferencingMatrix"]       # h5py (2,3)
    ncols, nrows = f["Grid/elevation"].shape        # h5py reversed -> MATLAB (169,166)
    dem = f["Grid/elevation"][()].T.astype(np.float64) / f["Grid/elevation"].attrs["divisor"][0]
print(f"fine grid (MATLAB rows x cols): {nrows} x {ncols}")

# pixel centers from the refmat: x = R[0,2]+R[0,1]*col, y = R[1,2]+R[1,0]*row
cols = np.arange(1, ncols + 1)
rows = np.arange(1, nrows + 1)
x = rm[0, 2] + rm[0, 1] * cols
y = rm[1, 2] + rm[1, 0] * rows
X, Y = np.meshgrid(x, y)                            # (nrows, ncols), MATLAB orientation

# CA Teale Albers (the topo attrs: eqaconicstd, parallels 34/40.5, lon0 -120,
# falsenorthing -4e6, WGS84 ellipsoid) -> MODIS sinusoidal (R=6371007.181)
teale = CRS.from_proj4("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 "
                       "+x_0=0 +y_0=-4000000 +ellps=WGS84 +units=m +no_defs")
sinu = CRS.from_proj4("+proj=sinu +R=6371007.181 +x_0=0 +y_0=0 +units=m +no_defs")
xs, ys = Transformer.from_crs(teale, sinu, always_xy=True).transform(X, Y)

with h5py.File(CCMAT, "r") as f:                    # v7.3 .mat = HDF5
    cc_tile = np.nan_to_num(f["cc"][()].T)          # transpose -> MATLAB (row,col), row 0 = north
    X0, Y0 = f["#refs#/e/TiePointWorld"][0]         # upper-left corner, world m

# fractional 0-based indices into the tile, then bilinear sample
cf = (xs - X0) / CELL_MODIS - 0.5
rf = (Y0 - ys) / CELL_MODIS - 0.5
assert cf.min() > 0 and rf.min() > 0 and cf.max() < 2399 and rf.max() < 2399, \
    "Shasta window falls outside tile h08v04"
r0 = np.floor(rf).astype(int); c0 = np.floor(cf).astype(int)
fr = rf - r0;                   fc = cf - c0
cc = (cc_tile[r0,     c0    ] * (1 - fr) * (1 - fc) +
      cc_tile[r0,     c0 + 1] * (1 - fr) * fc +
      cc_tile[r0 + 1, c0    ] * fr       * (1 - fc) +
      cc_tile[r0 + 1, c0 + 1] * fr       * fc).astype(np.float32)

Zveg = np.stack([np.zeros_like(cc), cc])            # (2, nrows, ncols): deciduous, coniferous

with h5py.File(OUT_LC, "w") as f:
    g = f.create_group("Grid")
    g.create_dataset("cc", data=cc.T)               # -> MATLAB (169,166)
    g.create_dataset("Z", data=Zveg.transpose(0, 2, 1))  # -> MATLAB (169,166,2)
    with h5py.File(TOPO, "r") as t:                 # carry the georeference along
        for k, v in t["Grid"].attrs.items():
            g.attrs[k] = v
print(f"cc range {cc.min():.3f}..{cc.max():.3f}, mean {cc.mean():.3f}, "
      f"forested fraction {(cc > 0).mean():.2f}")
print(f"wrote {OUT_LC}")

# Zdiff sanity: fine DEM vs the coarse cells it falls in
print(f"\nsanity: fine DEM mean {dem.mean():.0f} m (range {dem.min():.0f}.."
      f"{dem.max():.0f}), coarse Z range {Z.min():.0f}..{Z.max():.0f} m")
