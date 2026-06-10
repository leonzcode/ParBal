"""Stage-1 melt validation: Leon's run vs UCSB official, April 15 2019.

Converts Leon's Stage-1 melt energy M to mm exactly as reconstructSWE does
(mf = 3600/(rho_water*Lf)*1000, scaled by fsca) and compares to the official
melt cube for the same day.
"""
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DAY = 196  # 0-based index: Apr 15 2019 (Oct 1 = 0)

with h5py.File(r"D:\code\ucsb\data\ParBal\Shasta\energy\20190415.mat", "r") as f:
    M = f["M"][()].T.astype(np.float32)            # -> MATLAB (169,166)
M[M == 65535] = np.nan
mf = 3600.0 / (1000 * 3.34e5) * 1000
with h5py.File(r"D:\code\ucsb\data\ParBal\Shasta\inputs\fsca_Shasta_2019.h5", "r") as f:
    fsca = f["Grid/MODIS_GRID_500m/snow_fraction"][DAY, :, :].T.astype(np.float32) / 100.0
ours = M * mf * fsca

with h5py.File(r"D:\code\ucsb\data\ParBal\Shasta\reconstruction_Shasta_2019.h5", "r") as f:
    off = f["Grid/melt"][DAY, :, :].astype(np.float32).T
off[off == 65535] = np.nan

v = ~np.isnan(off) & ~np.isnan(ours) & ((off > 0) | (ours > 0))
print(f"pixels compared: {v.sum()}")
print(f"official: mean {np.nanmean(off[v]):.2f} max {np.nanmax(off):.1f} mm")
print(f"Leon's  : mean {np.nanmean(ours[v]):.2f} max {np.nanmax(ours):.1f} mm")
print(f"bias {np.mean(ours[v]-off[v]):+.2f} mm, RMSE {np.sqrt(np.mean((ours[v]-off[v])**2)):.2f} mm, "
      f"corr {np.corrcoef(off[v],ours[v])[0,1]:.3f}")

fig, ax = plt.subplots(1, 3, figsize=(16, 5.5))
vmax = max(np.nanpercentile(off[v], 99), np.nanpercentile(ours[v], 99))
for a, d, t in [(ax[0], off, "UCSB official melt"),
                (ax[1], ours, "Leon's Stage-1 melt (real GLDAS weather)")]:
    im = a.imshow(d, cmap="viridis", vmin=0, vmax=vmax)
    im.cmap.set_bad("0.88")
    a.set_title(t + " - 15 Apr 2019 (mm)")
    a.set_xticks([]); a.set_yticks([])
    fig.colorbar(im, ax=a, shrink=0.8)
ax[2].plot(off[v], ours[v], ".", ms=2, alpha=0.3)
lim = [0, float(vmax) * 1.3]
ax[2].plot(lim, lim, "r-", lw=1)
ax[2].set_xlim(lim); ax[2].set_ylim(lim)
ax[2].set_xlabel("official (mm)")
ax[2].set_ylabel("Leon's (mm)")
ax[2].set_title("pixel-by-pixel")
fig.tight_layout()
fig.savefig(r"D:\code\ucsb\data\ParBal\Shasta\outputs\compare_stage1_20190415.png", dpi=130)
print("saved compare_stage1_20190415.png")
