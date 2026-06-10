"""Final pipeline validation: our SWE reconstruction vs UCSB's official, WY2019.

Compares my_reconstruction_Shasta_2019.h5 (our full pipeline: GLDAS-only
forcings, rebuilt land cover / coarse topo, 2025-vintage grain+dust) against
reconstruction_Shasta_2019.h5 (UCSB official). April 1 headline map + melt-
season time series + per-day stats. SWE is only valid peak->melt-out, so
judgement window is Apr-Jul.
"""
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OURS = r"D:\code\ucsb\data\ParBal\Shasta\outputs\my_reconstruction_Shasta_2019_nocc.h5"
# no-canopy = the validated config (A/B: corr 0.978, bias -3.0% vs official;
# the canopy-ON file my_reconstruction_Shasta_2019.h5 came out +26.2%)
OFF = r"D:\code\ucsb\data\ParBal\Shasta\reconstruction_Shasta_2019.h5"
OUT = r"D:\code\ucsb\data\ParBal\Shasta\outputs\compare_swe_final.png"
FILL = 65535
APR1 = 182  # 0-based day index (Oct 1 = 0)

def load(path, var):
    with h5py.File(path, "r") as f:
        x = f["Grid/" + var][()].astype(np.float32)   # (365,166,169) h5py
    x[x == FILL] = np.nan
    return x

swe_o = load(OURS, "swe")
swe_u = load(OFF, "swe")

# April 1 maps (transpose -> north-up MATLAB orientation)
a, b = swe_o[APR1].T, swe_u[APR1].T
v = ~np.isnan(a) & ~np.isnan(b) & ((a > 0) | (b > 0))
corr = np.corrcoef(b[v], a[v])[0, 1]
bias = np.mean(a[v] - b[v])
rmse = np.sqrt(np.mean((a[v] - b[v]) ** 2))
print(f"April 1 2019 SWE: {v.sum()} px | official mean {np.mean(b[v]):.0f} max {np.nanmax(b):.0f} mm "
      f"| ours mean {np.mean(a[v]):.0f} max {np.nanmax(a):.0f} mm")
print(f"corr {corr:.3f}, bias {bias:+.0f} mm, RMSE {rmse:.0f} mm")

# domain-total SWE time series (km3 of water), Apr-Jul validity window
px_km2 = (463.3127 / 1000) ** 2
tot_o = np.nansum(np.nan_to_num(swe_o), axis=(1, 2)) * px_km2 / 1e6  # mm*km2 -> km3
tot_u = np.nansum(np.nan_to_num(swe_u), axis=(1, 2)) * px_km2 / 1e6
days = np.arange(365)

# per-day spatial correlation, Apr-Jul
corrs = []
for d in range(182, 304):  # Apr 1 .. Jul 31
    x, y = swe_o[d], swe_u[d]
    m = ~np.isnan(x) & ~np.isnan(y) & ((x > 0) | (y > 0))
    corrs.append(np.corrcoef(y[m], x[m])[0, 1] if m.sum() > 100 else np.nan)
corrs = np.array(corrs)
print(f"per-day spatial corr Apr-Jul: median {np.nanmedian(corrs):.3f}, "
      f"range {np.nanmin(corrs):.3f}..{np.nanmax(corrs):.3f}")

fig = plt.figure(figsize=(16, 9))
gs = fig.add_gridspec(2, 3, height_ratios=[1.4, 1])
vmax = np.nanpercentile(b[v], 99)
for i, (d, t) in enumerate([(b, "UCSB official SWE"), (a, "Leon's reconstruction")]):
    ax = fig.add_subplot(gs[0, i])
    im = ax.imshow(d, cmap="Blues", vmin=0, vmax=vmax)
    im.cmap.set_bad("0.88")
    ax.set_title(f"{t} - 1 Apr 2019 (mm)")
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, shrink=0.75)
ax = fig.add_subplot(gs[0, 2])
ax.plot(b[v], a[v], ".", ms=2, alpha=0.25)
lim = [0, float(vmax) * 1.4]
ax.plot(lim, lim, "r-", lw=1)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("official (mm)"); ax.set_ylabel("Leon's (mm)")
ax.set_title(f"Apr 1 pixel-by-pixel (r={corr:.2f})")

ax = fig.add_subplot(gs[1, :2])
ax.plot(days, tot_u, label="UCSB official", lw=2)
ax.plot(days, tot_o, label="Leon's (GLDAS-only)", lw=2)
ax.axvspan(182, 303, color="gold", alpha=0.15, label="validity window (Apr-Jul)")
ax.set_xticks([0, 31, 61, 92, 123, 151, 182, 212, 243, 273, 304, 334])
ax.set_xticklabels(["O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S"])
ax.set_ylabel("domain SWE (km$^3$ water)")
ax.set_title("Total snow water storage, WY2019")
ax.legend()
ax = fig.add_subplot(gs[1, 2])
ax.plot(np.arange(182, 304), corrs)
ax.set_ylim(0, 1)
ax.set_xticks([182, 212, 243, 273, 303])
ax.set_xticklabels(["Apr", "May", "Jun", "Jul", "Aug"])
ax.set_ylabel("spatial corr")
ax.set_title("Per-day map correlation")
fig.suptitle("ParBal Shasta WY2019: Leon's full pipeline vs UCSB official reconstruction", fontsize=14)
fig.tight_layout()
fig.savefig(OUT, dpi=120)
print(f"saved {OUT}")
