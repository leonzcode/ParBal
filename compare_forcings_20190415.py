"""Compare our Stage-1 hourly downscaled forcings vs UCSB's, Apr 15 2019.

Ours:   energy\20190415.mat (metvars_flag=true run; MATLAB v7.3 = HDF5)
Theirs: ucsb_stage1\20190415.h5 (extracted from ParBalForcingsWY2001to2019.zip)

Both are hourly (24 x grid). Goal: attribute the -27% melt bias to either
radiation source (GLDAS vs CERES) or snow albedo (2025- vs 2020-vintage
grain/dust) or temperature/wind.
"""
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OURS = r"E:\ucsb\data\ParBal\Shasta\energy\20190415.mat"
UCSB = r"E:\ucsb\data\ParBal\Shasta\ucsb_stage1\20190415.h5"
OUT = r"E:\ucsb\data\ParBal\Shasta\outputs\compare_forcings_20190415.png"

def ours_var(f, name, fill, scale=1.0, offset=0.0):
    x = f[name][()].astype(np.float32)          # h5py (24,166,169)
    x[x == fill] = np.nan
    return x * scale + offset                   # (24,166,169); transpose per-hour at use

def ucsb_var(f, name):
    ds = f["Grid/" + name]
    div = float(np.ravel(ds.attrs["divisor"])[0]) if "divisor" in ds.attrs else 1.0
    x = ds[()].astype(np.float32) / div
    return x

with h5py.File(OURS, "r") as f:
    ours = {
        "albedo":  ours_var(f, "albedo", 255),            # %
        "Lin":     ours_var(f, "Lin", 65535),             # W/m2
        "diffuse": ours_var(f, "diffuse", 65535),         # W/m2
        "direct":  ours_var(f, "direct", 65535),          # W/m2 (slope-incident)
        "Ta":      ours_var(f, "Ta", -128, 1.0, 273.15),  # int8 degC -> K
        "pres":    ours_var(f, "presZ", 255, 10.0),       # uint8 kPa -> mb
        "wind":    ours_var(f, "windspd", 127),           # m/s (fill=127, intmax int8)
    }
with h5py.File(UCSB, "r") as f:
    theirs = {
        "albedo":  ucsb_var(f, "SnowAlbedo"),
        "Lin":     ucsb_var(f, "DownwellingLongwaveRadiation"),
        "diffuse": ucsb_var(f, "DiffuseSolarRadiation"),
        "direct":  ucsb_var(f, "DirectNormalSolarRadiation"),  # beam-normal! see note
        "Ta":      ucsb_var(f, "AirTemperature"),
        "pres":    ucsb_var(f, "AirPressure"),
        "wind":    ucsb_var(f, "WindSpeed"),
    }

# common snow mask per hour: both albedos defined and positive
snow = (ours["albedo"] > 0) & (theirs["albedo"] > 0)
print(f"snow px in common (noon): {snow[12].sum()}")

print(f"\n{'variable':28s} {'ours':>9s} {'UCSB':>9s} {'bias':>8s} {'corr@noon':>9s}")
for k, label, day in [
    ("albedo", "snow albedo (%)", False),
    ("Lin", "longwave down (W/m2)", True),
    ("diffuse", "diffuse solar (W/m2)", True),
    ("direct", "direct solar* (W/m2)", True),
    ("Ta", "air temperature (K)", True),
    ("pres", "pressure (mb)", True),
    ("wind", "wind speed (m/s)", True),
]:
    o, t = ours[k], theirs[k]
    m = snow & ~np.isnan(o) & ~np.isnan(t)
    if k in ("diffuse", "direct"):                # daylight hours only
        m &= (t > 0) | (o > 0)
    om, tm = np.nanmean(o[m]), np.nanmean(t[m])
    h = 12                                        # local noon
    mh = m[h]
    corr = np.corrcoef(o[h][mh], t[h][mh])[0, 1] if mh.sum() > 10 else np.nan
    print(f"{label:28s} {om:9.2f} {tm:9.2f} {om-tm:+8.2f} {corr:9.3f}")
print("* their direct is beam-NORMAL, ours is slope-incident: bias not meaningful, corr is")

# albedo maps at noon
fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.2))
o12 = np.where(snow[12], ours["albedo"][12], np.nan).T   # -> (169,166) north-up
t12 = np.where(snow[12], theirs["albedo"][12], np.nan).T
d12 = o12 - t12
for a, d, t, cm, vr in [
    (ax[0], t12, "UCSB albedo (noon, %)", "viridis", (40, 90)),
    (ax[1], o12, "our albedo (noon, %)", "viridis", (40, 90)),
    (ax[2], d12, "ours - UCSB (%)", "RdBu_r", (-25, 25)),
]:
    im = a.imshow(d, cmap=cm, vmin=vr[0], vmax=vr[1])
    im.cmap.set_bad("0.88")
    a.set_title(t); a.set_xticks([]); a.set_yticks([])
    fig.colorbar(im, ax=a, shrink=0.8)
fig.suptitle("Snow albedo, 15 Apr 2019 noon - 2025-vintage grain/dust (ours) vs 2020 (UCSB)")
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print(f"\nsaved {OUT}")
