#!/usr/bin/env python3
"""Build lola_horizons.json (terrain horizon masks) from a NASA LOLA south-polar DEM.

Usage:
  pip install rasterio pyproj numpy
  python make_lola_horizons.py --dem LDEM_80S_20MPP_ADJ.TIF --html ../index.html --out ../lola_horizons.json

DEM: LRO LOLA south-pole digital elevation model (NASA PGDA, https://pgda.gsfc.nasa.gov/products/90).
Only sites inside the DEM (south of about 80 S) get a mask; others are skipped and the web page
keeps a smooth-sphere horizon for them. Rays stop at the DEM edge, so distant terrain beyond it is not seen.
"""
import argparse, json, math, re
import numpy as np, rasterio
from rasterio.windows import Window
from pyproj import Transformer

R = 1737400.0  # lunar reference radius (m), same sphere Horizons uses

def parse_sites(path):
    src = open(path, encoding="utf-8").read()
    return [(m[0], m[1].replace("\\'", "'"), float(m[2]), float(m[3]))
            for m in re.findall(r"\{id:'(\w+)',name:'((?:[^'\\]|\\.)*)',lat:(-?[\d.]+),lon:(-?[\d.]+)", src)]

def horizon(ds, tf, lat, lon, naz, max_m, scale):
    res = abs(ds.res[0])                      # metres per pixel (polar stereographic, metres)
    x0, y0 = tf.transform(lon, lat)
    col, row = ds.index(x0, y0)[1], ds.index(x0, y0)[0]
    if not (0 <= row < ds.height and 0 <= col < ds.width): return None
    half = int(max_m / res) + 2
    c0, r0 = max(0, col - half), max(0, row - half)
    w = Window(c0, r0, min(ds.width, col + half) - c0, min(ds.height, row + half) - r0)
    dem = ds.read(1, window=w).astype("float64") * scale
    nod = ds.nodata
    if nod is not None: dem[dem == nod * scale] = np.nan
    def height(xs, ys):
        rr, cc = rasterio.transform.rowcol(ds.transform, xs, ys)
        rr = np.asarray(rr) - r0; cc = np.asarray(cc) - c0
        ok = (rr >= 0) & (rr < dem.shape[0]) & (cc >= 0) & (cc < dem.shape[1])
        out = np.full(rr.shape, np.nan); out[ok] = dem[rr[ok], cc[ok]]
        return out
    h0 = float(height(np.array([x0]), np.array([y0]))[0])
    if not np.isfinite(h0): return None
    d = np.arange(res, max_m, res); dl = d / R
    p, l = math.radians(lat), math.radians(lon); out = []
    for k in range(naz):
        a = 2 * math.pi * k / naz
        p2 = np.arcsin(math.sin(p) * np.cos(dl) + math.cos(p) * np.sin(dl) * math.cos(a))
        l2 = l + np.arctan2(math.sin(a) * np.sin(dl) * math.cos(p), np.cos(dl) - math.sin(p) * np.sin(p2))
        xs, ys = tf.transform(np.degrees(l2), np.degrees(p2))
        h = height(xs, ys)
        bad = ~np.isfinite(h)
        if bad.any(): h, dl_k = h[:np.argmax(bad)], dl[:np.argmax(bad)]
        else: dl_k = dl
        if h.size == 0: out.append(0.0); continue
        ang = np.arctan2((R + h) * np.cos(dl_k) - (R + h0), (R + h) * np.sin(dl_k))
        out.append(round(math.degrees(float(ang.max())), 3))
    return out, h0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", required=True); ap.add_argument("--html", default="../index.html")
    ap.add_argument("--out", default="../lola_horizons.json"); ap.add_argument("--azimuths", type=int, default=360)
    ap.add_argument("--max-km", type=float, default=60.0, help="ray length")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply DEM values to get metres")
    a = ap.parse_args()
    ds = rasterio.open(a.dem)
    tf = Transformer.from_crs(ds.crs.geodetic_crs, ds.crs, always_xy=True)
    sites = {}
    for sid, name, lat, lon in parse_sites(a.html):
        r = horizon(ds, tf, lat, lon, a.azimuths, a.max_km * 1000, a.scale) if lat <= -79.9 else None
        if r is None: print(f"skip {sid} {name}: outside DEM"); continue
        hz, h0 = r; sites[sid] = {"name": name, "lat": lat, "lon": lon, "elev_m": round(h0, 1), "hz": hz}
        print(f"{sid} {name}: max horizon {max(hz):.2f} deg")
    json.dump({"source": "NASA LRO LOLA south-polar DEM", "dem": a.dem, "max_km": a.max_km, "azimuths": a.azimuths, "sites": sites},
              open(a.out, "w"), separators=(",", ":"))
    print("wrote", a.out)
if __name__ == "__main__": main()
