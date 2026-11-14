#!/usr/bin/env python3
"""Download Sun and Earth azimuth/elevation for every site in index.html from NASA JPL Horizons
and save them as horizons_data.json (hourly steps). Standard library only.

Browsers cannot call the Horizons API directly (it sends no CORS headers), so the page reads this file instead.
Usage:  python tools/fetch_horizons.py --html index.html --out horizons_data.json [--start 2026-11-01 --end 2028-01-01]
"""
import argparse, json, re, sys, time, urllib.parse, urllib.request
from datetime import datetime, timedelta, timezone

URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
MON = {m: i for i, m in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split())}

def parse_sites(path):
    src = open(path, encoding="utf-8").read()
    return [(m[0], m[1].replace("\\'", "'"), float(m[2]), float(m[3]))
            for m in re.findall(r"\{id:'(\w+)',name:'((?:[^'\\]|\\.)*)',lat:(-?[\d.]+),lon:(-?[\d.]+)", src)]

def parse(raw):
    i, j = raw.find("$$SOE"), raw.find("$$EOE")
    if i < 0 or j < 0: raise ValueError("no ephemeris rows: " + raw[:300].replace("\n", " "))
    def cells(line):
        c = [x.strip() for x in line.split(",")]
        while c and not c[-1]: c.pop()
        return c
    heads = [l for l in raw[:i].split("\n") if re.search("Azi", l) and re.search("Elev", l)]
    ai, ei = -2, -1
    if heads:
        h = cells(heads[-1])
        a = next((k for k, x in enumerate(h) if re.match(r"(?i)azi", x)), -1)
        e = next((k for k, x in enumerate(h) if re.match(r"(?i)elev", x)), -1)
        if a >= 0 and e >= 0: ai, ei = a, e
    out = {}
    for line in raw[i + 5:j].split("\n"):
        if not line.strip(): continue
        c = cells(line)
        m = re.match(r"(\d{4})-([A-Za-z]{3})-(\d{2}) (\d{2}):(\d{2})", c[0])
        if not m: continue
        try:
            az = float(c[ai if ai >= 0 else len(c) + ai]); el = float(c[ei if ei >= 0 else len(c) + ei])
        except (ValueError, IndexError): continue
        t = datetime(int(m[1]), MON[m[2]] + 1, int(m[3]), int(m[4]), int(m[5]), tzinfo=timezone.utc)
        out[t] = (az % 360.0, el)
    if not out: raise ValueError("could not read the Horizons rows")
    return out

def query(lon, lat, target, t0, t1):
    p = {"format": "json", "COMMAND": f"'{target}'", "OBJ_DATA": "NO", "MAKE_EPHEM": "YES", "EPHEM_TYPE": "OBSERVER",
         "CENTER": "coord@301", "COORD_TYPE": "GEODETIC", "SITE_COORD": f"{lon % 360:.6f},{lat:.6f},0",
         "START_TIME": f"'{t0:%Y-%m-%d %H:%M}'", "STOP_TIME": f"'{t1:%Y-%m-%d %H:%M}'", "STEP_SIZE": "'60 m'",
         "QUANTITIES": "'4'", "CSV_FORMAT": "YES", "ANG_FORMAT": "DEG", "EXTRA_PREC": "YES"}
    for attempt in range(4):
        try:
            with urllib.request.urlopen(URL + "?" + urllib.parse.urlencode(p), timeout=120) as r:
                j = json.load(r)
            if "error" in j: raise ValueError(j["error"])
            return parse(j.get("result", ""))
        except Exception as e:
            print(f"   retry {attempt + 1}: {e}", file=sys.stderr); time.sleep(5 * (attempt + 1))
    raise SystemExit("JPL Horizons request failed repeatedly")

def series(lon, lat, target, start, n, chunk_h):
    got = {}
    for k in range(0, n, chunk_h):
        t0 = start + timedelta(hours=k); t1 = start + timedelta(hours=min(n, k + chunk_h) - 1)
        got.update(query(lon, lat, target, t0, t1)); time.sleep(0.4)
    flat = []
    for k in range(n):
        v = got.get(start + timedelta(hours=k))
        if v is None: raise SystemExit(f"missing Horizons row at {start + timedelta(hours=k)}")
        flat += [round(v[0], 2), round(v[1], 2)]
    return flat

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", default="index.html"); ap.add_argument("--out", default="horizons_data.json")
    ap.add_argument("--start", default="2026-11-01"); ap.add_argument("--end", default="2028-01-01")
    ap.add_argument("--chunk-days", type=int, default=60)
    a = ap.parse_args()
    start = datetime.strptime(a.start, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    end = datetime.strptime(a.end, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    n = int((end - start).total_seconds() // 3600); chunk = a.chunk_days * 24
    sites = {}
    for sid, name, lat, lon in parse_sites(a.html):
        print(f"{sid} {name} ...", flush=True)
        sites[sid] = {"sun": series(lon, lat, "10", start, n, chunk), "ea": series(lon, lat, "399", start, n, chunk)}
    json.dump({"source": "NASA JPL Horizons API (observer table, Moon-fixed site, airless apparent az/el)",
               "start_ms": int(start.timestamp() * 1000), "step_h": 1, "n": n, "sites": sites},
              open(a.out, "w"), separators=(",", ":"))
    print("wrote", a.out, n, "hourly steps for", len(sites), "sites")
if __name__ == "__main__": main()
