# jwst_backgrounds — fast fork

A drop-in, much faster fork of the STScI [`jwst_backgrounds`](https://github.com/spacetelescope/jwst_backgrounds)
tool (JBT). It predicts the JWST background (zodiacal light, in-field galactic
ISM/CIB, scattered stray light, and thermal self-emission) at any sky position
and wavelength, versus calendar day.

It returns the **same numbers as the official package** (verified bit-for-bit —
see [Consistency](#consistency)), but is built for repeated/offline use and for
catalogs.

## What this fork adds

* **Pre-baked deep fields.** The background cache for 24 JWST extragalactic
  deep, lensing-cluster, wide-survey, and calibration fields is shipped *inside
  the package* (`jwst_backgrounds/field_cache/`, 111 healpix pixels, ~34 MB,
  **byte-identical to STScI**). Those positions need **no internet** — ever.
* **Snap-to-nearest.** A position that isn't itself pre-baked but lies within
  `snap_deg` of a baked field reuses that field's cache instead of downloading
  (**on by default at 1°**, ~1% approximation; see [Snapping](#snapping)).
* **Real caching.** Files downloaded at runtime are cached to disk
  (`$JBT_CACHE_DIR`, default `~/.cache/jwst_backgrounds`) and memory, keyed by
  healpix pixel. The original re-downloaded the data file *and* the `VERSION`
  file on **every** call.
* **Vectorized internals.** The binary cache is parsed with `numpy.frombuffer`
  (no per-day Python loop), and the bathtub curve is interpolated without
  rebuilding five `scipy.interp1d` objects per call.
* **Batch API** (`get_backgrounds`) that groups a catalog by healpix pixel — a
  whole survey field costs a single file read.
* **Field registry + re-baker**, `background.from_field("GOODS-N", 3.5)`,
  composable plots that return `(fig, ax)`, and CLI flags `--field`,
  `--list-fields`, `--cache`, `--snap`.
* **Optional `sl_cache_3.0`** (the regenerated 2025 cache) via
  `cache_url=jbt.CACHE_URLS["3.0"]`.

## Benchmarks

Single position = GOODS-N (a pre-baked field); catalog = 3,000 sources spread
across that field. Measured on a laptop over a normal connection (your network
latency dominates the official numbers).

| operation | official `jwst_backgrounds` | this fork | speed-up |
|---|---|---|---|
| one query, repeated (warm) | ~250 ms each (re-downloads every call) | **0.027 ms** | **~9,000×** |
| one query, first call in a process | ~250 ms | **3 ms** (reads local bundle, no network) | ~80× |
| 3,000-source catalog, one field | ~12 min (3,000 × 250 ms) | **0.05 s** | **~16,000×** |
| parse one cache file (CPU only) | ~1.0 ms | **0.02 ms** | ~45× |
| bathtub interpolation (CPU only) | ~0.20 ms | **0.007 ms** | ~30× |

The official tool makes two HTTPS requests on every single call and never
caches; ~95% of its wall-clock is network. For a survey field — which falls in
one or a few healpix pixels — it re-downloads the *same* file once per source.

## Pre-baked fields

```
jwst_backgrounds --list-fields
```

Wide/deep: **GOODS-N, GOODS-S, COSMOS, UDS, EGS, NEP-TDF**.
Lensing clusters: **Abell 2744, Abell 370, MACS 0416, MACS 0717, MACS 1149,
AS1063, MACS 1423, El Gordo, SMACS 0723**.
Wide-area / Euclid / Roman: **SSA22, SXDF, XMM-LSS, EDF-N, EDF-S, EDF-F,
Roman-TDF** — these multi-degree fields bake only a small exact core, so use
`snap_deg` (≈1°) to cover the rest of the footprint.
Calibration: **LMC, SMC**.

Field names are matched forgivingly — case-, space- and punctuation-insensitive,
with a large alias table, glob patterns, and fuzzy (typo) fallback:

```python
jbt.background.from_field("HUDF", 3.5)            # -> GOODS-S (so are CDF-S,
jbt.background.from_field("JADES deep field", 3.5) #    NGDEEP, "GOODS South", ...)
jbt.background.from_field("CEERS", 3.5)           # -> EGS
jbt.background.from_field("GOODS*N", 3.5)         # glob -> GOODS-N
fields.resolve("cosmso")                          # typo -> COSMOS
fields.search("MACS*")                            # -> ['MACS0416','MACS0717',...]
```

Ambiguous queries (`"MACS*"`, `"GOODS"`, `"Abell"`) raise with the candidate
list. Centres, radii, and the alias table live in
[`jwst_backgrounds/fields.py`](jwst_backgrounds/fields.py).

## Snapping

The STScI model is tabulated per healpix pixel at NSIDE=128 (~0.46° cells) with
**no** sub-pixel interpolation — the official tool already maps any query to its
containing pixel. The `snap_deg` option extends that idea: if a query's own
pixel isn't in the shipped bundle, reuse the **nearest baked pixel** within
`snap_deg` degrees (and flag it on the result).

How approximate is that? Worst-case difference between a pixel and one a given
distance away (max over N/S/E/W, on common calendar days):

| separation | @3.5 µm (NIR) | @21 µm (MIRI) | 1–25 µm (max) |
|---|---|---|---|
| 0.5° | ≲1.0% | ≲0.4% | ≲1.4% |
| **1.0°** | **≲1.5%** | **≲0.6%** | **≲2.5%** |
| 2.0° | ≲3.2% | ≲1.0% | ≲8.6% |

The MIRI range barely moves because it is thermal-dominated and thermal is
field-independent. **Snapping is on by default at 1°** (`snap_deg=1.0`): a query
within 1° of a baked field reuses it (and sets `.snapped`); queries farther than
1° from any field still download their own exact pixel. For results
byte-identical to STScI at every coordinate, pass `snap_deg=0`:

```python
bkg = jbt.background(189.9, 62.5, 7.7)              # within 1 deg of GOODS-N -> snaps
print(bkg.snapped)        # {'field': 'GOODS-N', 'pixel': ..., 'sep_deg': 0.7} or None
bkg = jbt.background(189.9, 62.5, 7.7, snap_deg=0)  # force the exact pixel (downloads)
```

## Quick start

```python
from jwst_backgrounds import jbt
import numpy as np

# A named deep field — served from the shipped bundle, no network:
bkg = jbt.background.from_field("GOODS-N", 3.5)
print(bkg.bathtub["good_days"], "good days at 3.5 um")

# Any position (downloads + caches on first use):
bkg = jbt.background(189.2, 62.2, 3.5)

# A whole catalog at once -> dict of numpy arrays
# (bkg_min/mean/max, good_days, ndays_observable, healpix, snap_sep_deg, ...):
ra  = np.array([189.20, 189.22, 189.25])
dec = np.array([ 62.22,  62.24,  62.26])
out = jbt.get_backgrounds(ra, dec, wavelength=3.5)

# Compose plots (returns fig, ax; pass ax= to place them):
import matplotlib.pyplot as plt
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5))
bkg.plot_background(ax=a, show=False)
bkg.plot_bathtub(ax=b, showsubbkgs=True, show=False)
```

Command line:

```
jwst_backgrounds --field GOODS-N 3.5          # named field, from the bundle
jwst_backgrounds 189.2 62.2 3.5 --cache 3.0   # explicit coords, newer cache
jwst_backgrounds 189.9 62.5 7.7 --snap 1.0    # snap to nearest baked field
```

## Re-baking the field bundle

To refresh the bundle, add fields, or bake from `sl_cache_3.0`:

```
python -m jwst_backgrounds.bake                          # re-bake the default set
python -m jwst_backgrounds.bake --fields GOODS-N COSMOS  # a subset
python -m jwst_backgrounds.bake --cache 3.0 --radius 0.3 # different cache / radius
```

Add your own fields by editing `DEEP_FIELDS` in
[`jwst_backgrounds/fields.py`](jwst_backgrounds/fields.py), then re-bake.

## Consistency

The default cache is the one the official package ships (`sl_cache_2.0`), and
the shipped `.bin` files are byte-identical to STScI. The test suite proves it:

```
pytest            # offline: vectorized parse == struct reference, bathtub == interp1d
pytest -m network # also checks every bundle file == a fresh STScI download
```

A direct comparison against the installed official package over all baked fields
and a wide range of wavelengths gives a maximum difference of **0.0** (bit-for-bit).

The default `snap_deg=1.0` approximates coordinates within 1° of a baked field by
~1% (less in MIRI); pass `snap_deg=0` for output identical to the official
package at every coordinate.

## Installation

```
pip install git+https://github.com/fengwusun/jwst_backgrounds
```

or from a clone:

```
git clone git@github.com:fengwusun/jwst_backgrounds.git
cd jwst_backgrounds
pip install -e .
```

Dependencies (installed automatically): `healpy>=1.10`, `numpy>=1.17`,
`scipy>=1.1`, `matplotlib>=3.1.1`.

## Troubleshooting

If matplotlib does not display the figures, choose a different backend in your
`~/.matplotlib/matplotlibrc` (`backend: MacOSX`, `TkAgg`, …), or use the
returned `(fig, ax)` to save them directly.

## Credits

The original `jwst_backgrounds` was written by **Jane Rigby** (GSFC) and
**Klaus Pontoppidan** (STScI); the background cache was prepared by Wayne Kinzel
at STScI and is the same one used by the JWST Exposure Time Calculator. This is a
performance fork maintained by **Fengwu Sun**; all background physics and data
are unchanged. Software is provided as-is, with no warranty — confirm
observability with APT and the ETC.
```

