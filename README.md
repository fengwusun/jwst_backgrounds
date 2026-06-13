jwst_backgrounds is a a simple program to predict the levels of background emission
in JWST observations, for use in proposal planning.

It accesses a precompiled background cache prepared by Space Telescope Science Institute. The background cache is hosted by the 
Mikulski Archive for Space Telescopes (MAST), so you need internet access to run the tool with the remote cache. It is possible to
download the full background cache to your local machine. Instructions for downloading the background cache can be found at http://archive.stsci.edu/archive_news/2017/08-Aug/index.html#article1

For a given target (RA, DEC), and wavelength, jwst_backgrounds does the following:
- Plot the spectrum of the background for that target on a given calendar day.
- Plot the total background for that target versus calendar day.
- Compute the number of days per year that the target is observable at low background,
  for a given wavelength and a selectable threshold.
- Save the retrieved background data to file.
  
This code was written by Jane Rigby (GSFC, Jane.Rigby@nasa.gov) and Klaus Pontoppidan (STScI, pontoppi@stsci.edu)
The background cache was prepared by Wayne Kinzel at STScI, and is the same as used by the JWST Exposure Time Calculator.

This software is provided as-is, with no warranty.


## Fast fork (this repository)

This is a performance fork of the STScI `jwst_backgrounds` tool. It returns the
**same numbers** as the official package (verified bit-for-bit; see
`tests/`), but is dramatically faster and works offline for the common case:

* **Pre-baked deep fields.** The background cache for the JWST extragalactic
  deep / calibration fields (GOODS-N/S, COSMOS, UDS, EGS, Abell 2744, Abell 370,
  NEP-TDF, LMC) is shipped *inside the package* (`jwst_backgrounds/field_cache/`,
  byte-identical to STScI). Those positions need **no internet** at all.
* **Real caching.** Every cache file downloaded at runtime is cached to disk
  (`$JBT_CACHE_DIR`, default `~/.cache/jwst_backgrounds`) and to memory, keyed by
  healpix pixel. The original re-downloaded the file *and* the `VERSION` file on
  every call.
* **Vectorized internals.** The binary cache is parsed with `numpy`
  (no per-day Python loop) and the bathtub curve is interpolated without
  rebuilding five `scipy.interp1d` objects per call.
* **Batch API** for catalogs (`get_backgrounds`) that groups sources by healpix
  pixel — a whole survey field costs a single file read.
* **Composable plots:** `plot_background`/`plot_bathtub` accept `ax=` and return
  `(fig, ax)` instead of only calling `plt.show()`.

Typical speedups: a single repeated query is ~10,000× faster (warm cache); a
catalog of a few thousand sources in one field goes from minutes to ~1 second.

Output stays consistent with the official package because the default cache is
the same one STScI ships (`sl_cache_2.0`). The newer `sl_cache_3.0` (regenerated
2025; identical at >2 µm, differs only at 0.5–1.2 µm) is available via
`cache_url=jbt.CACHE_URLS["3.0"]`.

### Quick start (fast fork)

```python
from jwst_backgrounds import jbt

# A named deep field — served from the shipped bundle, no network:
bkg = jbt.background.from_field("GOODS-N", 3.5)
print(bkg.bathtub["good_days"], "good days at 3.5 um")

# Any position (downloads + caches on first use):
bkg = jbt.background(189.2, 62.2, 3.5)

# A whole catalog at once -> dict of numpy arrays:
import numpy as np
ra  = np.array([189.20, 189.22, 189.25])
dec = np.array([ 62.22,  62.24,  62.26])
tab = jbt.get_backgrounds(ra, dec, wavelength=3.5)   # bkg_min/mean/max, good_days, ...
```

### Re-baking the field bundle

To refresh the bundle, add fields, or bake from `sl_cache_3.0`:

```
python -m jwst_backgrounds.bake                              # re-bake the default set
python -m jwst_backgrounds.bake --fields GOODS-N COSMOS      # a subset
python -m jwst_backgrounds.bake --cache 3.0 --radius 0.3     # different cache / radius
```

Field names and centres live in `jwst_backgrounds/fields.py`; many aliases
(`CEERS`→EGS, `A2744`→Abell2744, `CDFS`→GOODS-S, …) are accepted. List them with
`jwst_backgrounds --list-fields`.


INSTALLATION

Using pip:
----------
```
pip install jwst_backgrounds
```

Note: healpy (version >= 1.10) is a required dependency, so if you don't have it pip will install it automatically. 

Note: to upgrade the JBT with pip use `pip install jwst_background --upgrade`

Using Conda
-----------
First clone the repository

```
git clone git@github.com:fengwusun/jwst_backgrounds.git
cd jwst_backgrounds
conda create --name <env> --file requirements.txt
```

where `<env>` is the name of the environment you wish to create and requirements is the `requirements.txt` in the package directory.
To activate your JBT enter the following command:

```
source activate <env>
```

Manually
----------
Clone the repository from github and install using `easy_install`.

```
git clone git@github.com:fengwusun/jwst_backgrounds.git
cd jwst_backgrounds
easy_install .
```

   
RUNNING THE CODE:
```
python			# Start python.
from jwst_backgrounds import jbt 	# Import the background module
```

Below is an example that plots a background curve for a given RA, DEC, wavelength, threshold
```
jbt.get_background(261.6833333, -73.33222222, 2.15, thresh=1.1, \
                        plot_background=True, plot_bathtub=True, write_bathtub=True) 
```

Contributing
--------------
`jwst_backgrounds` follows the STScI ["forking workflow"](https://github.com/spacetelescope/style-guides/blob/master/guides/git-workflow.md#forking-workflow).


TROUBLESHOOTING:
-----------
If matplotlib does not display the images, then try editing your ~/.matplotlib/matplotlibrc file,
and choosing a different backend:  
```
backend: MacOSX
backend: TkAgg
backend: GTKCairo
```

Citation
--------
This code was written by Jane Rigby (GSFC, Jane.Rigby@nasa.gov) and Klaus Pontoppidan (STScI, pontoppi@stsci.edu) The background cache was prepared by Wayne Kinzel at STScI, and is the same as used by the JWST Exposure Time Calculator.

