"""
Predict JWST background levels for proposal planning -- fast fork.

This is a drop-in compatible, much faster fork of the STScI ``jwst_backgrounds``
tool.  It accesses the same precompiled STScI background cache and returns the
**same numbers** (see ``tests/test_consistency.py``), but:

* ships a **pre-baked cache** for the JWST extragalactic deep / calibration
  fields (GOODS-N/S, COSMOS, UDS, EGS, Abell 2744/370, NEP-TDF, LMC), so those
  positions need no network access at all;
* **caches** every downloaded cache file to disk and memory, keyed by healpix
  pixel, instead of re-downloading (and re-fetching ``VERSION``) on every call;
* parses the binary cache with vectorized ``numpy`` instead of a per-day Python
  loop, and interpolates the bathtub curve without rebuilding five
  ``interp1d`` objects per call;
* adds a **batch API** (:func:`get_backgrounds`) that groups a catalog by
  healpix pixel -- a whole survey field costs a single file read; and
* returns Matplotlib ``(fig, ax)`` from the plot helpers so figures compose.

The public API of the original (``background``, ``get_background``, and all of
``background``'s methods/attributes) is preserved.

Software is provided as-is, with no warranty.  Use the latest versions of APT
and ETC to confirm the observability of any JWST targets.
"""

import os
import json
import warnings
import urllib.request
import urllib.error

import healpy
import numpy as np
from scipy.interpolate import interp1d

from jwst_backgrounds.version import __version__

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
# Default cache == what the official STScI package uses, so default output
# matches it exactly.  sl_cache_3.0 (regenerated 2025) is also available; it is
# identical at >2 micron and differs only at 0.5-1.2 micron (the E. Han
# correction).  Select it with cache_url=CACHE_URLS["3.0"].
CACHE_URLS = {
    "2.0": "https://archive.stsci.edu/missions/jwst/simulations/straylight/sl_cache_2.0/",
    "3.0": "https://archive.stsci.edu/missions/jwst/simulations/straylight/sl_cache_3.0/",
}
DEFAULT_CACHE_URL = CACHE_URLS["2.0"]

NSIDE = 128                       # healpy tessellation of the background cache
WAVE_FILE = "std_spectrum_wavelengths.txt"
THERMAL_FILE = "thermal_curve_jwst_jrigby_v4.0.csv"
DEFAULT_SNAP_DEG = 0.0            # >0: snap an off-bundle query to the nearest baked
                                  # pixel within this many degrees (approximate)

_PKG_DIR = os.path.dirname(__file__)
_REFDATA = os.path.join(_PKG_DIR, "refdata")
_FIELD_CACHE = os.path.join(_PKG_DIR, "field_cache")   # pre-baked .bin bundle

# On-disk cache for files downloaded at runtime (override with $JBT_CACHE_DIR).
USER_CACHE_DIR = os.environ.get(
    "JBT_CACHE_DIR",
    os.path.join(os.path.expanduser("~"), ".cache", "jwst_backgrounds"),
)

# Process-lifetime memoization.
_STATIC_CACHE = {}   # thermal_file -> (wave_array, thermal_bg)
_BYTES_CACHE = {}    # (cache_url, healpix) -> raw bytes
_PARSE_CACHE = {}    # (cache_url, healpix, thermal_file) -> bkg_data dict
_VERSION_CACHE = {}  # cache_url -> version string
_BAKED_INDEX = None  # cached (pixels, unit vectors, field names) for snapping


# --------------------------------------------------------------------------- #
# Static reference data (loaded once per process)
# --------------------------------------------------------------------------- #
def read_static_data(thermal_file=THERMAL_FILE):
    """Return ``(wave_array, thermal_bg)``, memoized across all calls.

    ``thermal_bg`` is the constant thermal self-emission curve interpolated onto
    the standard wavelength grid (identical to the upstream computation).
    """
    if thermal_file in _STATIC_CACHE:
        return _STATIC_CACHE[thermal_file]
    wave_array = np.loadtxt(os.path.join(_REFDATA, WAVE_FILE))
    thermal = np.genfromtxt(os.path.join(_REFDATA, thermal_file), delimiter=",")
    # interp1d with these settings matches upstream: 0.0 outside the tabulated range
    thermal_bg = interp1d(thermal[:, 0], thermal[:, 1],
                          bounds_error=False, fill_value=0.0)(wave_array)
    _STATIC_CACHE[thermal_file] = (wave_array, thermal_bg)
    return wave_array, thermal_bg


# --------------------------------------------------------------------------- #
# Healpix bookkeeping (identical convention to upstream)
# --------------------------------------------------------------------------- #
def healpix_of(ra, dec):
    """Healpix RING pixel(s) for RA/DEC in decimal degrees (scalar or array)."""
    return healpy.pixelfunc.ang2pix(NSIDE, ra, dec, nest=False, lonlat=True)


def file_from_healpix(healpix):
    """Remote/relative path of the cache file for a healpix pixel."""
    s = str(int(healpix)).zfill(6)
    return s[0:4] + "/sl_pix_" + s + ".bin"


def _packaged_path(healpix):
    """Path of a pre-baked .bin inside the package (flat layout), or None."""
    p = os.path.join(_FIELD_CACHE, "sl_pix_%06d.bin" % int(healpix))
    return p if os.path.isfile(p) else None


def _bundle_cache_url():
    """cache_url the shipped field bundle was baked from (from its manifest)."""
    man = os.path.join(_FIELD_CACHE, "manifest.json")
    if os.path.isfile(man):
        try:
            with open(man) as f:
                return json.load(f).get("cache_url", DEFAULT_CACHE_URL)
        except (OSError, ValueError):
            pass
    return DEFAULT_CACHE_URL


def _baked_index():
    """Cached description of the baked pixels, for nearest-pixel snapping.

    Returns a dict with 'pixels' (int array), 'vecs' (N x 3 unit vectors of the
    pixel centres) and 'field' (list of field names), or None if no bundle.
    """
    global _BAKED_INDEX
    if _BAKED_INDEX is not None:
        return _BAKED_INDEX or None
    man = os.path.join(_FIELD_CACHE, "manifest.json")
    if not os.path.isfile(man):
        _BAKED_INDEX = {}
        return None
    with open(man) as f:
        manifest = json.load(f)
    pix_field = {}
    for name, info in manifest.get("fields", {}).items():
        for p in info["pixels"]:
            pix_field.setdefault(int(p), name)   # first field wins for shared pixels
    pixels = np.array(sorted(pix_field), dtype=int)
    if pixels.size == 0:
        _BAKED_INDEX = {}
        return None
    vecs = np.array(healpy.pix2vec(NSIDE, pixels, nest=False)).T   # (N, 3)
    _BAKED_INDEX = {"pixels": pixels, "vecs": vecs,
                    "field": [pix_field[int(p)] for p in pixels]}
    return _BAKED_INDEX


def nearest_baked(ra, dec):
    """Nearest baked pixel to a position: ``(healpix, separation_deg, field_name)``.

    Returns ``(None, inf, None)`` if no field bundle is present.
    """
    idx = _baked_index()
    if idx is None:
        return None, float("inf"), None
    v = np.asarray(healpy.ang2vec(ra, dec, lonlat=True))
    dots = np.clip(idx["vecs"] @ v, -1.0, 1.0)
    k = int(np.argmax(dots))
    sep = float(np.degrees(np.arccos(dots[k])))
    return int(idx["pixels"][k]), sep, idx["field"][k]


def _resolve_pixel(ra, dec, snap_deg):
    """Pick the healpix pixel to read for a query, applying snap if enabled.

    Snapping only kicks in when the query's own pixel is not in the shipped
    bundle and a baked pixel lies within ``snap_deg``; then the nearest baked
    pixel is used (approximate). Returns ``(healpix, snapped_info_or_None)``.
    """
    h = int(healpix_of(ra, dec))
    if snap_deg and snap_deg > 0 and _packaged_path(h) is None:
        p, sep, field = nearest_baked(ra, dec)
        if p is not None and sep <= snap_deg:
            return p, {"field": field, "pixel": p, "sep_deg": sep, "query_healpix": h}
    return h, None


# --------------------------------------------------------------------------- #
# Fetching + caching the binary cache files
# --------------------------------------------------------------------------- #
def get_cache_version(cache_url=DEFAULT_CACHE_URL, timeout=30):
    """Cache VERSION string (once per process).

    For the baked cache the shipped manifest records the authoritative version,
    so we use it directly (no network); otherwise we fetch it from STScI and
    fall back to the manifest when offline.
    """
    if cache_url in _VERSION_CACHE:
        return _VERSION_CACHE[cache_url]
    man = os.path.join(_FIELD_CACHE, "manifest.json")

    def _from_manifest():
        if cache_url == _bundle_cache_url() and os.path.isfile(man):
            try:
                with open(man) as f:
                    return json.load(f).get("version")
            except (OSError, ValueError):
                return None
        return None

    version = _from_manifest()
    if version is None:
        try:
            with urllib.request.urlopen(cache_url + "VERSION", timeout=timeout) as fh:
                version = fh.readlines()[0].decode("utf-8").rstrip("\n")
        except (urllib.error.URLError, OSError, IndexError):
            version = None
    _VERSION_CACHE[cache_url] = version
    return version


def fetch_bin(healpix, cache_url=DEFAULT_CACHE_URL, cache_dir=USER_CACHE_DIR,
              timeout=60, verbose=False):
    """Return the raw bytes of the cache file for ``healpix``.

    Resolution order: in-memory -> packaged bundle -> on-disk user cache ->
    download (written through to the on-disk cache).
    """
    healpix = int(healpix)
    key = (cache_url, healpix)
    if key in _BYTES_CACHE:
        return _BYTES_CACHE[key]

    # 1) pre-baked bundle shipped with the package (only for the baked cache_url)
    if cache_url == _bundle_cache_url():
        p = _packaged_path(healpix)
        if p is not None:
            with open(p, "rb") as f:
                data = f.read()
            if verbose:
                print("[jbt] bundle  sl_pix_%06d.bin" % healpix)
            _BYTES_CACHE[key] = data
            return data

    # 2) on-disk user cache, mirroring the remote directory layout
    rel = file_from_healpix(healpix)
    local = os.path.join(cache_dir, os.path.basename(cache_url.rstrip("/")), rel)
    if os.path.isfile(local):
        with open(local, "rb") as f:
            data = f.read()
        if verbose:
            print("[jbt] disk    %s" % rel)
        _BYTES_CACHE[key] = data
        return data

    # 3) download, then write through to disk (atomic) and memory
    if verbose:
        print("[jbt] download %s" % rel)
    with urllib.request.urlopen(cache_url + rel, timeout=timeout) as fh:
        data = fh.read()
    os.makedirs(os.path.dirname(local), exist_ok=True)
    tmp = local + ".part"
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, local)
    _BYTES_CACHE[key] = data
    return data


# --------------------------------------------------------------------------- #
# Vectorized binary parse
# --------------------------------------------------------------------------- #
def parse_bin(data, wave_array, thermal_bg):
    """Parse one cache .bin blob into the bkg_data dict (vectorized).

    Binary schema (native byte order), verified against the upstream
    ``generate_stray_light`` source:

        double RA, DEC, pos[3], nonzodi_bg[NWAVE]
        int32  date_map[366]
        for each observable day: double zodi_bg[NWAVE], stray_light_bg[NWAVE]
    """
    nwave = wave_array.size
    nA = 5 + nwave

    partA = np.frombuffer(data, dtype="f8", count=nA, offset=0)
    ra, dec = float(partA[0]), float(partA[1])
    pos = partA[2:5]
    nonzodi_bg = partA[5:5 + nwave]

    off = nA * 8
    date_map = np.frombuffer(data, dtype="i4", count=366, offset=off)
    calendar = np.where(date_map >= 0)[0]
    ndays = calendar.size

    partB = np.frombuffer(data, dtype="f8", offset=off + 366 * 4).reshape(ndays, 2, nwave)
    zodi_bg = partB[:, 0, :]
    stray_light_bg = partB[:, 1, :]

    total_bg = (nonzodi_bg + thermal_bg)[None, :] + stray_light_bg + zodi_bg

    return {
        "calendar": calendar, "ra": ra, "dec": dec, "pos": pos,
        "wave_array": wave_array, "nonzodi_bg": nonzodi_bg, "thermal_bg": thermal_bg,
        "zodi_bg": zodi_bg, "stray_light_bg": stray_light_bg, "total_bg": total_bg,
    }


def get_bkg_data_pixel(healpix, thermal_file=THERMAL_FILE, cache_url=DEFAULT_CACHE_URL,
                       cache_dir=USER_CACHE_DIR, verbose=False):
    """bkg_data dict for an explicit healpix pixel (cached)."""
    wave_array, thermal_bg = read_static_data(thermal_file)
    healpix = int(healpix)
    key = (cache_url, healpix, thermal_file)
    if key in _PARSE_CACHE:
        return _PARSE_CACHE[key]
    data = fetch_bin(healpix, cache_url=cache_url, cache_dir=cache_dir, verbose=verbose)
    parsed = parse_bin(data, wave_array, thermal_bg)
    _PARSE_CACHE[key] = parsed
    return parsed


def get_bkg_data(ra, dec, thermal_file=THERMAL_FILE, cache_url=DEFAULT_CACHE_URL,
                 cache_dir=USER_CACHE_DIR, verbose=False):
    """bkg_data dict for a position, using all caches (network at most once)."""
    return get_bkg_data_pixel(int(healpix_of(ra, dec)), thermal_file, cache_url,
                              cache_dir, verbose=verbose)


# --------------------------------------------------------------------------- #
# Vectorized single-wavelength interpolation ("bathtub")
# --------------------------------------------------------------------------- #
def _bracket(wave_array, wavelength):
    """Indices (j-1, j), x0, x1 for linear interp at a scalar wavelength.

    Mirrors ``interp1d(..., bounds_error=True)``: raises ValueError outside the
    grid, and reproduces its ``y0 + slope*(x-x0)`` arithmetic for equal results.
    """
    if wavelength < wave_array[0] or wavelength > wave_array[-1]:
        raise ValueError(
            "wavelength %g um is outside the cache range [%g, %g] um"
            % (wavelength, wave_array[0], wave_array[-1]))
    j = int(np.searchsorted(wave_array, wavelength))
    j = min(max(j, 1), wave_array.size - 1)
    return j, wave_array[j - 1], wave_array[j]


def make_bathtub_dict(bkg_data, wavelength, thresh=1.1):
    """Background vs calendar day at one wavelength (vectorized linear interp)."""
    wave_array = bkg_data["wave_array"]
    j, x0, x1 = _bracket(wave_array, wavelength)
    dx = wavelength - x0

    def lin(a):
        slope = (a[..., j] - a[..., j - 1]) / (x1 - x0)
        return a[..., j - 1] + slope * dx

    total = lin(bkg_data["total_bg"])
    themin = float(np.min(total))
    good_days = int(np.sum(total < themin * thresh))
    return {
        "wavelength": wavelength, "themin": themin, "good_days": good_days,
        "total_thiswave": total,
        "stray_thiswave": lin(bkg_data["stray_light_bg"]),
        "zodi_thiswave": lin(bkg_data["zodi_bg"]),
        "thermal_thiswave": lin(bkg_data["thermal_bg"]),
        "nonzodi_thiswave": lin(bkg_data["nonzodi_bg"]),
    }


# --------------------------------------------------------------------------- #
# background class -- API-compatible with the original
# --------------------------------------------------------------------------- #
class background():
    '''
    Main background class, initialized with all background data for a specific
    position (RA, DEC).  The wavelength at which the bathtub curve is calculated
    can be updated as needed via :meth:`make_bathtub`.

    Parameters
    ----------
    ra, dec : float
        Right ascension / declination in decimal degrees.
    wavelength : float
        Wavelength (micron) at which the bathtub curve is calculated.
    thresh : float
        Background threshold relative to the annual minimum (default 1.1).
    cache_url : str
        STScI cache to read; defaults to the same one the official package uses.
    cache_dir : str
        On-disk cache directory for files downloaded at runtime.

    Attributes
    ----------
    bkg_data : dict
        All background data for the input (RA, DEC).
    bathtub : dict
        Background vs calendar day, interpolated at ``wavelength``.
    '''

    def __init__(self, ra, dec, wavelength, thresh=1.1, thermal_file=THERMAL_FILE,
                 cache_url=DEFAULT_CACHE_URL, cache_dir=USER_CACHE_DIR,
                 snap_deg=DEFAULT_SNAP_DEG, verbose=False):
        # global attributes (kept for backwards compatibility)
        self.cache_url = cache_url
        self.local_path = _REFDATA
        self.wave_file = WAVE_FILE
        self.thermal_file = thermal_file
        self.cache_dir = cache_dir
        self.nside = NSIDE
        self.wave_array, self.thermal_bg = read_static_data(thermal_file)
        self.sl_nwave = self.wave_array.size

        # input parameters
        self.ra = ra
        self.dec = dec
        self.wavelength = wavelength
        self.thresh = thresh

        # load variable content (cached); optionally snap to nearest baked pixel
        self.snap_deg = snap_deg
        self.healpix = int(healpix_of(ra, dec))
        used_pixel, snapped = _resolve_pixel(ra, dec, snap_deg)
        self.used_healpix = used_pixel
        self.snapped = snapped
        if snapped is not None:
            warnings.warn(
                "Snapped (RA,DEC)=(%.5f,%.5f) to baked pixel %d near field '%s' "
                "(%.3f deg away); background is approximate (~1 percent in the NIR)."
                % (ra, dec, snapped["pixel"], snapped["field"], snapped["sep_deg"]),
                stacklevel=2)
        self.cache_file = file_from_healpix(used_pixel)
        self.cache_version = get_cache_version(cache_url)
        self.bkg_data = get_bkg_data_pixel(used_pixel, thermal_file, cache_url,
                                           cache_dir, verbose=verbose)

        # interpolate bathtub curve and package it
        self.make_bathtub(wavelength)

    @classmethod
    def from_field(cls, name, wavelength, **kwargs):
        """Construct from a named field (see :mod:`jwst_backgrounds.fields`)."""
        from jwst_backgrounds import fields
        _, ra, dec, _ = fields.resolve(name)
        return cls(ra, dec, wavelength, **kwargs)

    # --- kept for API compatibility -------------------------------------- #
    def myfile_from_healpix(self, ra, dec):
        return file_from_healpix(healpix_of(ra, dec))

    def read_static_data(self):
        return read_static_data(self.thermal_file)

    def read_bkg_data(self, cache_file=None, verbose=False):
        """Read and parse the cache file for this position (cached)."""
        return get_bkg_data(self.ra, self.dec, self.thermal_file, self.cache_url,
                            self.cache_dir, verbose=verbose)

    def interpolate_spec(self, wave, specin, new_wave, fill=np.nan):
        f = interp1d(wave, specin, bounds_error=False, fill_value=fill)
        return f(new_wave)

    def make_bathtub(self, wavelength):
        """Interpolate a bathtub curve at ``wavelength`` and count good days."""
        self.wavelength = wavelength
        self.bathtub = make_bathtub_dict(self.bkg_data, wavelength, self.thresh)
        return self.bathtub

    # --- plotting (now returns (fig, ax) and accepts ax/show) ------------ #
    def plot_background(self, fontsize=16, xrange=(0.6, 30), yrange=(1e-4, 1e4),
                        thisday=None, ax=None, show=True):
        """Plot the full background spectrum for one calendar day.  Returns (fig, ax)."""
        import matplotlib.pyplot as plt
        bd = self.bkg_data
        calendar = bd["calendar"]
        if thisday is None:
            thisday = int(calendar[calendar.size // 2])
        if thisday not in calendar:
            print("The input calendar day {} is not available".format(thisday))
            return None, None
        di = int(np.where(calendar == thisday)[0][0])

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        else:
            fig = ax.figure
        ax.plot(bd["wave_array"], bd["nonzodi_bg"], label="ISM")
        ax.plot(bd["wave_array"], bd["zodi_bg"][di, :], label="Zodi")
        ax.plot(bd["wave_array"], bd["stray_light_bg"][di, :], label="Stray light")
        ax.plot(bd["wave_array"], bd["thermal_bg"], label="Thermal")
        ax.plot(bd["wave_array"], bd["total_bg"][di, :], label="Total", color="black", lw=3)
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_yscale("log")
        ax.set_xlabel("wavelength (micron)", fontsize=fontsize)
        ax.set_ylabel("Equivalent in-field radiance (MJy/sr)", fontsize=fontsize)
        ax.set_title("Background for calendar day " + str(thisday))
        ax.legend()
        if show:
            plt.show()
        return fig, ax

    def plot_bathtub(self, showthresh=True, showplot=False, showsubbkgs=False,
                     showannotate=True, title=False, label=False, ax=None, show=True):
        """Plot background vs calendar day at ``self.wavelength``.  Returns (fig, ax)."""
        import matplotlib.pyplot as plt
        bathtub = self.bathtub
        if not label:
            label = "Total " + str(bathtub["wavelength"]) + " micron"
        calendar = self.bkg_data["calendar"]

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        else:
            fig = ax.figure
        ax.scatter(calendar, bathtub["total_thiswave"], s=20, label=label)
        ax.set_xlabel("Day of the year", fontsize=12)
        ax.set_xlim(0, 366)

        if showannotate:
            ax.set_title("{} good days out of {} days observable, for threshold {}".format(
                bathtub["good_days"], calendar.size, self.thresh))
            ax.set_ylabel("bkg at " + str(bathtub["wavelength"]) + " um (MJy/sr)", fontsize=12)
        else:
            ax.set_ylabel("bkg (MJy/SR)", fontsize=20)

        if showsubbkgs:
            ax.scatter(calendar, bathtub["zodi_thiswave"], s=20, label="Zodiacal")
            ax.scatter(calendar, bathtub["stray_thiswave"], s=20, label="Stray light")
            ax.scatter(calendar, bathtub["nonzodi_thiswave"] * np.ones_like(calendar),
                       s=20, label="ISM+CIB")
            ax.scatter(calendar, bathtub["thermal_thiswave"] * np.ones_like(calendar),
                       s=20, label="Thermal")
            ax.grid()
            ax.locator_params(axis="x", nbins=10)
            ax.locator_params(axis="y", nbins=10)
        if showthresh:
            ax.hlines((bathtub["themin"], bathtub["themin"] * self.thresh), 0, 365, color="black")
        if title:
            ax.set_title(title)
        ax.legend(fontsize=10, frameon=False, labelspacing=0)
        if show:
            plt.show()
        return fig, ax

    # --- text output ----------------------------------------------------- #
    def write_bathtub(self, bathtub_file="background_versus_day.txt"):
        with open(bathtub_file, "w") as f:
            f.write("# Output of JWST_backgrounds version " + str(__version__) + "\n")
            f.write("# background cache version " + str(self.cache_version) + "\n\n")
            f.write("# for RA=" + str(self.ra) + ", DEC=" + str(self.dec)
                    + " at wavelength=" + str(self.wavelength) + " micron \n")
            f.write("# Columns: \n# - Calendar day (Jan1=0) \n# - Total background (MJy/sr)\n")
            for i, calendar_day in enumerate(self.bkg_data["calendar"]):
                f.write("{0}    {1:5.4f}".format(calendar_day, self.bathtub["total_thiswave"][i]) + "\n")

    def write_background(self, background_file="background.txt", thisday=None):
        bd = self.bkg_data
        calendar = bd["calendar"]
        if thisday is None:
            thisday = int(calendar[calendar.size // 2])
        if thisday not in calendar:
            print("The input calendar day {} is not available".format(thisday))
            return
        di = int(np.where(calendar == thisday)[0][0])
        with open(background_file, "w") as f:
            f.write("# Output of JWST_backgrounds version " + str(__version__) + "\n")
            f.write("# background cache version " + str(self.cache_version) + "\n\n")
            f.write("# for RA=" + str(self.ra) + ", DEC=" + str(self.dec)
                    + " On calendar day " + str(thisday) + "\n")
            f.write("# Columns: \n# - Wavelength [micron] \n# - Total background (MJy/sr)\n"
                    "# - In-field zodiacal light (MJy/sr)\n# - In-field galactic light (MJy/sr)\n"
                    "# - Stray light (MJy/sr)\n# - Thermal self-emission (MJy/sr)\n")
            for i, wavelength in enumerate(bd["wave_array"]):
                f.write("{0:f}    {1:5.4f}    {2:5.4f}    {3:5.4f}    {4:5.4f}    {5:5.4f}".format(
                    wavelength, bd["total_bg"][di][i], bd["zodi_bg"][di][i], bd["nonzodi_bg"][i],
                    bd["stray_light_bg"][di][i], bd["thermal_bg"][i]) + "\n")


# --------------------------------------------------------------------------- #
# Batch API -- the big win for catalogs
# --------------------------------------------------------------------------- #
def get_backgrounds(ra, dec, wavelength, thresh=1.1, thermal_file=THERMAL_FILE,
                    cache_url=DEFAULT_CACHE_URL, cache_dir=USER_CACHE_DIR,
                    snap_deg=DEFAULT_SNAP_DEG, full=False, verbose=False):
    """Backgrounds for a whole catalog, grouped by healpix pixel.

    All sources sharing a pixel (e.g. an entire survey field) trigger a single
    file read; only the per-wavelength interpolation runs per source.

    Parameters
    ----------
    ra, dec : array-like
        Coordinates in decimal degrees.
    wavelength : float or array-like
        Wavelength(s) in micron.  Scalar applies to every source.
    full : bool
        If True, also return a list of per-source bathtub dicts.

    Returns
    -------
    dict of numpy arrays
        Keys: ra, dec, wavelength, healpix, ndays_observable, good_days,
        bkg_min, bkg_mean, bkg_max  (total background at ``wavelength``, MJy/sr).
        If ``full``, an extra key ``bathtubs`` holds the per-source dicts.
    """
    ra = np.atleast_1d(np.asarray(ra, dtype=float))
    dec = np.atleast_1d(np.asarray(dec, dtype=float))
    if ra.shape != dec.shape:
        raise ValueError("ra and dec must have the same shape")
    wl = np.atleast_1d(np.asarray(wavelength, dtype=float))
    if wl.size == 1:
        wl = np.full(ra.shape, wl[0])
    elif wl.shape != ra.shape:
        raise ValueError("wavelength must be scalar or match ra/dec shape")

    pix = healpix_of(ra, dec).astype(int)
    used = pix.copy()
    snap_sep = np.zeros(ra.size)
    snap_field = np.array([""] * ra.size, dtype=object)
    if snap_deg and snap_deg > 0:
        for i in range(ra.size):
            up, sn = _resolve_pixel(ra[i], dec[i], snap_deg)
            used[i] = up
            if sn is not None:
                snap_sep[i] = sn["sep_deg"]
                snap_field[i] = sn["field"]

    baths = [None] * ra.size
    wave_array, thermal_bg = read_static_data(thermal_file)
    uniq = np.unique(used)
    if verbose:
        nsnap = int(np.count_nonzero(snap_sep))
        print("[jbt] %d sources -> %d unique pixel(s)%s"
              % (ra.size, uniq.size, (", %d snapped" % nsnap) if nsnap else ""))
    for p in uniq:
        idx = np.where(used == p)[0]
        data = fetch_bin(p, cache_url=cache_url, cache_dir=cache_dir, verbose=verbose)
        parsed = parse_bin(data, wave_array, thermal_bg)
        for i in idx:
            baths[i] = make_bathtub_dict(parsed, float(wl[i]), thresh)

    tot = [b["total_thiswave"] for b in baths]
    out = {
        "ra": ra, "dec": dec, "wavelength": wl, "healpix": pix, "used_healpix": used,
        "snap_sep_deg": snap_sep, "snap_field": snap_field,
        "ndays_observable": np.array([t.size for t in tot]),
        "good_days": np.array([b["good_days"] for b in baths]),
        "bkg_min": np.array([b["themin"] for b in baths]),
        "bkg_mean": np.array([float(np.mean(t)) for t in tot]),
        "bkg_max": np.array([float(np.max(t)) for t in tot]),
    }
    if full:
        out["bathtubs"] = baths
    return out


def clear_cache(disk=False):
    """Clear in-memory caches; if ``disk`` also remove the runtime on-disk cache.

    The pre-baked package bundle is never removed.
    """
    global _BAKED_INDEX
    _BYTES_CACHE.clear()
    _PARSE_CACHE.clear()
    _VERSION_CACHE.clear()
    _BAKED_INDEX = None
    if disk and os.path.isdir(USER_CACHE_DIR):
        import shutil
        shutil.rmtree(USER_CACHE_DIR)


# --------------------------------------------------------------------------- #
# Convenience wrapper (unchanged behaviour from upstream)
# --------------------------------------------------------------------------- #
def get_background(ra, dec, wavelength, thresh=1.1, plot_background=True, plot_bathtub=True,
                   thisday=None, showsubbkgs=False, write_background=True, write_bathtub=True,
                   background_file="background.txt", bathtub_file="background_versus_day.txt",
                   cache_url=DEFAULT_CACHE_URL, snap_deg=DEFAULT_SNAP_DEG):
    """Get the background data and create plots/outputs with one command."""
    bkg = background(ra, dec, wavelength, thresh=thresh, cache_url=cache_url, snap_deg=snap_deg)
    calendar = bkg.bkg_data["calendar"]
    print("These coordinates are observable by JWST", len(calendar), "days per year.")
    print("For", bkg.bathtub["good_days"], "of those days, the background is <", thresh,
          "times the minimum, at wavelength", wavelength, "micron")

    if thisday not in calendar:
        ndays = calendar.size
        if ndays > 0:
            thisday_input = thisday
            thisday = calendar[int(ndays / 2)]
            print("Warning: The input calendar day {}".format(thisday_input)
                  + " is not available, assuming the middle day: {} instead".format(thisday))
        else:
            print("No valid days")
            return

    if plot_background:
        bkg.plot_background(thisday=thisday)
    if write_background:
        bkg.write_background(thisday=thisday, background_file=background_file)
    if plot_bathtub:
        bkg.plot_bathtub(showsubbkgs=showsubbkgs)
    if write_bathtub:
        bkg.write_bathtub(bathtub_file=bathtub_file)
    return bkg
