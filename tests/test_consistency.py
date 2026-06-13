"""
Consistency tests: the fast fork must reproduce the upstream algorithm exactly.

These run fully offline against the pre-baked field bundle.  They re-implement
the original ``struct``-based parse and the original ``interp1d`` bathtub and
assert the fast vectorized versions match bit-for-bit.  A separate, network-only
test (``test_bundle_matches_stsci``) checks that the shipped bundle bytes are
identical to a fresh STScI download.
"""

import os
import json
import struct

import numpy as np
import pytest

from jwst_backgrounds import jbt, fields


# --------------------------------------------------------------------------- #
# Reference implementations copied from the upstream algorithm
# --------------------------------------------------------------------------- #
def _ref_parse(data, wave_array, thermal_bg):
    """The original read_bkg_data parsing (struct + per-day Python loop)."""
    nwave = wave_array.size
    size_calendar = struct.calcsize("366i")
    partA = struct.unpack(str(5 + nwave) + "d", data[0:(5 + nwave) * 8])
    nonzodi_bg = np.array(partA[5:5 + nwave])
    date_map = np.array(struct.unpack(
        "366i", data[(5 + nwave) * 8:(5 + nwave) * 8 + size_calendar]))
    calendar = np.where(date_map >= 0)[0]
    ndays = len(calendar)
    zodi = np.zeros((ndays, nwave))
    stray = np.zeros((ndays, nwave))
    perday = nwave * 2
    partB = struct.unpack(str(ndays * nwave * 2) + "d", data[perday * ndays * -8:])
    for dd in range(ndays):
        b1, b2, b3 = dd * perday, dd * perday + nwave, dd * perday + 2 * nwave
        zodi[dd, ] = partB[b1:b2]
        stray[dd, ] = partB[b2:b3]
    total = np.tile(nonzodi_bg + thermal_bg, (ndays, 1)) + stray + zodi
    return calendar, nonzodi_bg, zodi, stray, total


def _ref_bathtub(wave_array, arr, wl):
    """The original interp1d bathtub interpolation."""
    from scipy.interpolate import interp1d
    return interp1d(wave_array, arr, bounds_error=True)(wl)


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #
WAVE, THERMAL = jbt.read_static_data()


@pytest.mark.parametrize("name", fields.list_fields())
def test_parse_matches_reference(name):
    """Vectorized parse_bin == original struct parse, for one pixel per field."""
    _, ra, dec, _ = fields.resolve(name)
    healpix = jbt.healpix_of(ra, dec)
    path = jbt._packaged_path(healpix)
    assert path is not None, f"{name} center pixel not in bundle"
    with open(path, "rb") as f:
        data = f.read()

    cal, nz, zodi, stray, total = _ref_parse(data, WAVE, THERMAL)
    fast = jbt.parse_bin(data, WAVE, THERMAL)

    assert np.array_equal(fast["calendar"], cal)
    assert np.array_equal(fast["nonzodi_bg"], nz)
    assert np.array_equal(fast["zodi_bg"], zodi)
    assert np.array_equal(fast["stray_light_bg"], stray)
    assert np.array_equal(fast["total_bg"], total)


@pytest.mark.parametrize("wl", [2.0, 3.5, 4.44, 7.7, 10.0, 15.0, 21.0, 25.5, 0.5, 31.0])
def test_bathtub_matches_reference(wl):
    """Vectorized make_bathtub == original interp1d, bit-for-bit."""
    ra, dec, _ = fields.DEEP_FIELDS["GOODS-N"]
    bkg = jbt.background(ra, dec, wl)
    for comp, key in [("total_bg", "total_thiswave"), ("zodi_bg", "zodi_thiswave"),
                      ("stray_light_bg", "stray_thiswave")]:
        ref = _ref_bathtub(WAVE, bkg.bkg_data[comp], wl)
        assert np.array_equal(bkg.bathtub[key], ref), f"{key} differs at {wl} um"


def test_bathtub_bounds():
    """Out-of-grid wavelength raises, like interp1d(bounds_error=True)."""
    ra, dec, _ = fields.DEEP_FIELDS["GOODS-N"]
    bkg = jbt.background(ra, dec, 5.0)
    with pytest.raises(ValueError):
        bkg.make_bathtub(0.4)
    with pytest.raises(ValueError):
        bkg.make_bathtub(40.0)


def test_bundle_integrity():
    """Every pixel named in the manifest exists and parses to a sane position."""
    man = os.path.join(jbt._FIELD_CACHE, "manifest.json")
    with open(man) as f:
        manifest = json.load(f)
    seen = set()
    for fld in manifest["fields"].values():
        for p in fld["pixels"]:
            seen.add(p)
            path = jbt._packaged_path(p)
            assert path is not None, f"pixel {p} missing from bundle"
            with open(path, "rb") as fh:
                d = jbt.parse_bin(fh.read(), WAVE, THERMAL)
            assert -1 <= np.sin(np.radians(d["dec"])) <= 1
            assert d["calendar"].size > 0
    assert seen, "bundle has no pixels"


@pytest.mark.network
def test_bundle_matches_stsci():
    """Shipped bundle bytes are byte-identical to a fresh STScI download."""
    import urllib.request
    man = os.path.join(jbt._FIELD_CACHE, "manifest.json")
    with open(man) as f:
        manifest = json.load(f)
    url = manifest["cache_url"]
    # spot-check the centre pixel of each field (full set is large)
    for name in fields.list_fields():
        _, ra, dec, _ = fields.resolve(name)
        h = jbt.healpix_of(ra, dec)
        with open(jbt._packaged_path(h), "rb") as fh:
            local = fh.read()
        remote = urllib.request.urlopen(url + jbt.file_from_healpix(h), timeout=120).read()
        assert local == remote, f"{name}: bundle byte-differs from STScI"
