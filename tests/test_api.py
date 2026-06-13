"""Tests for the fork's new conveniences: field lookup, batch API, caching."""

import numpy as np
import pytest

from jwst_backgrounds import jbt, fields


def test_field_registry_and_aliases():
    canon, ra, dec, rad = fields.resolve("goods n")
    assert canon == "GOODS-N"
    assert fields.resolve("CEERS")[0] == "EGS"
    assert fields.resolve("a2744")[0] == "Abell2744"
    with pytest.raises(KeyError):
        fields.resolve("not-a-field")


@pytest.mark.parametrize("name", fields.list_fields())
def test_from_field_offline(name):
    """Every baked field is served from the bundle (works without network)."""
    bkg = jbt.background.from_field(name, 5.0)
    assert bkg.bkg_data["calendar"].size > 0
    assert np.isfinite(bkg.bathtub["themin"])
    assert 0 <= bkg.bathtub["good_days"] <= bkg.bkg_data["calendar"].size


def test_batch_matches_single():
    """get_backgrounds agrees with per-source background() calls."""
    ras, decs = [], []
    for name in ["GOODS-N", "GOODS-S", "COSMOS"]:
        _, ra, dec, _ = fields.resolve(name)
        ras.append(ra)
        decs.append(dec)
    wls = [3.5, 7.7, 21.0]
    tab = jbt.get_backgrounds(ras, decs, wls)
    for i in range(len(ras)):
        b = jbt.background(ras[i], decs[i], wls[i])
        assert tab["good_days"][i] == b.bathtub["good_days"]
        assert tab["bkg_min"][i] == b.bathtub["themin"]
        assert tab["healpix"][i] == b.healpix


def test_batch_groups_by_pixel():
    """Many sources in one field collapse to a single unique pixel."""
    _, ra, dec, _ = fields.resolve("GOODS-N")
    rng = np.random.default_rng(0)
    n = 500
    ras = ra + rng.normal(0, 0.02, n) / np.cos(np.radians(dec))
    decs = dec + rng.normal(0, 0.02, n)
    tab = jbt.get_backgrounds(ras, decs, 4.0)
    assert tab["good_days"].size == n
    assert np.unique(tab["healpix"]).size <= 2


def test_memory_cache_reuse():
    """Second call for the same pixel reuses the parsed dict object."""
    jbt.clear_cache()
    _, ra, dec, _ = fields.resolve("COSMOS")
    a = jbt.background(ra, dec, 5.0)
    b = jbt.background(ra, dec, 9.0)
    assert a.bkg_data is b.bkg_data        # same cached object, no re-parse


def test_scalar_wavelength_broadcast():
    out = jbt.get_backgrounds([189.2286, 53.1228], [62.2389, -27.8051], 5.0)
    assert np.all(out["wavelength"] == 5.0)


# --- snap-to-nearest-baked-pixel ---------------------------------------- #
def test_snap_off_by_default_keeps_exact():
    """At a baked field centre, snapping never triggers (exact pixel is baked)."""
    _, ra, dec, _ = fields.resolve("GOODS-N")
    bkg = jbt.background(ra, dec, 5.0, snap_deg=1.0)
    assert bkg.snapped is None
    assert bkg.used_healpix == bkg.healpix


def test_snap_uses_nearest_baked_offline():
    """A point ~0.7 deg off GOODS-N snaps to a baked GOODS-N pixel (no network)."""
    _, ra, dec, _ = fields.resolve("GOODS-N")
    with pytest.warns(UserWarning):
        bkg = jbt.background(ra, dec + 0.7, 5.0, snap_deg=1.0)
    assert bkg.snapped is not None
    assert bkg.snapped["field"] == "GOODS-N"
    assert 0.0 < bkg.snapped["sep_deg"] <= 1.0
    assert bkg.used_healpix != bkg.healpix
    assert bkg.bkg_data["calendar"].size > 0


def test_snap_respects_radius():
    """nearest_baked reports a real separation; tiny snap_deg does not snap."""
    _, ra, dec, _ = fields.resolve("GOODS-N")
    p, sep, field = jbt.nearest_baked(ra, dec + 0.7)
    assert field == "GOODS-N" and sep > 0
    px, info = jbt._resolve_pixel(ra, dec + 0.7, snap_deg=0.01)
    assert info is None              # too far to snap at 0.01 deg


def test_batch_snap_columns():
    """get_backgrounds reports snap separations per source when snap_deg>0."""
    _, ra, dec, _ = fields.resolve("GOODS-N")
    out = jbt.get_backgrounds([ra, ra], [dec, dec + 0.7], 5.0, snap_deg=1.0)
    assert out["snap_sep_deg"][0] == 0.0            # centre: exact
    assert out["snap_sep_deg"][1] > 0.0             # offset: snapped
    assert out["snap_field"][1] == "GOODS-N"
