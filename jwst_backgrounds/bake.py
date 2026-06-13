"""
Pre-bake the STScI background cache for named JWST fields.

For each field this downloads every healpix cache file covering a disc around the
field centre and stores it (byte-for-byte, so output stays identical to the
official package) in ``jwst_backgrounds/field_cache/``, together with a
``manifest.json`` describing what was baked.

Re-bake the shipped default set::

    python -m jwst_backgrounds.bake

Bake a custom set / cache / radius::

    python -m jwst_backgrounds.bake --fields GOODS-N COSMOS --radius 0.3 --cache 2.0
"""

import os
import json
import argparse
import urllib.request

import numpy as np
import healpy

from jwst_backgrounds import fields as _fields
from jwst_backgrounds.jbt import (NSIDE, CACHE_URLS, DEFAULT_CACHE_URL, _FIELD_CACHE,
                                  file_from_healpix, get_cache_version)


def pixels_for(ra, dec, radius_deg):
    """All healpix RING pixels overlapping a disc of ``radius_deg`` about (ra, dec)."""
    vec = healpy.ang2vec(ra, dec, lonlat=True)
    px = healpy.query_disc(NSIDE, vec, np.radians(radius_deg), inclusive=True, nest=False)
    return sorted(int(p) for p in px)


def bake(field_names=None, cache_url=DEFAULT_CACHE_URL, out_dir=_FIELD_CACHE,
         radius_override=None, verbose=True):
    """Download and store the cache files covering the requested fields.

    Parameters
    ----------
    field_names : list of str or None
        Fields to bake; ``None`` bakes the whole :data:`fields.DEEP_FIELDS` set.
    cache_url : str
        STScI cache to bake from (default = the one the official package uses).
    out_dir : str
        Destination directory for the ``.bin`` files and ``manifest.json``.
    radius_override : float or None
        If set, use this disc radius (deg) for every field instead of the
        per-field default in the registry.
    """
    if field_names is None:
        field_names = _fields.list_fields()

    os.makedirs(out_dir, exist_ok=True)
    version = get_cache_version(cache_url)

    manifest = {"cache_url": cache_url, "version": version, "nside": NSIDE, "fields": {}}
    all_pixels = set()
    for name in field_names:
        canon, ra, dec, rad = _fields.resolve(name)
        if radius_override is not None:
            rad = radius_override
        px = pixels_for(ra, dec, rad)
        manifest["fields"][canon] = {"ra": ra, "dec": dec, "radius_deg": rad, "pixels": px}
        all_pixels.update(px)
        if verbose:
            print("%-10s (%.4f, %.4f) r=%.2f deg -> %d pixel(s)" % (canon, ra, dec, rad, len(px)))

    total = 0
    for i, p in enumerate(sorted(all_pixels), 1):
        dest = os.path.join(out_dir, "sl_pix_%06d.bin" % p)
        if os.path.isfile(dest):
            total += os.path.getsize(dest)
            continue
        url = cache_url + file_from_healpix(p)
        with urllib.request.urlopen(url, timeout=120) as fh:
            data = fh.read()
        tmp = dest + ".part"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, dest)
        total += len(data)
        if verbose:
            print("  [%d/%d] sl_pix_%06d.bin  (%.0f KB)" % (i, len(all_pixels), p, len(data) / 1024))

    manifest["n_pixels"] = len(all_pixels)
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)

    if verbose:
        print("\nBaked %d field(s), %d unique pixel(s), %.1f MB into %s"
              % (len(field_names), len(all_pixels), total / 1024 / 1024, out_dir))
    return manifest


def main():
    ap = argparse.ArgumentParser(description="Pre-bake JWST background cache for named fields.")
    ap.add_argument("--fields", nargs="+", default=None,
                    help="field names to bake (default: all known fields)")
    ap.add_argument("--cache", default="2.0", choices=sorted(CACHE_URLS),
                    help="STScI cache version to bake from (default: 2.0, matches the package)")
    ap.add_argument("--radius", type=float, default=None,
                    help="override disc radius in degrees for every field")
    ap.add_argument("--out", default=_FIELD_CACHE, help="output directory")
    args = ap.parse_args()
    bake(args.fields, cache_url=CACHE_URLS[args.cache], out_dir=args.out,
         radius_override=args.radius)


if __name__ == "__main__":
    main()
