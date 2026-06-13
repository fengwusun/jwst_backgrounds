1.3.0+fast (June 13, 2026)
==========================

Fast fork of jwst_backgrounds. Output is verified identical to the official
1.3.0 package (the default cache is unchanged: sl_cache_2.0).

- Pre-baked background cache shipped for JWST deep fields (GOODS-N/S, COSMOS,
  UDS, EGS, Abell 2744, Abell 370, NEP-TDF, LMC): those positions work offline.
- Disk + in-memory caching of cache files keyed by healpix pixel; VERSION is
  fetched once per process instead of on every call.
- Vectorized binary parse (numpy.frombuffer) and single-wavelength bathtub
  interpolation (no per-call interp1d rebuilds).
- New batch API jbt.get_backgrounds() that groups a catalog by healpix pixel.
- background.from_field(), a field registry (jwst_backgrounds/fields.py), and a
  re-baking tool (python -m jwst_backgrounds.bake).
- CLI gains --field, --list-fields and --cache; plot_background/plot_bathtub
  return (fig, ax) and accept an existing ax=.
- Optional access to the regenerated sl_cache_3.0 via jbt.CACHE_URLS["3.0"].

1.3.0 (July 18, 2024)
=====================

- Updates for JWST cycle 4
- Added v4.0 of version file
- CI/CD testing updates

1.2.0 (Nov 28, 2022)
====================

- Updates for JWST cycle 2
- New thermal file.
- Updated cache location for updated stray light files.

1.1.2 (Mar 25, 2020)
====================

- New thermal file was added to package. Updated documentation and applied the level 1 INS JWST Community Standard to the repository.