1.3.0+fast (June 13, 2026)
==========================

Fast fork of jwst_backgrounds. Output is verified identical to the official
1.3.0 package (the default cache is unchanged: sl_cache_2.0).

- Pre-baked background cache shipped for 17 JWST extragalactic deep,
  lensing-cluster and calibration fields (GOODS-N/S, COSMOS, UDS, EGS, NEP-TDF,
  Abell 2744/370, MACS 0416/0717/1149, AS1063, MACS 1423, El Gordo, SMACS 0723,
  LMC, SMC): those positions work offline, byte-identical to STScI.
- Forgiving field lookup: case/space/punctuation-insensitive, a large alias
  table (HUDF/CDF-S/NGDEEP/"JADES deep field" -> GOODS-S, CEERS -> EGS, ...),
  glob patterns ("MACS*"), and fuzzy/typo fallback, with disambiguation.
- Optional snap-to-nearest: an off-bundle query within snap_deg of a baked field
  reuses it (off by default; ~1% at 1 deg in the NIR, less in the MIRI range).
- Disk + in-memory caching of cache files keyed by healpix pixel; VERSION is
  read from the bundle manifest for the baked cache (no per-call network).
- Vectorized binary parse (numpy.frombuffer) and single-wavelength bathtub
  interpolation (no per-call interp1d rebuilds).
- New batch API jbt.get_backgrounds() that groups a catalog by healpix pixel.
- background.from_field(), a field registry (jwst_backgrounds/fields.py), and a
  re-baking tool (python -m jwst_backgrounds.bake).
- CLI gains --field, --list-fields, --cache and --snap; plot_background/
  plot_bathtub return (fig, ax) and accept an existing ax=.
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