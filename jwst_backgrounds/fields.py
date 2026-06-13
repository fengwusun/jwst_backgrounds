"""
Registry of JWST extragalactic deep / calibration fields.

Each entry is (RA, DEC, radius_deg) in J2000 decimal degrees, where ``radius``
is the disc (around the field centre) that the pre-baked background cache should
cover so that any pointing within the survey footprint is served locally.

Use :func:`resolve` to look up a field by name (case-insensitive, aliases
honoured).  These coordinates are convenience centres for caching, not official
pointing definitions -- always confirm exact pointings with APT.
"""

# name -> (ra_deg, dec_deg, cache_radius_deg)
DEEP_FIELDS = {
    "GOODS-N":   (189.228621,  62.238884, 0.20),
    "GOODS-S":   ( 53.122751, -27.805089, 0.20),
    "COSMOS":    (150.116321,   2.207641, 0.55),   # COSMOS-Web is ~0.54 deg^2
    "UDS":       ( 34.349577,  -5.201059, 0.30),   # SXDS / PRIMER-UDS
    "EGS":       (214.825000,  52.825000, 0.60),   # AEGIS / CEERS strip is long
    "Abell2744": (  3.586250, -30.400190, 0.25),   # UNCOVER + parallels
    "Abell370":  ( 39.968000,  -1.586500, 0.20),   # Frontier Field
    "NEP-TDF":   (260.699600,  65.822500, 0.25),   # JWST N. Ecliptic Pole TDF
    "LMC":       ( 80.894167, -69.756111, 0.20),   # LMC astrometric field
}

# alternative spellings -> canonical key in DEEP_FIELDS
ALIASES = {
    "GOODSN": "GOODS-N", "GOODS_N": "GOODS-N", "JADES-GN": "GOODS-N",
    "HDFN": "GOODS-N", "GDN": "GOODS-N",
    "GOODSS": "GOODS-S", "GOODS_S": "GOODS-S", "JADES-GS": "GOODS-S",
    "CDFS": "GOODS-S", "CDF-S": "GOODS-S", "ECDFS": "GOODS-S", "GDS": "GOODS-S",
    "COSMOS-WEB": "COSMOS", "COSMOSWEB": "COSMOS", "PRIMER-COSMOS": "COSMOS",
    "SXDS": "UDS", "PRIMER-UDS": "UDS",
    "AEGIS": "EGS", "CEERS": "EGS",
    "A2744": "Abell2744", "ABELL-2744": "Abell2744", "UNCOVER": "Abell2744",
    "A370": "Abell370", "ABELL-370": "Abell370",
    "NEP": "NEP-TDF", "TDF": "NEP-TDF", "JWST-NEP": "NEP-TDF",
}


def _norm(name):
    return str(name).strip().upper().replace(" ", "")


def resolve(name):
    """Return ``(canonical_name, ra, dec, radius_deg)`` for a field name.

    Lookup is case- and separator-insensitive and honours :data:`ALIASES`.
    Raises ``KeyError`` with the list of known fields if the name is unknown.
    """
    key = _norm(name)
    # direct (normalised) match against canonical names
    norm_map = {_norm(k): k for k in DEEP_FIELDS}
    if key in norm_map:
        canon = norm_map[key]
    elif key in {_norm(a): a for a in ALIASES}:
        canon = ALIASES[{_norm(a): a for a in ALIASES}[key]]
    else:
        raise KeyError(
            f"Unknown field {name!r}. Known fields: "
            f"{', '.join(sorted(DEEP_FIELDS))} (+ aliases)."
        )
    ra, dec, rad = DEEP_FIELDS[canon]
    return canon, ra, dec, rad


def list_fields():
    """Return the sorted list of canonical field names."""
    return sorted(DEEP_FIELDS)
