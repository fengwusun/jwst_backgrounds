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
    # --- wide extragalactic survey / legacy fields ---
    "GOODS-N":   (189.228621,  62.238884, 0.20),   # JADES / FRESCO / CONGRESS
    "GOODS-S":   ( 53.122751, -27.805089, 0.20),   # JADES / CDF-S
    "HUDF":      ( 53.162500, -27.791400, 0.10),   # Hubble Ultra Deep Field / XDF / NGDEEP
    "COSMOS":    (150.116321,   2.207641, 0.55),   # COSMOS-Web is ~0.54 deg^2
    "UDS":       ( 34.349577,  -5.201059, 0.30),   # SXDS / PRIMER-UDS
    "EGS":       (214.825000,  52.825000, 0.60),   # AEGIS / CEERS strip is long
    "NEP-TDF":   (260.699600,  65.822500, 0.25),   # JWST N. Ecliptic Pole TDF
    # --- lensing-cluster fields (high-z behind massive clusters) ---
    "Abell2744": (  3.586250, -30.400190, 0.25),   # UNCOVER / GLASS + parallels
    "Abell370":  ( 39.968000,  -1.586500, 0.20),   # Frontier Field
    "MACS0416":  ( 64.038100, -24.067500, 0.15),   # MACS J0416.1-2403 (Frontier Field)
    "MACS0717":  (109.383800,  37.755800, 0.15),   # MACS J0717.5+3745 (Frontier Field)
    "MACS1149":  (177.398700,  22.398500, 0.15),   # MACS J1149.5+2223 (Frontier Field)
    "AS1063":    (342.183200, -44.530900, 0.15),   # Abell S1063 / RXC J2248.7-4431
    "MACS1423":  (215.949200,  24.078300, 0.15),   # MACS J1423.8+2404 (CANUCS)
    "ElGordo":   ( 15.729200, -49.254700, 0.15),   # ACT-CL J0102-4915
    "SMACS0723": (110.831300, -73.454300, 0.15),   # SMACS J0723.3-7327 (JWST ERO)
    # --- Magellanic calibration fields ---
    "LMC":       ( 80.894167, -69.756111, 0.20),   # LMC astrometric field
    "SMC":       ( 13.158300, -72.800300, 0.20),   # SMC
}

# alternative spellings -> canonical key in DEEP_FIELDS
ALIASES = {
    "GOODSN": "GOODS-N", "GOODS_N": "GOODS-N", "JADES-GN": "GOODS-N",
    "HDFN": "GOODS-N", "GDN": "GOODS-N",
    "GOODSS": "GOODS-S", "GOODS_S": "GOODS-S", "JADES-GS": "GOODS-S",
    "CDFS": "GOODS-S", "CDF-S": "GOODS-S", "ECDFS": "GOODS-S", "GDS": "GOODS-S",
    "XDF": "HUDF", "UDF": "HUDF", "HXDF": "HUDF", "NGDEEP": "HUDF",
    "COSMOS-WEB": "COSMOS", "COSMOSWEB": "COSMOS", "PRIMER-COSMOS": "COSMOS",
    "SXDS": "UDS", "PRIMER-UDS": "UDS",
    "AEGIS": "EGS", "CEERS": "EGS",
    "A2744": "Abell2744", "ABELL-2744": "Abell2744", "UNCOVER": "Abell2744",
    "A370": "Abell370", "ABELL-370": "Abell370",
    "MACSJ0416": "MACS0416", "MACS-0416": "MACS0416", "MACSJ0416.1-2403": "MACS0416",
    "MACSJ0717": "MACS0717", "MACS-0717": "MACS0717",
    "MACSJ1149": "MACS1149", "MACS-1149": "MACS1149",
    "ABELLS1063": "AS1063", "RXCJ2248": "AS1063", "RXJ2248": "AS1063",
    "MACSJ1423": "MACS1423", "MACS-1423": "MACS1423",
    "ACTJ0102": "ElGordo", "ACT-CLJ0102-4915": "ElGordo", "EL-GORDO": "ElGordo",
    "SMACSJ0723": "SMACS0723", "SMACS-0723": "SMACS0723", "SMACS": "SMACS0723", "ERO": "SMACS0723",
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
