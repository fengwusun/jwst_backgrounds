"""
Registry of JWST extragalactic deep / calibration fields, with forgiving lookup.

Each entry is ``(RA, DEC, radius_deg)`` in J2000 decimal degrees, where
``radius`` is the disc (around the field centre) that the pre-baked background
cache covers so any pointing in the survey footprint is served locally.

:func:`resolve` is case-, space-, and punctuation-insensitive, honours a large
alias table, and -- when there is no exact match -- falls back to glob patterns
(``"GOODS*N"``, ``"MACS*"``) and fuzzy/typo matching, raising on ambiguity.  So
most reasonable names for a field just work, e.g. all of ``GOODS-S``, ``CDF-S``,
``HUDF``, ``"JADES deep field"`` resolve to the same field.

These coordinates are convenience caching centres, not official pointing
definitions -- always confirm exact pointings with APT.
"""

import fnmatch
import difflib

# name -> (ra_deg, dec_deg, cache_radius_deg)
DEEP_FIELDS = {
    # --- wide extragalactic survey / legacy fields ---
    "GOODS-N":   (189.228621,  62.238884, 0.20),   # JADES / FRESCO / CONGRESS / HDF-N
    "GOODS-S":   ( 53.122751, -27.805089, 0.20),   # JADES / CDF-S / HUDF / NGDEEP
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

# Any spelling -> canonical key.  Normalisation (see _norm) strips case and every
# non-alphanumeric character, so "CDF-S", "CDF S", "cdfs" are equivalent and
# multi-word names like "JADES deep field" match without separate entries.
ALIASES = {
    # GOODS-N
    "GOODS N": "GOODS-N", "GOODS NORTH": "GOODS-N", "GDN": "GOODS-N",
    "JADES-GN": "GOODS-N", "JADES GOODS-N": "GOODS-N",
    "HDF-N": "GOODS-N", "HDFN": "GOODS-N", "HUBBLE DEEP FIELD": "GOODS-N",
    "CDF-N": "GOODS-N", "CDFN": "GOODS-N", "CHANDRA DEEP FIELD NORTH": "GOODS-N",
    "CONGRESS": "GOODS-N", "FRESCO-N": "GOODS-N",
    # GOODS-S (incl. HUDF / CDF-S / NGDEEP -- same field for background purposes)
    "GOODS S": "GOODS-S", "GOODS SOUTH": "GOODS-S", "GDS": "GOODS-S",
    "JADES-GS": "GOODS-S", "JADES GOODS-S": "GOODS-S",
    "JADES DEEP FIELD": "GOODS-S", "JADES DEEP": "GOODS-S",
    "CDF-S": "GOODS-S", "CDFS": "GOODS-S", "ECDFS": "GOODS-S",
    "CHANDRA DEEP FIELD SOUTH": "GOODS-S",
    "HUDF": "GOODS-S", "UDF": "GOODS-S", "XDF": "GOODS-S", "HXDF": "GOODS-S",
    "HUDF09": "GOODS-S", "HUDF12": "GOODS-S", "NGDEEP": "GOODS-S",
    "HUBBLE ULTRA DEEP FIELD": "GOODS-S", "ULTRA DEEP FIELD": "GOODS-S",
    # COSMOS
    "COSMOS-WEB": "COSMOS", "COSMOS FIELD": "COSMOS", "PRIMER-COSMOS": "COSMOS",
    # UDS
    "SXDS": "UDS", "PRIMER-UDS": "UDS", "XMM-UDS": "UDS",
    "SUBARU XMM DEEP SURVEY": "UDS",
    # EGS
    "AEGIS": "EGS", "CEERS": "EGS", "EXTENDED GROTH STRIP": "EGS", "GROTH": "EGS",
    # NEP
    "NEP": "NEP-TDF", "TDF": "NEP-TDF", "JWST-NEP": "NEP-TDF",
    "NORTH ECLIPTIC POLE": "NEP-TDF", "TIME DOMAIN FIELD": "NEP-TDF",
    # clusters
    "A2744": "Abell2744", "ABELL 2744": "Abell2744", "UNCOVER": "Abell2744",
    "GLASS": "Abell2744", "PANDORA": "Abell2744", "PANDORAS CLUSTER": "Abell2744",
    "A370": "Abell370", "ABELL 370": "Abell370",
    "MACS J0416.1-2403": "MACS0416", "MACSJ0416": "MACS0416",
    "MACS J0717.5+3745": "MACS0717", "MACSJ0717": "MACS0717",
    "MACS J1149.5+2223": "MACS1149", "MACSJ1149": "MACS1149",
    "ABELL S1063": "AS1063", "RXC J2248.7-4431": "AS1063", "RXJ2248": "AS1063",
    "MACS J1423.8+2404": "MACS1423", "MACSJ1423": "MACS1423",
    "ACT-CL J0102-4915": "ElGordo", "EL GORDO": "ElGordo", "ACTJ0102": "ElGordo",
    "SMACS J0723.3-7327": "SMACS0723", "SMACS": "SMACS0723", "ERO": "SMACS0723",
    "WEBBS FIRST DEEP FIELD": "SMACS0723", "FIRST DEEP FIELD": "SMACS0723",
    # Magellanic
    "LARGE MAGELLANIC CLOUD": "LMC", "SMALL MAGELLANIC CLOUD": "SMC",
}


def _norm(name):
    """Uppercase, alphanumeric-only key (drops spaces, dashes, dots, etc.)."""
    return "".join(ch for ch in str(name).upper() if ch.isalnum())


def _norm_glob(name):
    """Like _norm but keep glob metacharacters so patterns survive."""
    return "".join(ch for ch in str(name).upper() if ch.isalnum() or ch in "*?[]")


def _table():
    """Mapping of normalised key -> canonical name for all names and aliases."""
    t = {_norm(k): k for k in DEEP_FIELDS}
    for alias, canon in ALIASES.items():
        t[_norm(alias)] = canon
    return t


def _pick(hits, raw):
    hits = sorted(set(hits))
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise KeyError("Field %r is ambiguous; matches: %s. Please be more specific."
                       % (raw, ", ".join(hits)))
    raise KeyError("Unknown field %r. Known fields: %s (+ many aliases; "
                   "wildcards like 'MACS*' also work)." % (raw, ", ".join(list_fields())))


def resolve(name):
    """Return ``(canonical_name, ra, dec, radius_deg)`` for a field.

    Accepts canonical names, aliases, glob patterns ('GOODS*N', 'MACS*'), and
    near-miss spellings.  Raises ``KeyError`` (listing candidates) if the name is
    ambiguous or unknown.
    """
    raw = str(name).strip()
    table = _table()
    key = _norm(raw)

    if key in table:                                   # exact (normalised) match
        canon = table[key]
    elif any(c in raw for c in "*?["):                 # glob pattern
        pat = _norm_glob(raw)
        canon = _pick((c for k, c in table.items() if fnmatch.fnmatchcase(k, pat)), raw)
    else:                                              # substring, then fuzzy
        hits = [c for k, c in table.items()
                if len(key) >= 4 and len(k) >= 4 and (key in k or k in key)]
        if not hits:
            close = difflib.get_close_matches(key, list(table), n=5, cutoff=0.82)
            hits = [table[k] for k in close]
        canon = _pick(hits, raw)

    ra, dec, rad = DEEP_FIELDS[canon]
    return canon, ra, dec, rad


def search(pattern):
    """Return all canonical field names matching a glob or substring (sorted)."""
    table = _table()
    if any(c in pattern for c in "*?["):
        pat = _norm_glob(pattern)
        hits = (c for k, c in table.items() if fnmatch.fnmatchcase(k, pat))
    else:
        key = _norm(pattern)
        hits = (c for k, c in table.items() if key and (key in k or k in key))
    return sorted(set(hits))


def list_fields():
    """Return the sorted list of canonical field names."""
    return sorted(DEEP_FIELDS)
