"""Small deterministic primitives shared by the workflow."""

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

AA = set("ACDEFGHIKLMNPQRSTVWY")


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def sequence(value):
    if not isinstance(value, str):
        raise ValueError("Sequence must be a string")
    result = re.sub(r"\s+", "", value).upper()
    if not result or set(result) - AA:
        raise ValueError("Exact sequence required: no gaps, stops, X or nonstandard residues")
    return result


def interval(feature, length):
    start, end = feature.get("start"), feature.get("end")
    if feature.get("boundary_status", "exact") != "exact":
        raise ValueError("Fuzzy boundaries require review and an exact replacement annotation")
    if type(start) is not int or type(end) is not int or not 1 <= start <= end <= length:
        raise ValueError("Coordinates must be 1-based inclusive within the reference sequence")
    return start, end


def overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def contains(a, b):
    return a[0] <= b[0] and b[1] <= a[1]


def sourced(value):
    e = value.get("evidence", {})
    return bool(e.get("source") not in {None, "", "unknown", "pending"} and e.get("version") not in {None, "", "unknown", "pending"} and e.get("kind") in {
        "curated_annotation", "prediction", "local_experiment", "synthetic_fixture"
    })


def tsv(path, rows, columns):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: canonical(v) if isinstance(v, (dict, list)) else v
                             for k, v in row.items()})


def unique(items, key):
    values = [i[key] for i in items]
    if len(set(values)) != len(values):
        raise ValueError(f"Duplicate {key}: identities must not be silently merged")


def translate(dna):
    dna = re.sub(r"\s+", "", dna).upper()
    if not dna or set(dna) - set("ACGT") or len(dna) % 3:
        raise ValueError("CDS must be unambiguous and divisible by three")
    bases = "TCAG"
    letters = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = [a+b+c for a in bases for b in bases for c in bases]
    code = dict(zip(codons, letters))
    protein = "".join(code[dna[i:i+3]] for i in range(0, len(dna), 3))
    if "*" in protein.rstrip("*") or protein.endswith("**"):
        raise ValueError("Internal stop in CDS")
    return protein.removesuffix("*")
