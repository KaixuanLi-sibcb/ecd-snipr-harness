"""Lossless workbook fields with conservative, versioned layout interpretation."""

import re
import unicodedata
from collections import Counter
from .common import digest

DICTIONARY_VERSION = "1.0"
ALIASES = {
    "AgID": "agid", "gene": "gene", "Entry": "accession",
    "Ag-SNIPR": "receiver_label", "Ag-SNPIR": "receiver_label",
    "Ag-SNIPR plasmid": "receiver_plasmid", "Ag-SNPIR plasmid": "receiver_plasmid",
    "SNIPR Assay(Y/N)": "snipr_assay_status", "SNPIR Assay(Y/N)": "snipr_assay_status",
    "self-activation of Ag-SNIPR": "self_activation",
    "self-activation of Ag-SNPIR": "self_activation",
    "Surface_level(MFI)": "surface_mfi",
    "BFP+/myc+ (MFI)": "bfp_myc_mfi", "BFP+/myc+ | Freq. of Parent": "bfp_myc_frequency",
    "BFP+/myc+/mRuby3+": "bfp_myc_mruby_unspecified",
    "BFP+/myc+/mRuby3+ | Freq. of Parent": "bfp_myc_mruby_frequency",
    "Single Cells/BFP+": "bfp_unspecified",
    "Single Cells/BFP+ | Freq. of Parent": "bfp_frequency",
    "ECD_aa": "ecd_aa", "ECD_nt": "ecd_nt", "ECD": "ecd_unspecified",
    "length ECD_aa": "ecd_length_aa", "length ECD_nt": "ecd_length_nt",
    "length ECD": "ecd_length_unspecified", "ECD_identity": "ecd_identity",
    "FULL_nt": "full_nt", "loc.": "location", "note": "note", "transmem_note": "topology_note",
    "mRNA-LNP": "mrna_lnp_label", "number (mRNA-LNP)": "mrna_lnp_label",
    "mRNA-LNP(Y/N)": "mrna_lnp_status", "mRNA plasmid (Y/N)": "mrna_plasmid_status",
    "Ag mRNA plasmid": "mrna_plasmid_label", "antigen": "antigen_label",
    "number": "number_label", "Sample:": "sample_label", "BCID": "bcid",
    "Construct engineering strategy": "engineering_strategy", "name": "name",
    "Source": "source_label", "Serological level(EC50)": "serology_ec50",
    "Cross-reference (CCDS)": "ccds", "state": "state", "Resistance": "resistance",
}


def header_key(value):
    text = unicodedata.normalize("NFKC", str(value or "")).lower()
    text = text.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-")
    return re.sub(r"[\s_]", "", text)


LOOKUP = {header_key(k): v for k, v in ALIASES.items()}
ASSAY_FIELDS = {v for v in ALIASES.values() if v.startswith("bfp_")} | {
    "self_activation", "surface_mfi", "snipr_assay_status", "serology_ec50"}
IDENTITY_FIELDS = {"agid", "gene", "accession", "receiver_label", "antigen_label", "sample_label"}
ACCESSION = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?")


def nonblank(value):
    return value is not None and str(value).strip() != ""


def column(coordinate):
    return re.sub(r"[0-9]+$", "", coordinate)


def field_values(row, field):
    return [f["raw_value"] for f in row["fields"] if f["field"] == field and nonblank(f["raw_value"])]


def extract(data, layout=None):
    """Map cells, never collapse duplicate columns or fill merged cells implicitly."""
    layout = layout or {}
    regions = layout.get("regions", [])
    if layout and layout.get("source_sha256") != data["source_sha256"]:
        raise ValueError("Layout workbook hash mismatch")
    sheets = {s["name"] for s in data["sheets"]}
    for i, region in enumerate(regions):
        if (region.get("sheet") not in sheets or not region.get("rationale") or
                type(region.get("start_row")) is not int or type(region.get("end_row")) is not int or
                not 1 <= region["start_row"] <= region["end_row"] or not region.get("columns")):
            raise ValueError("Layout needs exact sheet, row interval, columns and rationale")
        if any(not re.fullmatch("[A-Z]+", c) or f not in set(ALIASES.values()) for c, f in region["columns"].items()):
            raise ValueError("Unknown explicit layout column/field")
        if any(r["sheet"] == region["sheet"] and r["start_row"] <= region["end_row"] and
               region["start_row"] <= r["end_row"] for r in regions[:i]):
            raise ValueError("Overlapping layout regions")
    raw_cells, records, headers = [], [], []
    for sheet in data["sheets"]:
        active, header_row = {}, None
        for row in sheet["rows"]:
            detected = {column(c["coordinate"]): {"field": LOOKUP.get(header_key(c["value"]), "unmapped"),
                        "raw_header": c["value"], "header_cell": c["coordinate"]}
                        for c in row["cells"] if nonblank(c["value"])}
            known = {c["field"] for c in detected.values()} - {"unmapped"}
            is_header = len(known) >= 3 and bool(known & IDENTITY_FIELDS)
            if is_header:
                active, header_row = detected, row["row"]
                headers.append({"sheet": sheet["name"], "row": header_row, "columns": active})
            override = next((r for r in regions if r["sheet"] == sheet["name"] and
                             r["start_row"] <= row["row"] <= r["end_row"]), None)
            mapping = ({c: {"field": f, "raw_header": None, "header_cell": None}
                        for c, f in override["columns"].items()} if override else active)
            row_id = digest([data["source_sha256"], sheet["name"], row["row"]])
            source = {"workbook": data["source_name"], "workbook_sha256": data["source_sha256"],
                      "sheet": sheet["name"], "row": row["row"]}
            values = {column(c["coordinate"]): c for c in row["cells"]}
            fields = []
            for col in dict.fromkeys(list(values) + [c for c, h in mapping.items() if h["field"] in ASSAY_FIELDS]):
                cell = values.get(col, {"coordinate": f"{col}{row['row']}", "value": None, "formula": None})
                h = mapping.get(col, {"field": "unmapped", "raw_header": None, "header_cell": None})
                item = dict(source, cell_id=digest([data["source_sha256"], sheet["name"], cell["coordinate"]]),
                            row_id=row_id, cell=cell["coordinate"], raw_value=cell["value"],
                            raw_formula=cell.get("formula"), cell_type=cell.get("cell_type"),
                            formula_present=cell.get("formula_present", cell.get("formula") is not None),
                            formula_attributes=cell.get("formula_attributes", {}),
                            shared_formula_anchors=cell.get("shared_formula_anchors", []),
                            style_index=cell.get("style_index"), header_row=header_row,
                            **h, mapping_basis="explicit_layout" if override else "header_dictionary",
                            cell_role="header" if is_header else "data", present_in_source=col in values)
                fields.append(item)
                if col in values:
                    raw_cells.append(item)
            if is_header or not any(nonblank(c["value"]) or c.get("formula_present") or c.get("formula") for c in row["cells"]):
                continue
            record = dict(source, row_id=row_id, header_row=header_row, fields=fields,
                          layout_status="mapped" if mapping else "no_header", issues=[])
            accessions = field_values(record, "accession")
            elsewhere = [c for c in row["cells"] if ACCESSION.fullmatch(str(c["value"] or "").strip())
                         and mapping.get(column(c["coordinate"]), {}).get("field") != "accession"]
            if (accessions and any(not ACCESSION.fullmatch(str(a).strip()) for a in accessions)) or (elsewhere and not accessions):
                record["layout_status"] = "suspect_or_composite_identity"
                record["issues"].append("identity_or_shifted_layout_requires_review")
            records.append(record)
    counts = Counter(c["field"] for c in raw_cells if c["cell_role"] == "data" and nonblank(c["raw_value"]))
    return {"rows": records, "cells": raw_cells, "headers": headers,
            "counts": {"source_cells": len(raw_cells), "source_sheets": len(data["sheets"]),
                       "data_rows": len(records), "header_sections": len(headers),
                       "nonblank_data_cells_by_field": dict(sorted(counts.items())),
                       "layout_status_rows": dict(Counter(r["layout_status"] for r in records))}}
