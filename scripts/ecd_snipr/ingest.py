"""Read-only tabular inventory and explicit field mapping; no biological inference."""

import csv
import posixpath
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from .common import file_hash

NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def inventory(path):
    path = Path(path)
    result = {"source_name": path.name, "source_sha256": file_hash(path), "sheets": [],
              "interpretation": "raw cells only; blanks and Y/N are not outcome labels"}
    if path.suffix.lower() in {".tsv", ".csv"}:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            values = list(csv.reader(handle, delimiter="\t" if path.suffix.lower() == ".tsv" else ","))
        rows = [{"row": i+1, "cells": [{"coordinate": f"{column_name(j+1)}{i+1}", "value": v, "formula": None} for j, v in enumerate(row)]} for i, row in enumerate(values)]
        result["sheets"].append({"name": path.stem, "rows": rows, "merged_cells": []})
        return result
    if path.suffix.lower() != ".xlsx":
        raise ValueError("Supported read-only inputs: XLSX, TSV, CSV")
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            tree = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = ["".join(t.text or "" for t in si.findall(".//s:t", NS)) for si in tree.findall("s:si", NS)]
        relations = {n.attrib["Id"]: n.attrib["Target"] for n in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))}
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        for sheet in workbook.findall("s:sheets/s:sheet", NS):
            target = relations[sheet.attrib[REL]]
            target = target.lstrip("/") if target.startswith("/") else posixpath.normpath("xl/" + target)
            tree = ET.fromstring(archive.read(target))
            rows = []
            for r in tree.findall("s:sheetData/s:row", NS):
                cells = []
                for cell in r.findall("s:c", NS):
                    v = cell.find("s:v", NS)
                    f = cell.find("s:f", NS)
                    value = v.text if v is not None else None
                    if cell.get("t") == "s" and value is not None:
                        value = strings[int(value)]
                    elif cell.get("t") == "inlineStr":
                        value = "".join(n.text or "" for n in cell.findall(".//s:t", NS))
                    cells.append({"coordinate": cell.attrib["r"], "value": value,
                                  "formula": f.text if f is not None else None,
                                  "cell_type": cell.get("t", "n"), "style_index": cell.get("s", "")})
                rows.append({"row": int(r.attrib["r"]), "cells": cells})
            result["sheets"].append({"name": sheet.attrib["name"], "rows": rows,
                                     "merged_cells": [m.attrib["ref"] for m in tree.findall("s:mergeCells/s:mergeCell", NS)]})
    return result


def column_name(number):
    result = ""
    while number:
        number, rem = divmod(number-1, 26)
        result = chr(65+rem) + result
    return result


def map_rows(data, mapping):
    sheet = next((s for s in data["sheets"] if s["name"] == mapping["sheet"]), None)
    if sheet is None:
        raise ValueError("Exact sheet name not found (including trailing spaces)")
    if type(mapping.get("header_row")) is not int:
        raise ValueError("header_row required")
    fields = mapping.get("fields", {})
    if not fields or any(not re.fullmatch("[A-Z]+", c) for c in fields.values()):
        raise ValueError("Explicit field -> Excel column letters mapping required")
    result = []
    for row in sheet["rows"]:
        if row["row"] <= mapping["header_row"]:
            continue
        cells = {re.sub(r"\d+$", "", c["coordinate"]): c for c in row["cells"]}
        result.append({"values": {k: cells.get(c, {}).get("value") for k, c in fields.items()},
                       "source": {"workbook": data["source_name"], "sha256": data["source_sha256"], "sheet": sheet["name"], "row": row["row"],
                                  "cells": {k: cells.get(c, {"coordinate": f"{c}{row['row']}", "value": None}) for k, c in fields.items()}},
                       "interpretation": "unconfirmed; explicit semantic contract required before assay normalization"})
    return result
