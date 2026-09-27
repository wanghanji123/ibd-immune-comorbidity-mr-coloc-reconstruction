"""Select one Zotero-exported BibTeX record per audited DOI.

Records remain Zotero-derived.  Duplicate DOI entries are ranked by metadata
completeness; unsupported/manual web records are not invented here.
"""
from pathlib import Path
import csv
import re

import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROUND / "manuscript_rebuild_2026-09-28"
FULL_BIB = MANUSCRIPT / "zotero_library_export.bib"
AUDIT = ROUND / "tables/existing_reference_crossref_audit.csv"


def entries(text):
    starts = [m.start() for m in re.finditer(r"(?m)^@", text)]
    for i, start in enumerate(starts):
        end = starts[i+1] if i+1 < len(starts) else len(text)
        block = text[start:end].strip()
        key = re.match(r"@\w+\{([^,]+),", block).group(1)
        doi = re.search(r"(?im)^\s*doi\s*=\s*[\{\"]([^}\"]+)", block)
        yield key, (doi.group(1).lower().strip() if doi else ""), block


def score(block):
    wanted = ["title", "author", "year", "journal", "volume", "number", "pages", "doi"]
    return sum(bool(re.search(rf"(?im)^\s*{field}\s*=", block)) for field in wanted), len(block)


def main():
    audited = pd.read_csv(AUDIT).fillna("")
    wanted = {d.lower() for d in audited.doi if d}
    full = FULL_BIB.read_text(encoding="utf-8")
    by_doi = {}
    for key, doi, block in entries(full):
        if doi in wanted:
            by_doi.setdefault(doi, []).append((key, block))

    selected = []
    mapping = []
    for doi in sorted(wanted):
        candidates = by_doi.get(doi, [])
        if not candidates:
            mapping.append(dict(doi=doi, zotero_bibtex_key="", n_duplicates=0,
                                status="NOT_FOUND_IN_ZOTERO_EXPORT"))
            continue
        key, block = max(candidates, key=lambda x: score(x[1]))
        selected.append(block)
        mapping.append(dict(doi=doi, zotero_bibtex_key=key,
                            n_duplicates=len(candidates), status="SELECTED"))

    (MANUSCRIPT / "references_curated_from_zotero.bib").write_text(
        "\n\n".join(selected) + "\n", encoding="utf-8")
    pd.DataFrame(mapping).to_csv(
        ROUND / "tables/zotero_reference_key_map.csv", index=False,
        quoting=csv.QUOTE_MINIMAL)
    print({"requested_dois": len(wanted), "selected": len(selected),
           "missing": sum(x["status"] != "SELECTED" for x in mapping),
           "duplicate_doi_records": sum(max(0, x["n_duplicates"]-1) for x in mapping)})


if __name__ == "__main__":
    main()
