"""Audit DOI-bearing references in the rejected manuscript against Crossref.

The output is diagnostic: it never silently edits citations or substitutes
records.  Network failures are reported as unverified rather than as matches.
"""
from pathlib import Path
import json
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

from docx import Document
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
DOCX = Path("C:/Users/Lenovo/IBD_rebuild_runtime/source_project/reruns/"
            "rejection_audit_2026-09-10/supplied_package/current_submission/"
            "01_Manuscript_Scientific_Reports.docx")
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", re.I)


def clean_doi(value):
    return value.rstrip(".,;)")


def norm(value):
    text = unicodedata.normalize("NFKD", value or "").encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={
        "User-Agent": "IBD-methods-audit/1.0 (mailto:wangyuspine@sina.com)"})
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.load(response)["message"]


def main():
    doc = Document(DOCX)
    start = next(i for i, p in enumerate(doc.paragraphs)
                 if p.text.strip().lower() == "references")
    refs = []
    for p in doc.paragraphs[start+1:]:
        text = p.text.strip()
        m = re.match(r"^(\d+)\.\s+(.*)$", text)
        if not m:
            if refs:
                break
            continue
        number = int(m.group(1)); raw = m.group(2)
        dm = DOI_RE.search(raw)
        refs.append(dict(number=number, raw_reference=raw,
                         doi=clean_doi(dm.group(0)) if dm else ""))

    rows = []
    for item in refs:
        row = dict(item, verification_status="NO_DOI_MANUAL_SOURCE_CHECK",
                   crossref_title="", crossref_year="", crossref_journal="",
                   title_token_overlap="")
        if item["doi"]:
            try:
                meta = crossref(item["doi"])
                title = (meta.get("title") or [""])[0]
                journal = (meta.get("container-title") or [""])[0]
                year_parts = (meta.get("published-print") or
                              meta.get("published-online") or
                              meta.get("issued") or {}).get("date-parts", [[None]])
                year = year_parts[0][0]
                title_tokens = set(norm(title).split())
                raw_tokens = set(norm(item["raw_reference"]).split())
                overlap = len(title_tokens & raw_tokens) / max(1, len(title_tokens))
                row.update(crossref_title=title, crossref_year=year,
                           crossref_journal=journal, title_token_overlap=overlap,
                           verification_status=("DOI_METADATA_MATCH" if overlap >= .8
                                                else "DOI_RESOLVES_TITLE_REVIEW_REQUIRED"))
            except (urllib.error.URLError, TimeoutError, KeyError, ValueError) as exc:
                row.update(verification_status="DOI_LOOKUP_FAILED_UNVERIFIED",
                           lookup_error=repr(exc))
            time.sleep(.08)
        rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(ROUND / "tables/existing_reference_crossref_audit.csv", index=False)
    report = {
        "status": "PASS_WITH_MANUAL_ITEMS" if len(result) else "FAIL_NO_REFERENCES",
        "n_references": len(result),
        "n_with_doi": int(result.doi.ne("").sum()),
        "status_counts": result.verification_status.value_counts().to_dict(),
        "scope": "Bibliographic DOI/title metadata only; sentence-level support still requires full-text claim checking.",
    }
    (ROUND / "reports/existing_reference_crossref_audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
