"""Audit the final Zotero-derived bibliography against Crossref metadata."""
from pathlib import Path
import json
import re
import subprocess
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
M = ROUND / "manuscript_rebuild_2026-09-28"


def norm(value):
    text = unicodedata.normalize("NFKD", value or "").encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "IBD-round03-reference-audit/1.0 (mailto:wangyuspine@sina.com)"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)["message"]


def main():
    csl_json = M / "references_submission.csl.json"
    subprocess.run(
        ["pandoc", str(M / "references_submission.bib"), "-f", "biblatex",
         "-t", "csljson", "-o", str(csl_json)], check=True
    )
    records = json.loads(csl_json.read_text(encoding="utf-8"))
    manuscript = (M / "manuscript_source.md").read_text(encoding="utf-8")
    cited = set(re.findall(r"@([A-Za-z0-9_:.\-]+)", manuscript))

    rows = []
    for item in records:
        doi = (item.get("DOI") or "").strip().lower()
        local_title = item.get("title") or ""
        row = {
            "citekey": item["id"],
            "cited_in_manuscript": item["id"] in cited,
            "doi": doi,
            "local_title": local_title,
            "local_year": ((item.get("issued") or {}).get("date-parts") or [[None]])[0][0],
            "verification_status": "NO_DOI_UNVERIFIED",
            "crossref_title": "",
            "crossref_year": "",
            "crossref_journal": "",
            "title_token_overlap": "",
            "lookup_error": "",
        }
        if doi:
            try:
                meta = crossref(doi)
                title = (meta.get("title") or [""])[0]
                journal = (meta.get("container-title") or [""])[0]
                date = meta.get("published-print") or meta.get("published-online") or meta.get("issued") or {}
                year = (date.get("date-parts") or [[None]])[0][0]
                title_tokens = set(norm(title).split())
                local_tokens = set(norm(local_title).split())
                overlap = len(title_tokens & local_tokens) / max(1, len(title_tokens))
                row.update(
                    crossref_title=title,
                    crossref_year=year,
                    crossref_journal=journal,
                    title_token_overlap=overlap,
                    verification_status=("DOI_METADATA_MATCH" if overlap >= 0.8
                                         else "DOI_RESOLVES_TITLE_REVIEW_REQUIRED"),
                )
            except (urllib.error.URLError, TimeoutError, KeyError, ValueError) as exc:
                row.update(verification_status="DOI_LOOKUP_FAILED_UNVERIFIED", lookup_error=repr(exc))
            time.sleep(0.08)
        rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(ROUND / "tables/final_reference_crossref_audit.csv", index=False)
    missing_keys = sorted(cited - set(result.citekey))
    duplicate_dois = sorted(result.loc[result.doi.ne("") & result.doi.duplicated(False), "doi"].unique())
    report = {
        "status": "PASS" if not missing_keys and not duplicate_dois else "FAIL",
        "n_bibliography_records": len(result),
        "n_cited_records": int(result.cited_in_manuscript.sum()),
        "n_unique_cited_keys_in_manuscript": len(cited),
        "n_unique_doi": int(result.doi.replace("", pd.NA).nunique()),
        "missing_citekeys": missing_keys,
        "duplicate_dois": duplicate_dois,
        "status_counts_all": result.verification_status.value_counts().to_dict(),
        "status_counts_cited": result.loc[result.cited_in_manuscript, "verification_status"].value_counts().to_dict(),
        "metadata_boundary": "Crossref agreement validates bibliographic identity, not whether each cited paper supports every surrounding sentence. Claim-level support remains a manuscript review task.",
        "manual_correction_log": "tables/reference_metadata_corrections.csv",
    }
    (ROUND / "reports/final_reference_crossref_audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
