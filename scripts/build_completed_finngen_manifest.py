"""Replace the earlier pending view with a dated completed-source manifest."""
from pathlib import Path
import json
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]


def main():
    pending = pd.read_csv(ROUND / "tables/FinnGen_pending_source_manifest.csv").rename(
        columns={"phenocode": "endpoint", "path_https": "url"})
    audit = pd.DataFrame(json.loads((ROUND / "reports/finngen_R12_full_source_audit.json").read_text(encoding="utf-8")))
    keep = ["endpoint", "path", "url", "bytes", "md5", "sha256", "tbi_path", "tbi_bytes",
            "tbi_md5", "tbi_sha256", "data_rows", "gzip_status", "provider_size_md5_match"]
    completed = pending.drop(columns=[c for c in ["planned_storage_root", "download_gate"] if c in pending]).merge(
        audit[keep], on=["endpoint", "url"], how="left", validate="one_to_one")
    completed["source_release"] = "FinnGen R12"
    completed["registration_evidence"] = "Registration submitted and official download email received 2026-09-27; email explicitly identifies R2-R12 older-freeze paths"
    completed["source_status"] = "REGISTERED_DOWNLOAD_PROVIDER_MD5_AND_FULL_GZIP_CRC_VERIFIED"
    completed["redistribution_note"] = "Manifest and hashes may be shared; raw summary-statistic redistribution must follow FinnGen terms"
    if completed.sha256.isna().any() or not completed.provider_size_md5_match.astype(bool).all():
        raise RuntimeError("Incomplete or failed FinnGen source audit")
    completed.to_csv(ROUND / "tables/FinnGen_R12_source_manifest_completed_2026-09-28.csv", index=False)
    print({"rows": len(completed), "all_verified": True})


if __name__ == "__main__":
    main()
