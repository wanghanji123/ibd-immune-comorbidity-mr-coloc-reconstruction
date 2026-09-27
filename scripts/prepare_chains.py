"""Fetch only public UCSC coordinate chains; no participant or outcome data."""
from pathlib import Path
import gzip
import hashlib
import json
import requests
from build_instruments import ROUND


def main():
    records = []
    session = requests.Session()
    for source, dest in [('hg19', 'Hg38'), ('hg38', 'Hg19')]:
        filename = source+'To'+dest+'.over.chain.gz'
        url = 'https://hgdownload.soe.ucsc.edu/goldenPath/'+source+'/liftOver/'+filename
        response = session.get(url, timeout=(20, 120))
        response.raise_for_status()
        assert len(response.content) < 10000000
        content = gzip.decompress(response.content)
        assert content.startswith(b'chain ')
        target = ROUND/'metadata'/filename
        if target.exists():
            assert target.read_bytes() == response.content, 'Remote chain changed; version explicitly'
        else:
            target.write_bytes(response.content)
        (ROUND/'metadata'/filename[:-3]).write_bytes(content)
        records.append(dict(url=url, compressed_bytes=len(response.content), sha256=hashlib.sha256(response.content).hexdigest(),
                            uncompressed_sha256=hashlib.sha256(content).hexdigest()))
    (ROUND/'metadata/coordinate_chain_sources.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    print(json.dumps(records), flush=True)


if __name__ == '__main__':
    main()
