"""Resume and verify the three missing FinnGen R12 outcomes and their tabix indexes."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import os
import time
import requests

DEST = Path('C:/Users/Lenovo/IBD_rebuild_runtime/finngen_R12')
BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/'
FILES = {
  'finngen_R12_H7_IRIDOCYCLITIS.gz': (813852600, '2c0a65138be5b6a5558e6466ce46c6bd'),
  'finngen_R12_H7_IRIDOCYCLITIS.gz.tbi': (1775602, '6693c4a7ba74f20af988ccb8c6eda5bc'),
  'finngen_R12_L12_PSORIASIS.gz': (814553655, 'f933068402790c9483218585a3967ba4'),
  'finngen_R12_L12_PSORIASIS.gz.tbi': (1775473, '1c75fda09aef12ef6fb768c2b9ca58b3'),
  'finngen_R12_M13_PSORIARTH.gz': (811909342, 'e43dd5dd071eca2b56a27c401edf7c3c'),
  'finngen_R12_M13_PSORIARTH.gz.tbi': (1775199, '276fd24a02c873a3929145d14bb2031f'),
}


def digest(path, algorithm='md5'):
    h = hashlib.new(algorithm)
    with path.open('rb') as inp:
        for block in iter(lambda: inp.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def download(name, expected_size, expected_md5):
    DEST.mkdir(parents=True, exist_ok=True)
    target = DEST/name
    partial = DEST/(name+'.part')
    if target.exists():
        if target.stat().st_size == expected_size and digest(target) == expected_md5:
            return dict(file=name, status='EXISTING_VERIFIED', bytes=expected_size, md5=expected_md5)
        raise RuntimeError(f'Existing completed file does not match provider metadata: {target}')
    for attempt in range(1, 11):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {'Range': f'bytes={offset}-'} if offset else {}
        try:
            with requests.get(BASE+name, headers=headers, stream=True, timeout=(20,120)) as response:
                response.raise_for_status()
                if offset and response.status_code != 206:
                    stale = partial.with_suffix(partial.suffix+f'.server_no_range_{int(time.time())}')
                    partial.replace(stale)
                    offset = 0
                mode = 'ab' if offset and response.status_code == 206 else 'wb'
                with partial.open(mode) as out:
                    for block in response.iter_content(8*1024*1024):
                        if block:
                            out.write(block)
                if partial.stat().st_size != expected_size:
                    raise IOError(f'incomplete {partial.stat().st_size}/{expected_size}')
                actual_md5 = digest(partial)
                if actual_md5 != expected_md5:
                    bad = partial.with_suffix(partial.suffix+f'.bad_md5_{int(time.time())}')
                    partial.replace(bad)
                    raise IOError(f'MD5 mismatch {actual_md5} != {expected_md5}; retained {bad}')
                os.replace(partial,target)
                return dict(file=name,status='DOWNLOADED_VERIFIED',bytes=expected_size,md5=actual_md5,
                            sha256=digest(target,'sha256'),url=BASE+name)
        except Exception as exc:
            if attempt == 10:
                raise
            print(json.dumps({'file':name,'attempt':attempt,'partial_bytes':partial.stat().st_size if partial.exists() else 0,
                              'retry_reason':str(exc)}),flush=True)
            time.sleep(min(30,2**attempt))


if __name__ == '__main__':
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(download,name,*meta):name for name,meta in FILES.items()}
        for future in as_completed(futures):
            result=future.result()
            records.append(result)
            print(json.dumps(result),flush=True)
    records.sort(key=lambda x:x['file'])
    (DEST/'download_verification.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    assert len(records)==len(FILES) and all(x['status'].endswith('VERIFIED') for x in records)
    print('PASS: all missing R12 outcomes and indexes match provider MD5',flush=True)
