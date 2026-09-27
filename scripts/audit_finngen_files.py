"""Full compressed-stream audit of all frozen FinnGen R12 outcome inputs."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import gzip
import hashlib
import json

ROUND = Path(__file__).resolve().parents[1]
PROJECT = ROUND.parents[1]
RUNTIME = Path('C:/Users/Lenovo/IBD_rebuild_runtime/finngen_R12')
BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/'
EXPECTED_HEADER = ['#chrom','pos','ref','alt','rsids','nearest_genes','pval','mlogp','beta','sebeta','af_alt','af_alt_cases','af_alt_controls']
SOURCES = {
 'M13_ANKYLOSPON': (PROJECT/'data_external/finngen_as_r12/finngen_R12_M13_ANKYLOSPON.gz',811239151,'eed0bdc572ec59c55945f99ab43d5d26',1774653,'daa901317e768583c944826277e3f76f'),
 'M13_ANKYLOSPON_STRICT': (PROJECT/'data_external/finngen_as_r12_strict/finngen_R12_M13_ANKYLOSPON_STRICT.gz',804365077,'9c51bae1f7aa6179014acb6a99e13f60',1773171,'6cacd9c5e6d108615a329f281ea3671a'),
 'H7_IRIDOCYCLITIS': (RUNTIME/'finngen_R12_H7_IRIDOCYCLITIS.gz',813852600,'2c0a65138be5b6a5558e6466ce46c6bd',1775602,'6693c4a7ba74f20af988ccb8c6eda5bc'),
 'L12_PSORIASIS': (RUNTIME/'finngen_R12_L12_PSORIASIS.gz',814553655,'f933068402790c9483218585a3967ba4',1775473,'1c75fda09aef12ef6fb768c2b9ca58b3'),
 'M13_PSORIARTH': (RUNTIME/'finngen_R12_M13_PSORIARTH.gz',811909342,'e43dd5dd071eca2b56a27c401edf7c3c',1775199,'276fd24a02c873a3929145d14bb2031f'),
}


def hashes(path):
    md5,sha=hashlib.md5(),hashlib.sha256()
    with path.open('rb') as inp:
        for block in iter(lambda:inp.read(8*1024*1024),b''):
            md5.update(block); sha.update(block)
    return md5.hexdigest(),sha.hexdigest()


def audit(item):
    endpoint,(path,size,expected_md5,tbi_size,tbi_md5)=item
    tbi=Path(str(path)+'.tbi')
    assert path.exists() and tbi.exists(), (path,tbi)
    assert path.stat().st_size==size and tbi.stat().st_size==tbi_size
    actual_md5,actual_sha=hashes(path)
    actual_tbi_md5,actual_tbi_sha=hashes(tbi)
    assert actual_md5==expected_md5 and actual_tbi_md5==tbi_md5
    lines=0; uncompressed=0; first=b''; final=b''
    with gzip.open(path,'rb') as inp:
        while True:
            block=inp.read(8*1024*1024)
            if not block: break
            if not first: first=block.splitlines()[0]
            final=(final+block)[-4096:]
            lines+=block.count(b'\n'); uncompressed+=len(block)
    header=first.decode('ascii').split('\t')
    assert header==EXPECTED_HEADER
    assert lines>1 and final.endswith(b'\n')
    return dict(endpoint=endpoint,path=str(path),url=BASE+path.name,bytes=size,md5=actual_md5,sha256=actual_sha,
      tbi_path=str(tbi),tbi_bytes=tbi_size,tbi_md5=actual_tbi_md5,tbi_sha256=actual_tbi_sha,
      header=header,data_rows=lines-1,uncompressed_bytes=uncompressed,
      gzip_status='FULL_EOF_CRC_PASS',provider_size_md5_match=True)


if __name__=='__main__':
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        tasks={pool.submit(audit,item):item[0] for item in SOURCES.items()}
        for task in as_completed(tasks):
            result=task.result(); records.append(result)
            print(json.dumps({k:result[k] for k in ['endpoint','data_rows','gzip_status']}),flush=True)
    records.sort(key=lambda x:x['endpoint'])
    (ROUND/'reports/finngen_R12_full_source_audit.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    assert len(records)==5
    print('PASS: five outcome files and indexes match provider metadata; all gzip streams reach CRC-checked EOF',flush=True)
