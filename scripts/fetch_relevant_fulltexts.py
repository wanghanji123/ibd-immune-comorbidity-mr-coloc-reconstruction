"""Retrieve primary full-text JATS for critical prior-study comparison, not training."""
import hashlib
import json
import requests
import xml.etree.ElementTree as ET
from build_instruments import ROUND

STUDIES = {
    'PMC9475439':'10.1001/jamadermatol.2022.3682',
    'PMC11881308':'10.1186/s12876-024-03566-4',
    'PMC11193636':'10.1053/j.gastro.2024.02.026',
    'PMC10290755':None,
}


def main():
    folder=ROUND/'metadata/literature'
    folder.mkdir(exist_ok=True)
    records=[]
    for pmc,expected in STUDIES.items():
        url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'
        try:
            response=requests.get(url,timeout=(15,60))
            response.raise_for_status()
            xml=ET.fromstring(response.content)
            meta=xml.find('./front/article-meta')
            assert meta is not None
            dois=[x.text for x in meta.findall("article-id[@pub-id-type='doi']")]
            assert expected is None or expected in dois,(pmc,dois)
            title=''.join(meta.find('title-group/article-title').itertext())
            (folder/f'{pmc}.xml').write_bytes(response.content)
            sections=[]
            for node in xml.findall('./body/sec'):
                sections.append('\n'.join(''.join(x.itertext()) for x in node.iter() if x.tag in {'title','p','table-wrap'}))
            (folder/f'{pmc}_body.txt').write_text('\n\n'.join(sections),encoding='utf-8')
            licenses=[''.join(x.itertext()) for x in meta.findall('./permissions/license')]
            records.append(dict(pmc=pmc,url=url,doi=dois,title=title,status='FULLTEXT_RETRIEVED_NOT_AUTOMATICALLY_REVIEWED',
                                licenses=licenses,sha256=hashlib.sha256(response.content).hexdigest()))
        except Exception as exc:
            records.append(dict(pmc=pmc,url=url,status='NOT_RETRIEVED',error=str(exc)))
        print(pmc,records[-1]['status'],flush=True)
    (ROUND/'metadata/literature_sources.json').write_text(json.dumps(records,indent=2),encoding='utf-8')


if __name__ == '__main__':
    main()
