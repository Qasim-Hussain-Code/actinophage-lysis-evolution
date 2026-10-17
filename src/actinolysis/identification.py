"""Independent curated-reference homology, separated from annotation labels."""
from __future__ import annotations
from collections import Counter
import re
import subprocess
import urllib.parse
from Bio import SeqIO
from .common import ROOT, NA, budget_check, read_tsv, sha256, utcnow, write_json, write_tsv
from .acquisition import download

def annotation_category(product):
    p=product.lower().replace('_',' ')
    if re.search(r'\b(lysin[ -]?b|lysb)\b',p) or ('mycol' in p and ('esterase' in p or 'hydrolase' in p)):
        return 'lysin_b'
    if re.search(r'\b(lysin[ -]?a|lysa|endolysin)\b',p): return 'lysin_a'
    if 'holin' in p: return 'holin_like'  # A deposited label is not an experimentally established holin.
    return 'unresolved'

def build_references(cfg):
    specifications=read_tsv(ROOT/'references/protein_controls.tsv')
    frozen=ROOT/'data/manifests/reference_manifest.tsv'
    prior=read_tsv(frozen) if frozen.exists() else []
    prior_by_id={r['reference_id']:r for r in prior}
    pinned_sources={r['source_accession']:r['resolved_source_accession'] for r in prior}
    rows=[]; sequences=[]; source_records={}
    for spec in specifications:
        source=spec['source_accession']
        if source not in source_records:
            requested=pinned_sources.get(source,source)
            url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'nuccore','id':requested,'rettype':'gbwithparts','retmode':'text','tool':'actinolysis'})
            path=download(url,ROOT/f'data/raw/references/{source}.gb',cfg,cfg.get('offline',False),maximum=3_000_000)
            record=SeqIO.read(path,'genbank')
            if source in pinned_sources and record.id != requested: raise ValueError(f'Pinned reference source version mismatch: {record.id} vs {requested}')
            source_records[source]=(record,path,url)
        record,path,url=source_records[source]
        # Retain the original resolution URL in an existing freeze. Replay uses
        # its resolved accession version above, never the current base record.
        provenance_url=prior_by_id.get(spec['reference_id'],{}).get('source_url')
        if not provenance_url:
            provenance_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'nuccore','id':record.id,'rettype':'gbwithparts','retmode':'text','tool':'actinolysis'})
        target=spec.get('protein_accession','')
        matches=[]
        for f in record.features:
            if f.type!='CDS': continue
            accession=f.qualifiers.get('protein_id',[''])[0]
            gene=f.qualifiers.get('gene',[''])[0]
            product=f.qualifiers.get('product',[''])[0].lower()
            gene_match=bool(spec['gene']) and (gene.lower() in {spec['gene'].lower(),'gp'+spec['gene'].lower()}
                       or re.fullmatch(r'(?:gp|gene|protein)[ _-]?'+re.escape(spec['gene']),product) is not None)
            capsid_match=spec['reference_id']=='D29_capsid_negative' and 'major capsid' in product
            accession_match = accession == target if '.' in target else accession.split('.')[0] == target
            if (target and accession_match) or (not target and (gene_match or capsid_match)):
                matches.append(f)
        if len(matches)!=1: raise ValueError(f"Reference {spec['reference_id']} resolves to {len(matches)} CDS in {record.id}")
        f=matches[0]; protein=f.qualifiers.get('translation',[''])[0].replace(' ','').replace('\n','')
        if not protein or '*' in protein: raise ValueError(f"Invalid characterized reference sequence: {spec['reference_id']}")
        rows.append({**spec,'resolved_source_accession':record.id,'resolved_protein_accession':f.qualifiers['protein_id'][0],
            'length_aa':len(protein),'start':int(f.location.start),'end':int(f.location.end),
            'source_url':provenance_url,'source_sha256':sha256(path),'protein_sha256':__import__('hashlib').sha256(protein.encode()).hexdigest()})
        sequences.append(f">{spec['reference_id']}\n{protein}\n")
    path=ROOT/'data/processed/references.faa'; path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(''.join(sequences),encoding='ascii')
    if frozen.exists():
        old=read_tsv(frozen)
        normalized=[{k:str(v) for k,v in r.items()} for r in rows]
        if old!=normalized:
            raise ValueError('Reference sequence or curated scientific provenance changed against frozen reference manifest')
    else: write_tsv(frozen,rows)
    return rows

def run(cfg):
    build_references(cfg)
    output=ROOT/'data/processed/reference_hits.tsv'; output.parent.mkdir(parents=True,exist_ok=True)
    tmp=ROOT/'.tmp/reference_search'; tmp.mkdir(parents=True,exist_ok=True)
    exe=str((ROOT/cfg['tools']['mmseqs']).resolve()) if '/' in cfg['tools']['mmseqs'] else cfg['tools']['mmseqs']
    cmd=[exe,'easy-search',str(ROOT/'data/processed/proteins.faa'),str(ROOT/'data/processed/references.faa'),str(output),str(tmp),
        '--threads',str(min(int(cfg.get('threads',4)),4)),'-s','7.5','-k','5','--alignment-mode','3','-e','0.001','--max-seqs','100','--min-seq-id','0.25',
        '-c','0.5','--cov-mode','0','--split-memory-limit',cfg['resources']['mmseqs_memory'],
        '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits','--remove-tmp-files','1']
    budget_check(cfg,600_000_000)
    log=ROOT/'results/logs/reference_search.log'; log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w') as stream: subprocess.run(cmd,check=True,stdout=stream,stderr=subprocess.STDOUT,cwd=ROOT)
    lines=output.read_text(encoding='utf-8').splitlines()
    output.write_text(''.join(line+'\n' for line in sorted(lines)),encoding='utf-8')
    write_json(ROOT/'results/reports/reference_search.json',{'command':cmd,'hit_file_sha256':sha256(output),'reference_count':len(read_tsv(ROOT/'data/manifests/reference_manifest.tsv'))})

def evaluate_hit(hit,ref,cfg):
    identity=float(hit['identity']); qcov=float(hit['query_coverage']); tcov=float(hit['target_coverage']); evalue=float(hit['evalue'])
    if ref['category']=='negative_control': return 'unresolved'
    if ref['evidence_type']!='experimental': return 'low'
    limits=cfg['identification']
    if identity>=0.999999 and qcov>=0.999999 and tcov>=0.999999 and evalue<=limits['max_evalue']: return 'high'
    if identity>=limits['min_identity'] and qcov>=limits['min_query_coverage'] and tcov>=limits['min_target_coverage'] and evalue<=limits['max_evalue']:
        return 'moderate'
    return 'low'

def evidence(cfg):
    refs={r['reference_id']:r for r in read_tsv(ROOT/'data/manifests/reference_manifest.tsv')}
    cds={r['protein_id']:r for r in read_tsv(ROOT/'data/processed/cds.tsv')}
    fields=['protein_id','reference_id','identity','alignment_length','query_coverage','target_coverage','evalue','bits']
    import csv
    with (ROOT/'data/processed/reference_hits.tsv').open() as stream: hits=list(csv.DictReader(stream,fieldnames=fields,delimiter='\t'))
    scores={}; negative=[]
    for hit in hits:
        if hit['protein_id'] not in cds or hit['reference_id'] not in refs: raise ValueError('Unexpected MMseqs hit identifiers')
        ref=refs[hit['reference_id']]
        if ref['category']=='negative_control': negative.append(hit); continue
        key=(hit['protein_id'],ref['category'])
        confidence=evaluate_hit(hit,ref,cfg)
        rank=(-{'high':3,'moderate':2,'low':1,'unresolved':0}[confidence],float(hit['evalue']),-float(hit['bits']),hit['reference_id'])
        if key not in scores or rank<scores[key][0]: scores[key]=(rank,hit)
    result=[]
    for p in cds.values():
        annot=annotation_category(p['product']+' '+p.get('note',''))
        existing=[]
        for (pid,category),(_,hit) in scores.items():
            if pid!=p['protein_id']: continue
            ref=refs[hit['reference_id']]; confidence=evaluate_hit(hit,ref,cfg)
            if p['partial']=='true': confidence='low'
            existing.append(category)
            result.append({**{k:p[k] for k in ['protein_id','protein_accession','accession','start','end','strand','product','partial']},
                'annotation_note':p.get('note',''),
                'category':category,'confidence':confidence,'annotation_category':annot,'source':'curated_reference_homology',
                **{k:hit[k] for k in ['reference_id','identity','query_coverage','target_coverage','evalue']},
                'reference_evidence_type':ref['evidence_type'],'evidence_url':ref['reference_url'],'context_used':'false'})
        if annot!='unresolved' and annot not in existing:
            result.append({**{k:p[k] for k in ['protein_id','protein_accession','accession','start','end','strand','product','partial']},
                'annotation_note':p.get('note',''),
                'category':annot,'confidence':'low','annotation_category':annot,'source':'genbank_annotation',
                'reference_id':NA,'identity':NA,'query_coverage':NA,'target_coverage':NA,'evalue':NA,
                'reference_evidence_type':'deposited_annotation','evidence_url':f"https://www.ncbi.nlm.nih.gov/nuccore/{p['accession']}",'context_used':'false'})
    schema=['protein_id','protein_accession','accession','start','end','strand','product','partial','annotation_note','category','confidence','annotation_category','source','reference_id','identity','query_coverage','target_coverage','evalue','reference_evidence_type','evidence_url','context_used']
    write_tsv(ROOT/'data/processed/lysis_evidence.tsv',result,schema)
    write_tsv(ROOT/'results/tables/lysis_evidence.tsv',result,schema)
    coverage=[]
    for r in read_tsv(ROOT/'data/processed/genomes.tsv'):
        local=[p for p in cds.values() if p['accession']==r['accession']]
        assignments=[p for p in result if p['accession']==r['accession']]
        coverage.append({'accession':r['accession'],'phage':r['phage'],'host_genus':r['host_genus'],'cluster':r['cluster'],'cds':len(local),
            'with_specific_annotation':sum(annotation_category(p['product']+' '+p.get('note',''))!='unresolved' for p in local),
            'supported_lysis_proteins':len({p['protein_id'] for p in assignments if p['confidence'] in {'high','moderate'}}),
            'annotation_only_lysis_proteins':len({p['protein_id'] for p in assignments if p['confidence']=='low' and p['source']=='genbank_annotation'}),
            'hypothetical_proteins':sum('hypothetical' in p['product'].lower() for p in local)})
    write_tsv(ROOT/'results/tables/annotation_coverage.tsv',coverage)
    write_tsv(ROOT/'results/tables/negative_control_hits.tsv',negative,fields)
    write_json(ROOT/'results/reports/functional_identification.json',{'assignments':len(result),'by_confidence':dict(Counter(p['confidence'] for p in result)),
        'limitations':['Only a small experimentally characterized reference panel; divergent functions can be undetected.',
        'No transmembrane predictor used; membrane features alone cannot establish holin function.',
        'Deposited annotation alone is low confidence; neighborhood never upgrades function.']})
