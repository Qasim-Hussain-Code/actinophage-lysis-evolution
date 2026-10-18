"""GenBank validation and explicit CDS extraction without reannotation."""
from __future__ import annotations
import hashlib
import json
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition
from Bio.SeqUtils import gc_fraction
from .common import ROOT, NA, read_tsv, sha256, write_json, write_tsv

def canonical_hash(sequence):
    sequence = str(sequence).upper()
    reverse = str(Seq(sequence).reverse_complement())
    return hashlib.sha256(min(sequence, reverse).encode()).hexdigest()

def sketch(sequence, k=21, size=1000):
    sequence = str(sequence).upper(); hashes = set()
    # bounded max ~200k integers per phage, discard after retaining bottom1000.
    for i in range(len(sequence)-k+1):
        word = sequence[i:i+k]
        if set(word) <= set('ACGT'):
            word = min(word,str(Seq(word).reverse_complement()))
            hashes.add(int.from_bytes(hashlib.blake2b(word.encode(),digest_size=8).digest(),'little'))
    return set(sorted(hashes)[:size])

def translation_check(feature, record):
    given = ''.join(feature.qualifiers.get('translation',[''])[0].split()).removesuffix('*').upper()
    if not given: return '', 'missing_translation'
    if feature.qualifiers.get('pseudo') is not None or 'pseudogene' in feature.qualifiers:
        return given,'pseudogene'
    sequence = feature.extract(record.seq)
    offset = int(feature.qualifiers.get('codon_start',['1'])[0])-1
    sequence = sequence[offset:]
    if len(sequence)%3: return given,'incomplete_codon'
    table = int(feature.qualifiers.get('transl_table',['11'])[0])
    computed = str(sequence.translate(table=table)).removesuffix('*').upper()
    if '*' in given: return given,'translation_mismatch'
    if 'transl_except' in feature.qualifiers and computed != given:
        return given,'translation_exception_unresolved'
    if '*' in computed: return given,'translation_mismatch'
    alternative_start = str(sequence[:3]).upper() in CodonTable.unambiguous_dna_by_id[table].start_codons
    if computed == given or (alternative_start and len(computed)==len(given) and given.startswith('M') and computed[1:]==given[1:]):
        return given,'verified'
    if 'transl_except' in feature.qualifiers: return given,'translation_exception_unresolved'
    return given,'translation_mismatch'

def extract_record(record, metadata):
    length = len(record.seq)
    if length < 10000 or length > 500000: raise ValueError('genome_length_out_of_bounds')
    if not set(str(record.seq).upper()) <= set('ACGTRYSWKMBDHVN'): raise ValueError('invalid_nucleotide_symbols')
    if sum(str(record.seq).upper().count(s) for s in 'RYSWKMBDHVN') / length > 0.01:
        raise ValueError('more_than_one_percent_ambiguous_bases')
    # RefSeq titles may use "virus" while ORGANISM or the PHG division
    # explicitly identifies the record as a bacteriophage.
    phage_names = (record.description, record.annotations.get('organism', ''))
    if not (any('phage' in name.lower() for name in phage_names)
            or record.annotations.get('data_file_division') == 'PHG'):
        raise ValueError('not_phage_record')
    if not any(w in record.description.lower() for w in ['complete genome','complete sequence']):
        raise ValueError('not_explicitly_complete')
    features = sorted([f for f in record.features if f.type=='CDS'],key=lambda f:(int(f.location.start),int(f.location.end),f.location.strand or 0))
    if not features: raise ValueError('no_cds_annotations')
    cds = []
    for index,feature in enumerate(features):
        parts = [[int(p.start),int(p.end)] for p in feature.location.parts]
        if any(not (0<=start<end<=length) for start,end in parts): raise ValueError('invalid_cds_coordinates')
        partial = any(not isinstance(pos,ExactPosition) for part in feature.location.parts for pos in [part.start,part.end])
        protein,status = translation_check(feature,record)
        if status=='translation_mismatch': raise ValueError(f'translation_mismatch_cds_{index}')
        if '*' in protein or not set(protein)<=set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
            raise ValueError(f'invalid_protein_cds_{index}')
        if status!='verified': partial=True
        q=feature.qualifiers
        cds.append({'protein_id':f'{record.id}__cds{index:04d}','protein_accession':q.get('protein_id',[NA])[0],
            'accession':record.id,'cds_index':index,'start':int(feature.location.start),'end':int(feature.location.end),
            'location_parts':json.dumps(parts,separators=(',',':')),'strand':feature.location.strand or 0,'partial':str(partial).lower(),
            'product':q.get('product',[NA])[0].replace('\t',' ').replace('\n',' '),
            'note':' | '.join(q.get('note',[])).replace('\t',' ').replace('\n',' '),
            'translation':protein,'gene':q.get('gene',[NA])[0],'locus_tag':q.get('locus_tag',[NA])[0],
            'translation_status':status})
    source = next((f for f in record.features if f.type=='source'),None)
    source_host = NA if source is None else source.qualifiers.get('host',source.qualifiers.get('lab_host',[NA]))[0]
    result={**metadata,'length':length,'gc_percent':round(100*gc_fraction(record.seq),6),'cds_count':len(cds),
        'topology':record.annotations.get('topology',NA),'annotation_source':'NCBI_GenBank',
        'source_host':source_host,'sequence_sha256':canonical_hash(record.seq),'inclusion_status':'included','exclusion_reason':'not_applicable'}
    return result,cds

def run(cfg):
    manifest=read_tsv(ROOT/f"data/manifests/{cfg['mode']}_inputs.tsv")
    genomes=[]; proteins=[]; exclusions=[]; sketches={}; seen={}
    for row in manifest:
        try:
            if sha256(ROOT/row['path']) != row['sha256']: raise RuntimeError(f"Raw checksum mismatch: {row['accession']}")
            record=SeqIO.read(ROOT/row['path'],'genbank')
            if record.id!=row['accession']: raise ValueError('accession_version_mismatch')
            genome,cds=extract_record(record,row)
            digest=genome['sequence_sha256']
            if digest in seen: raise ValueError(f'exact_sequence_duplicate_of:{seen[digest]}')
            seen[digest]=record.id; genomes.append(genome); proteins.extend(cds)
            sketches[record.id]=sketch(record.seq)
        except ValueError as error:
            exclusions.append({**row,'inclusion_status':'excluded','exclusion_reason':str(error)})
    if not genomes: raise ValueError('No genomes pass QC; see exclusions')
    # Deterministic connected components are screening groups, not inferred clades or ANI.
    parent={g['accession']:g['accession'] for g in genomes}
    def find(a):
        while parent[a]!=a: a=parent[a]
        return a
    near=[]
    accessions=sorted(parent)
    for i,a in enumerate(accessions):
        for b in accessions[i+1:]:
            j=len(sketches[a]&sketches[b])/max(1,len(sketches[a]|sketches[b]))
            if j>=0.95:
                x,y=sorted([find(a),find(b)]); parent[y]=x
                near.append({'accession_a':a,'accession_b':b,'sketch_jaccard':j,'threshold':0.95,'k':21,'sketch_size':1000})
    for g in genomes:
        g['similarity_group']=find(g['accession']); g['representative']=str(g['accession']==find(g['accession'])).lower()
    write_tsv(ROOT/'data/processed/genomes.tsv',genomes)
    write_tsv(ROOT/'data/processed/cds.tsv',proteins)
    write_tsv(ROOT/'results/tables/dataset_manifest.tsv',genomes)
    fields=list(manifest[0]); write_tsv(ROOT/'results/tables/exclusions.tsv',exclusions,fields)
    write_tsv(ROOT/'results/tables/near_identical_pairs.tsv',near,['accession_a','accession_b','sketch_jaccard','threshold','k','sketch_size'])
    fasta=ROOT/'data/processed/proteins.faa'; fasta.parent.mkdir(parents=True,exist_ok=True)
    fasta.write_text(''.join(f">{p['protein_id']}\n{p['translation']}\n" for p in proteins if p['translation'] and p['partial']=='false'),encoding='ascii')
    write_json(ROOT/'results/reports/quality_control.json',{'selected':len(manifest),'included':len(genomes),'excluded':len(exclusions),
        'cds':len(proteins),'translation_status':{s:sum(p['translation_status']==s for p in proteins) for s in sorted({p['translation_status'] for p in proteins})},
        'near_identical_screen':'Bottom1000 canonical 21mer sketch Jaccard >=0.95; proxy, not ANI; no inferred phylogenetic independence'})
