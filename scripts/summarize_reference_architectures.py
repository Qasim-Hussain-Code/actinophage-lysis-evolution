#!/usr/bin/env python3
"""Plot functional coverage and measured distances in assay-source controls."""
from __future__ import annotations
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import NA, read_tsv, sha256, write_json, write_tsv

SOURCE = ROOT / 'results/reference_architecture_controls'
OUT = SOURCE / 'comparison'

def run():
    manifest = json.loads((SOURCE / 'artifact_manifest.json').read_text())
    for name,digest in manifest['artifacts'].items():
        if sha256(SOURCE / name) != digest: raise ValueError('Source-survey artifact changed')
    evidence = read_tsv(SOURCE / 'tables/lysis_evidence.tsv')
    controls = {row['source_accession']:row for row in read_tsv(SOURCE / 'assayed_source_recovery.tsv')}
    rows = []
    for row in sorted(read_tsv(SOURCE / 'architecture_summary.tsv'), key=lambda r:r['phage'].lower()):
        accepted_a = [e for e in evidence if e['accession'] == row['source_accession'] and e['category'] == 'lysin_a'
                      and e['confidence'] in {'high','moderate'} and e['partial'] == 'false']
        if len(accepted_a) != int(row['lysin_a_count']): raise ValueError('A coverage differs from the architecture table')
        a_conf = ('high' if any(e['confidence']=='high' for e in accepted_a) else 'moderate') if accepted_a else 'unresolved'
        b_conf = controls[row['source_accession']]['confidence']
        rows.append({**row,'a_confidence':a_conf,'assayed_b_confidence':b_conf})
    write_tsv(OUT / 'source_coverage_and_distances.tsv', rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    import numpy as np
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    colors = ['#e6e8ea','#cc8a35','#326c96']
    encoded = {'unresolved':0,'moderate':1,'high':2}
    matrix = np.array([[encoded[r['a_confidence']], encoded[r['assayed_b_confidence']]] for r in rows])
    fig,axes = plt.subplots(1,2,figsize=(9.4,7.2),gridspec_kw={'width_ratios':[1,1.25]})
    axes[0].imshow(matrix,vmin=0,vmax=2,cmap=ListedColormap(colors),aspect='auto')
    axes[0].set(yticks=np.arange(len(rows)),yticklabels=[r['phage'] for r in rows],
                xticks=[0,1],xticklabels=['Lysin A','Assayed Lysin B'], title='A  Operational functional support')
    axes[0].tick_params(length=0)
    for y,r in enumerate(rows):
        for x,key in [(0,'a_confidence'),(1,'assayed_b_confidence')]:
            val=r[key]
            axes[0].text(x,y,{'unresolved':'?','moderate':'M','high':'H'}[val],ha='center',va='center',
                         color='white' if val=='high' else '#222222')
    for y,r in enumerate(rows):
        if r['intervening_cds'] == NA:
            axes[1].text(.05,y,'Unresolved A/B pair',va='center',color='#686868',fontsize=9)
        else:
            distance=int(r['intervening_cds'])
            axes[1].plot(distance,y,'o',color='#326c96',markersize=7)
            axes[1].text(distance+.15,y,str(distance),va='center')
    axes[1].axvline(4,color='#555555',ls='--',lw=1,label='Separated when >4')
    axes[1].set(xlim=(-.35,5.3),ylim=(len(rows)-.5,-.5),xticks=range(6),yticks=[],
                xlabel='Intervening CDS in supported A/B pair',title='B  CDS-order separation')
    axes[1].spines[['top','right','left']].set_visible(False)
    axes[1].legend(loc='lower right',frameon=False,fontsize=9)
    fig.legend(handles=[Patch(color=colors[2],label='H: high'),Patch(color=colors[1],label='M: moderate'),
                        Patch(color=colors[0],label='?: unresolved')],loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.075))
    fig.text(.025,.025,'Control-selected source genomes; B self-matches are expected. Unresolved is not absence. No population inference.',fontsize=9)
    fig.tight_layout(rect=[0,.14,1,1],w_pad=3)
    for suffix in ('png','pdf'):
        metadata = {'Software':'actinolysis'} if suffix=='png' else {'Creator':'actinolysis','CreationDate':None,'ModDate':None}
        fig.savefig(OUT / f'source_architectures.{suffix}',dpi=300,metadata=metadata)
    plt.close(fig)
    names=['source_coverage_and_distances.tsv','source_architectures.png','source_architectures.pdf']
    write_json(OUT / 'figure_manifest.json', {'caption':'Operational A/B support and CDS-order separation in thirteen complete sources of assayed B proteins. Five sources have a unique supported A/B pair; all five are closely linked. The other eight lack supported A assignments. The Ms6 partial locus is excluded. B source self-matches are expected positive controls, not independent validation or population sensitivity.',
        'source_manifest_sha256':sha256(SOURCE/'artifact_manifest.json'),'summary_script_sha256':sha256(__file__),
        'artifacts':{name:sha256(OUT/name) for name in names}})
    print(json.dumps({'sources':len(rows),'supported_pairs':sum(r['architecture']=='closely_linked' for r in rows),'status':'pass'}))

if __name__=='__main__':run()
