"""Descriptive neighborhood contrasts conditional on deposited lysis labels."""
from __future__ import annotations
from collections import Counter, defaultdict
from itertools import combinations
from . import architecture
from . import neighborhoods as nb

SCOPE = 'deposited_annotation_only'
PAIR_FIELDS = ['evidence_scope', 'pair_id', 'accession_a', 'accession_b', 'cluster', 'family_identity', 'window',
    'architecture_threshold', 'genome_subset', 'pair_architecture', 'metric', 'lysin_b_a', 'lysin_b_b', 'anchor_family',
    'target_score', 'control_mean', 'delta', 'n_controls', 'n_eligible_controls', 'status', 'reason']

def analyze(cfg, genomes, cds, evidence, memberships, annotation_architectures):
    genome_map = {g['accession']: g for g in genomes}
    if len(genome_map) != len(genomes) or len({g['protein_id'] for g in cds}) != len(cds):
        raise ValueError('Duplicate genome or CDS identifiers')
    labels, _, _ = architecture._assignments(cds, evidence, 'annotation')
    all_lysis = {row['protein_id'] for row in evidence if row['category'] in nb.LYSIS}
    all_lysis |= {pid for pid, categories in labels.items() if categories}
    per_genome = defaultdict(list)
    clusters = defaultdict(list)
    for row in cds:
        if row['accession'] not in genome_map: raise ValueError('CDS has no matching genome')
        per_genome[row['accession']].append(row)
    for row in genomes:
        if row['cluster'].lower() not in nb.MISSING: clusters[row['cluster']].append(row['accession'])
    annos = {(row['accession'], int(row['threshold'])): row for row in annotation_architectures if row['scope'] == 'annotation'}
    thresholds = sorted(set(int(value) for value in cfg['architecture']['thresholds']))
    minimum, maximum = max(5, int(cfg['statistics']['min_controls'])), int(cfg['statistics']['matched_controls'])
    if maximum < minimum: raise ValueError('Control maximum is below the unchanged minimum')
    pairs, controls, neighbors = [], [], []
    for identity in sorted(set(float(row['family_identity']) for row in memberships)):
        selected_members = [row for row in memberships if float(row['family_identity']) == identity]
        assignments = {row['protein_id']: row['family_id'] for row in selected_members}
        prevalence = {row['family_id']: float(row['prevalence']) for row in selected_members}
        genes = {acc: [{**row, 'family': assignments.get(row['protein_id'], 'UNCLUSTERED::' + row['protein_id'])} for row in rows]
                 for acc, rows in per_genome.items()}
        excluded = {assignments[pid] for pid in all_lysis if pid in assignments}
        for window in sorted(set(int(value) for value in cfg['neighborhood']['windows'])):
            windows = {row['protein_id']: nb.neighborhood(genes[acc], row['protein_id'], window,
                       str(genome_map[acc].get('topology', '')).lower() == 'circular')
                       for acc, rows in genes.items() for row in rows}
            for acc, rows in genes.items():
                for anchor in rows:
                    if 'lysin_b' not in labels.get(anchor['protein_id'], {}): continue
                    local = windows[anchor['protein_id']]
                    neighbors.append({'evidence_scope': SCOPE, 'family_identity': identity, 'window': window,
                        'accession': acc, 'anchor': anchor['protein_id'], 'anchor_family': anchor['family'],
                        'neighbor_proteins': ';'.join(row['protein_id'] for row in local),
                        'neighbor_families': ';'.join(row['family'] for row in local),
                        'signed_offsets': ';'.join(str(row['offset']) for row in local),
                        'relative_strands': ';'.join(str(row['relative_strand']) for row in local)})
            for cluster, accessions in sorted(clusters.items()):
                for acc_a, acc_b in combinations(sorted(accessions), 2):
                    pair_id = acc_a + '|' + acc_b
                    ga, gb = genome_map[acc_a], genome_map[acc_b]
                    aa = [g for g in genes[acc_a] if 'lysin_b' in labels.get(g['protein_id'], {})]
                    bb = [g for g in genes[acc_b] if 'lysin_b' in labels.get(g['protein_id'], {})]
                    reason, anchors = '', None
                    if len(aa) != 1 or len(bb) != 1:
                        reason = 'requires_one_deposited_b_label_per_genome'
                    else:
                        anchors = a, b = aa[0], bb[0]
                        if any(g['protein_id'] not in assignments for g in anchors): reason = 'anchor_excluded_from_family_clustering'
                        elif a['family'] != b['family']: reason = 'annotated_b_anchors_not_in_same_family'
                        elif any(sum(g['family'] == anchor['family'] for g in genes[acc]) != 1 for acc, anchor in [(acc_a, a), (acc_b, b)]):
                            reason = 'annotated_b_family_has_multiple_copies'
                        elif any(not nb.strand_value(g['strand']) for g in anchors): reason = 'anchor_strand_unresolved'
                        elif str(ga.get('similarity_group', '')).lower() not in nb.MISSING and ga.get('similarity_group') == gb.get('similarity_group'):
                            reason = 'same_near_identical_similarity_group'
                    chosen, available, target = [], 0, {metric: None for metric in nb.METRICS}
                    if not reason:
                        by_cut = [[annos.get((acc, cut)) for acc in (acc_a, acc_b)] for cut in thresholds]
                        if any(row is None for rows in by_cut for row in rows): raise ValueError('Missing deposited-annotation architecture')
                        valid_cuts = [rows for rows in by_cut if all(row['architecture'] in {'closely_linked', 'separated'} for row in rows)]
                        if not valid_cuts:
                            reason = 'annotation_architecture_complex_or_unresolved'
                        elif any(row['lysin_b'] != anchor['protein_id'] for rows in valid_cuts for row, anchor in zip(rows, anchors)):
                            reason = 'annotation_architecture_and_anchor_disagree'
                    if not reason:
                        a, b = anchors
                        target = nb.scores(windows[a['protein_id']], windows[b['protein_id']])
                        if any(target[metric] is None for metric in nb.METRICS): reason = 'empty_anchor_neighborhoods'
                        else:
                            chosen, available = nb.matched_controls(genes[acc_a], genes[acc_b], a, b, windows, prevalence,
                                                                   excluded, ga, gb, maximum)
                            controls.extend({'evidence_scope': SCOPE, 'pair_id': pair_id, 'family_identity': identity,
                                             'window': window, **row} for row in chosen)
                            if len(chosen) < minimum: reason = 'fewer_than_minimum_matched_controls'
                    for threshold in thresholds:
                        archs = [annos.get((acc, threshold), {}).get('architecture', 'unresolved') for acc in (acc_a, acc_b)]
                        pair_architecture = archs[0] if archs[0] == archs[1] else 'mixed'
                        for subset in ('all', 'representatives'):
                            failure = reason
                            if not failure and any(category not in {'closely_linked', 'separated'} for category in archs):
                                failure = 'annotation_architecture_complex_or_unresolved'
                            if subset == 'representatives' and not all('representative' in g for g in (ga, gb)):
                                failure = 'qc_representative_metadata_unavailable'
                            elif subset == 'representatives' and not all(nb._truth(g['representative']) for g in (ga, gb)):
                                failure = 'not_both_qc_representatives'
                            for metric in nb.METRICS:
                                values = [row[metric] for row in chosen if row[metric] is not None]
                                mean = sum(values) / len(values) if values else None
                                metric_failure = failure or ('fewer_than_minimum_scored_controls' if len(values) < minimum else '')
                                pairs.append({'evidence_scope': SCOPE, 'pair_id': pair_id, 'accession_a': acc_a, 'accession_b': acc_b,
                                    'cluster': cluster, 'family_identity': identity, 'window': window, 'architecture_threshold': threshold,
                                    'genome_subset': subset, 'pair_architecture': pair_architecture, 'metric': metric,
                                    'lysin_b_a': anchors[0]['protein_id'] if anchors else '', 'lysin_b_b': anchors[1]['protein_id'] if anchors else '',
                                    'anchor_family': anchors[0]['family'] if anchors else '', 'target_score': target[metric] if target[metric] is not None else '',
                                    'control_mean': mean if mean is not None else '',
                                    'delta': target[metric] - mean if not metric_failure else '', 'n_controls': len(values), 'n_eligible_controls': available,
                                    'status': 'descriptive_eligible' if not metric_failure else 'unavailable', 'reason': metric_failure})
    return pairs, controls, neighbors

def summarize(pairs):
    groups = defaultdict(list)
    fields = ['evidence_scope', 'family_identity', 'window', 'architecture_threshold', 'genome_subset', 'pair_architecture', 'metric']
    for row in pairs: groups[tuple(row[key] for key in fields)].append(row)
    summary = []
    for key, rows in sorted(groups.items()):
        clusters = defaultdict(list)
        for row in rows:
            if row['status'] == 'descriptive_eligible': clusters[row['cluster']].append(float(row['delta']))
        means = [sum(values) / len(values) for values in clusters.values()]
        summary.append({**dict(zip(fields, key)), 'n_candidate_pairs': len(rows),
            'n_available_pairs': sum(len(values) for values in clusters.values()), 'n_cluster_labels': len(clusters),
            'equal_cluster_mean_delta': sum(means) / len(means) if means else '',
            'interpretation': 'Descriptive conditional contrast; deposited labels do not independently establish lysis function',
            'formal_inference': 'not_performed'})
    return summary
