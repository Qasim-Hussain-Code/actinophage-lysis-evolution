#!/usr/bin/env python3
"""Run only the five prespecified stages of a prepared exploratory screen."""
from __future__ import annotations

import argparse
import copy
import importlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import budget_check, load_config, read_tsv, sha256, utcnow, write_json, write_tsv

STAGES = {1: ('acquisition', 'run'), 2: ('quality', 'run'),
          3: ('identification', 'run'), 4: ('identification', 'evidence'),
          5: ('architecture', 'run')}
OUTPUTS = {
    2: ['data/processed/genomes.tsv', 'data/processed/cds.tsv', 'data/processed/proteins.faa',
        'results/tables/dataset_manifest.tsv', 'results/tables/exclusions.tsv',
        'results/tables/near_identical_pairs.tsv', 'results/reports/quality_control.json'],
    3: ['data/processed/reference_hits.tsv', 'data/manifests/reference_manifest.tsv',
        'data/processed/references.faa', 'results/reports/reference_search.json'],
    4: ['data/processed/lysis_evidence.tsv', 'results/tables/lysis_evidence.tsv',
        'results/tables/annotation_coverage.tsv', 'results/tables/negative_control_hits.tsv',
        'results/reports/functional_identification.json'],
    5: ['data/processed/architectures.tsv', 'data/processed/architecture_pairs.tsv',
        'results/tables/architectures.tsv', 'results/tables/architecture_pairs.tsv',
        'results/tables/architecture_sensitivity.tsv'],
}


def digest_object(value):
    import hashlib
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def inside(root, relative):
    path = (Path(root) / relative).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError('Input path escapes its declared root')
    return path


def tree_hashes(root, directories=(), files=()):
    paths = [inside(root, name) for name in files]
    paths += [p for directory in directories for p in inside(root, directory).rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc']
    return {p.relative_to(root).as_posix(): sha256(p) for p in sorted(set(paths))}


def verify_preparation(root, workspace):
    experiments = inside(root, '.experiments')
    prepared = inside(experiments, workspace)
    if prepared == experiments:
        raise ValueError('A named preparation workspace is required')
    path = prepared / 'preparation.json'
    report = json.loads(path.read_text(encoding='utf-8'))
    if (report.get('status') != 'awaiting_sequence_acquisition'
            or report.get('planned_selected_records') != 50
            or report.get('planned_frozen_pilot_records') != 20):
        raise ValueError('Expected a frozen 50-record exploratory preparation')
    fixed = {report['configuration_file']: report['configuration_file_sha256'],
             report['eligibility_snapshot']: report['eligibility_snapshot_sha256'],
             'data/manifests/phagesdb_snapshot.json': report['phagesdb_snapshot_sha256'],
             'data/manifests/pilot_inputs.tsv': report['pilot_manifest_sha256'],
             'docs/exploratory_screen_protocol.md': report['protocol_sha256']}
    fixed.update(report['code_sha256'])
    for name, expected in fixed.items():
        if sha256(inside(root, name)) != expected:
            raise ValueError(f'Prepared source changed: {name}')
    if tree_hashes(root, ('src',)) != report['source_tree_sha256']:
        raise ValueError('Prepared source tree changed')
    if report.get('runner_sha256', sha256(__file__)) != sha256(__file__):
        raise ValueError('Prepared exploratory runner changed')
    panel = inside(root, report['panel_root'])
    if sha256(panel / 'artifact_manifest.json') != report['panel_artifact_manifest_sha256']:
        raise ValueError('Frozen panel artifact manifest changed')
    artifacts = json.loads((panel / 'artifact_manifest.json').read_text(encoding='utf-8'))['artifacts']
    for name, expected in report['panel_files_sha256'].items():
        if artifacts.get(name) != expected or sha256(inside(panel, name)) != expected:
            raise ValueError('Frozen reference panel changed')
    cfg = report['effective_configuration']
    configured = load_config(inside(root, report['configuration_file']))
    configured['mode'] = 'intermediate'
    if (digest_object(cfg) != report['effective_configuration_sha256'] or cfg != configured
            or cfg['mode'] != 'intermediate' or cfg['sample_sizes']['intermediate'] != 50):
        raise ValueError('Prepared 50-record configuration changed')
    ledger_path = prepared / 'selection.tsv'
    if sha256(ledger_path) != report['selection_sha256']:
        raise ValueError('Prepared selection ledger changed')
    ledger = read_tsv(ledger_path)
    if (len(ledger) != 50 or [int(r['rank']) for r in ledger] != list(range(1, 51))
            or len({r['accession'].split('.')[0] for r in ledger}) != 50
            or [r['pilot'] for r in ledger] != ['true'] * 20 + ['false'] * 30):
        raise ValueError('Invalid prepared selection order or pilot nesting')
    return prepared, report, ledger


def pilot_hashes(root):
    files = ['data/manifests/pilot_inputs.tsv']
    files += [r['path'] for r in read_tsv(root / files[0])]
    return tree_hashes(root, ('data/processed', 'results/tables', 'results/reports', 'figures'), files)


def reference_cache(root, row):
    relative = Path('data/raw/references') / (row['source_accession'] + '.gb')
    prefixes = ('', '.experiments/lysterase_2026', '.experiments/lysterase_2026_replay',
                '.experiments/amendment_profile_v1', '.experiments/amendment_experimental_v1')
    candidates = [root / prefix / relative for prefix in prefixes]
    candidates.append(root / 'data/raw/lysterase_update' / (row['resolved_source_accession'] + '.gb'))
    for path in candidates:
        if path.is_file():
            if sha256(path) != row['source_sha256']:
                raise ValueError('Reference cache differs from the frozen 39-reference manifest')
            return path
    raise FileNotFoundError('A checksum-verified reference source is unavailable')


def root_guard(cfg, root, proposed=0):
    inventory = root / 'results/reports/software_storage.json'
    if inventory.exists():
        frozen = json.loads(inventory.read_text(encoding='utf-8'))['immutable_directories']
        if any(inside(root, '.experiments').is_relative_to(inside(root, name)) for name in frozen):
            raise ValueError('Immutable storage inventory cannot exclude the exploratory workspace')
    return budget_check(cfg, proposed, root=root)


class AcquisitionLedger:
    """Freeze every successful source before unchanged acquisition can proceed."""
    def __init__(self, module, runtime, selection):
        self.module, self.runtime, self.selection = module, runtime, selection
        self.path = runtime / 'data/manifests/acquisition_ledger.tsv'
        self.download = module.download
        self.rows = read_tsv(self.path) if self.path.exists() else [
            {'rank': str(i), 'phage': r['phage'], 'requested_accession': r['accession'],
             'pilot': r['pilot'], 'status': 'not_attempted', 'resolved_accession': '',
             'source_sha256': '', 'cache_path': '', 'error_type': '', 'failure_reason': ''}
            for i, r in enumerate(selection, 1)]
        expected = [(str(i), r['phage'], r['accession'], r['pilot']) for i, r in enumerate(selection, 1)]
        if [(r['rank'], r['phage'], r['requested_accession'], r['pilot']) for r in self.rows] != expected:
            raise ValueError('Private acquisition ledger differs from the frozen selection')
        self.by_base = {r['requested_accession'].split('.')[0]: r for r in self.rows}
        self.save()
        self.validate_cache()

    def save(self):
        temporary = self.path.with_suffix('.tsv.next')
        write_tsv(temporary, self.rows)
        temporary.replace(self.path)

    def fail(self, row, reason, error_type='SourceIntegrityError'):
        row.update(status='failed', error_type=error_type, failure_reason=reason)
        self.save()

    def validate_cache(self):
        allowed = {}
        for selected, row in zip(self.selection, self.rows):
            if row['status'] not in {'not_attempted', 'failed', 'acquired'}:
                raise ValueError('Invalid private acquisition status')
            if selected['pilot'] == 'true':
                for accession in (selected['accession'], selected['frozen_pilot_accession']):
                    allowed['data/raw/genbank/' + accession + '.gb'] = (row, selected['frozen_pilot_sha256'])
            if row['resolved_accession']:
                accession = row['resolved_accession']
                expected = 'data/raw/genbank/' + accession + '.gb'
                if (not re.fullmatch(r'[A-Z]{1,6}_?\d{5,10}\.\d+', accession)
                        or accession.split('.')[0] != selected['accession'].split('.')[0]
                        or ('.' in selected['accession'] and accession != selected['accession'])
                        or row['cache_path'] != expected
                        or not re.fullmatch(r'[0-9a-f]{64}', row['source_sha256'])
                        or (selected['pilot'] == 'true' and (accession, row['source_sha256']) !=
                            (selected['frozen_pilot_accession'], selected['frozen_pilot_sha256']))):
                    self.fail(row, 'invalid_frozen_source_identity')
                    raise ValueError('Invalid frozen acquisition source identity')
                allowed[expected] = (row, row['source_sha256'])
                if not inside(self.runtime, expected).is_file():
                    self.fail(row, 'registered_source_cache_missing')
                    raise ValueError('Registered acquisition source cache is missing')
            elif row['status'] == 'acquired':
                self.fail(row, 'acquired_source_has_no_frozen_identity')
                raise ValueError('Acquired source has no frozen identity')
        for path in (self.runtime / 'data/raw/genbank').glob('*.gb'):
            name = path.relative_to(self.runtime).as_posix()
            if name not in allowed:
                row = self.by_base.get(path.stem.split('.')[0])
                if row:
                    self.fail(row, 'unregistered_source_cache')
                raise ValueError('Unregistered acquisition source cache; preserve it for investigation')
            row, expected = allowed[name]
            if sha256(path) != expected:
                self.fail(row, 'source_cache_checksum_mismatch')
                raise ValueError('Acquisition source cache checksum mismatch')

    def guarded_download(self, url, path, cfg, offline=False, maximum=None):
        requested = urllib.parse.parse_qs(urllib.parse.urlparse(url).query).get('id', [''])[0]
        row = self.by_base.get(requested.split('.')[0])
        if row is None:
            raise ValueError('Acquisition requested a source outside the frozen selection')
        selected = self.selection[int(row['rank']) - 1]
        try:
            self.validate_cache()
            if row['resolved_accession']:
                source = inside(self.runtime, row['cache_path'])
            elif selected['pilot'] == 'true':
                candidates = [self.runtime / 'data/raw/genbank' / (acc + '.gb')
                              for acc in (selected['accession'], selected['frozen_pilot_accession'])]
                source = next((p for p in candidates if p.is_file()), None)
                if source is None:
                    raise ValueError('frozen_pilot_cache_missing')
            else:
                if requested != selected['accession']:
                    raise ValueError('requested_accession_differs_from_selection')
                source = self.download(url, path, cfg, offline, maximum=maximum)
            record = self.module.SeqIO.read(source, 'genbank')
            resolved = record.id
            if (not re.fullmatch(r'[A-Z]{1,6}_?\d{5,10}\.\d+', resolved)
                    or resolved.split('.')[0] != selected['accession'].split('.')[0]
                    or ('.' in selected['accession'] and resolved != selected['accession'])):
                raise ValueError('resolved_accession_version_mismatch')
            digest = sha256(source)
            if selected['pilot'] == 'true' and (resolved, digest) != (selected['frozen_pilot_accession'], selected['frozen_pilot_sha256']):
                raise ValueError('frozen_pilot_source_mismatch')
            relative = 'data/raw/genbank/' + resolved + '.gb'
            target = inside(self.runtime, relative)
            if source != target:
                if target.exists():
                    if sha256(target) != digest:
                        raise ValueError('conflicting_versioned_source_cache')
                else:
                    source.replace(target)
            row.update(status='acquired', resolved_accession=resolved, source_sha256=digest,
                       cache_path=relative, error_type='', failure_reason='')
            self.save()
            return target
        except Exception as error:
            reason = 'http_' + str(error.code) if isinstance(error, urllib.error.HTTPError) else (
                str(error) if isinstance(error, ValueError) and re.fullmatch(r'[a-z_]+', str(error)) else 'source_acquisition_failed')
            self.fail(row, reason, type(error).__name__)
            raise RuntimeError('Exploratory source acquisition stopped; inspect the private acquisition ledger') from None

    def __enter__(self):
        self.module.download = self.guarded_download
        return self

    def __exit__(self, *error):
        self.module.download = self.download


def prepare_runtime(root, prepared, plan, ledger, offline):
    import yaml
    runtime = inside(prepared, 'runtime')
    if runtime.exists():
        raise ValueError('Runtime already exists; preserve it and use --resume-from')
    cfg = copy.deepcopy(plan['effective_configuration'])
    cfg['offline'] = offline
    baseline = root / 'results/reports/resource_baseline.json'
    if baseline.exists():
        effective = json.loads(baseline.read_text(encoding='utf-8'))['effective_budget_bytes']
        cfg['resources']['effective_budget_bytes'] = min(cfg['resources'].get('effective_budget_bytes', effective), effective)
    binary = inside(root, cfg['tools']['mmseqs'])
    if not binary.is_file():
        raise FileNotFoundError('Configured MMseqs executable is unavailable')
    cfg['tools']['mmseqs'] = str(binary)
    panel = inside(root, plan['panel_root'])
    references = read_tsv(panel / 'data/manifests/reference_manifest.tsv')
    if len(references) != 39:
        raise ValueError('Expected the frozen 39-reference catalogue')
    copies = []
    for row in references:
        copies.append((reference_cache(root, row), Path('data/raw/references') / (row['source_accession'] + '.gb')))
    pilot = read_tsv(root / 'data/manifests/pilot_inputs.tsv')
    for selected, frozen in zip(ledger[:20], pilot):
        if (selected['frozen_pilot_accession'] != frozen['accession']
                or selected['frozen_pilot_sha256'] != frozen['sha256']
                or sha256(inside(root, frozen['path'])) != frozen['sha256']):
            raise ValueError('Frozen pilot cache or prepared nesting changed')
        copies.append((inside(root, frozen['path']), Path('data/raw/genbank') / (selected['accession'] + '.gb')))
    root_guard(cfg, root, sum(p.stat().st_size for p, _ in copies) + 100_000_000)
    original = pilot_hashes(root)
    shutil.copytree(root / 'src', runtime / 'src', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ('config/analysis.yaml', 'references/protein_controls.tsv',
                 'data/manifests/reference_manifest.tsv'):
        source = panel / name if name != 'config/analysis.yaml' else root / name
        target = runtime / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    for name in ('data/manifests/eligibility_snapshot.tsv', 'data/manifests/phagesdb_snapshot.json'):
        shutil.copyfile(root / name, runtime / name)
    target_driver = runtime / 'scripts/run_exploratory_screen.py'
    target_driver.parent.mkdir()
    shutil.copyfile(__file__, target_driver)
    (runtime / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    for source, relative in copies:
        target = runtime / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and sha256(target) != sha256(source):
            raise ValueError('Conflicting cached source aliases')
        shutil.copyfile(source, target)
    frozen = tree_hashes(runtime, ('src', 'config', 'references', 'scripts'),
                         ('data/manifests/eligibility_snapshot.tsv', 'data/manifests/phagesdb_snapshot.json',
                          'data/manifests/reference_manifest.tsv'))
    write_json(runtime / 'execution_inputs.json', {'status': 'exploratory_prepared',
        'preparation_sha256': sha256(prepared / 'preparation.json'), 'selection_sha256': plan['selection_sha256'],
        'driver_sha256': sha256(__file__), 'protocol_sha256': plan['protocol_sha256'],
        'planned_configuration_sha256': plan['effective_configuration_sha256'],
        'execution_configuration_sha256': digest_object(cfg), 'snapshot_sha256': frozen,
        'mmseqs_sha256': sha256(binary), 'primary_pilot_sha256': original,
        'completed_stages': [], 'public_export_performed': False})
    root_guard(cfg, root)
    return runtime


def verify_runtime(root, prepared, runtime, plan):
    inputs = json.loads((runtime / 'execution_inputs.json').read_text(encoding='utf-8'))
    if (inputs['preparation_sha256'] != sha256(prepared / 'preparation.json')
            or inputs['selection_sha256'] != plan['selection_sha256']
            or inputs['driver_sha256'] != sha256(__file__)):
        raise ValueError('Execution preparation or driver changed')
    for name, expected in inputs['snapshot_sha256'].items():
        if sha256(inside(runtime, name)) != expected:
            raise ValueError(f'Isolated execution input changed: {name}')
    if tree_hashes(runtime, ('src',)) != plan['source_tree_sha256']:
        raise ValueError('Isolated source tree differs from preparation')
    cfg = load_config(runtime / 'config/analysis.yaml')
    if digest_object(cfg) != inputs['execution_configuration_sha256']:
        raise ValueError('Isolated execution configuration changed')
    if sha256(Path(cfg['tools']['mmseqs'])) != inputs['mmseqs_sha256']:
        raise ValueError('Pinned MMseqs executable changed')
    if pilot_hashes(root) != inputs['primary_pilot_sha256']:
        raise ValueError('Original pilot artifacts changed')
    return inputs, cfg


def stage_outputs(runtime, number):
    if number == 1:
        manifest = 'data/manifests/intermediate_inputs.tsv'
        return [manifest, 'data/manifests/acquisition_ledger.tsv', 'data/manifests/eligibility_snapshot.tsv', 'data/manifests/phagesdb_snapshot.json'] + [
            r['path'] for r in read_tsv(runtime / manifest)]
    paths = list(OUTPUTS[number])
    if number == 3:
        paths += [p.relative_to(runtime).as_posix() for p in sorted((runtime / 'data/raw/references').glob('*.gb'))]
    return paths


def check_acquisition(runtime, ledger):
    manifest = read_tsv(runtime / 'data/manifests/intermediate_inputs.tsv')
    identity = lambda row: (row['phage'], row['accession'].split('.')[0])
    if [identity(r) for r in manifest] != [identity(r) for r in ledger]:
        raise ValueError('Acquired cohort differs from the frozen 50-record selection')
    for actual, selected in zip(manifest[:20], ledger[:20]):
        if (actual['accession'], actual['sha256']) != (selected['frozen_pilot_accession'], selected['frozen_pilot_sha256']):
            raise ValueError('An original pilot accession version or source checksum changed')
    private = runtime / 'data/manifests/acquisition_ledger.tsv'
    if private.exists():
        sources = read_tsv(private)
        if len(sources) != 50 or any(r['status'] != 'acquired' for r in sources):
            raise ValueError('All 50 selected sources must be acquired before this screen can complete')
        for actual, source in zip(manifest, sources):
            if (actual['accession'], actual['sha256'], actual['path']) != (
                    source['resolved_accession'], source['source_sha256'], source['cache_path']):
                raise ValueError('Final acquisition manifest differs from the immediate source ledger')


def check_controls(root, runtime, plan):
    panel = inside(root, plan['panel_root'])
    artifacts = json.loads((panel / 'artifact_manifest.json').read_text(encoding='utf-8'))['artifacts']
    for relative, phage, scope in [('results/tables/architectures.tsv', 'D29', 'high_moderate'),
                                    ('results/tables/architecture_sensitivity.tsv', 'Yvonnetastic', 'annotation')]:
        if sha256(panel / relative) != artifacts[relative]:
            raise ValueError('Frozen control architecture source changed')
        select = lambda rows: [r for r in rows if r['phage'] == phage and r['scope'] == scope and r['threshold'] == '4']
        expected, observed = select(read_tsv(panel / relative)), select(read_tsv(runtime / relative))
        keys = ('accession', 'architecture', 'confidence', 'lysin_a', 'lysin_b', 'intervening_cds', 'gap_bp',
                'lysin_a_count', 'lysin_b_count')
        if len(expected) != 1 or len(observed) != 1 or any(expected[0][k] != observed[0][k] for k in keys):
            raise ValueError(f'Frozen {phage} control architecture is contradicted')
        pairs = 'results/tables/architecture_pairs.tsv'
        if sha256(panel / pairs) != artifacts[pairs]:
            raise ValueError('Frozen control component-pair source changed')
        select_pairs = lambda rows: sorted((r for r in rows if r['phage'] == phage and r['scope'] == scope), key=digest_object)
        expected_pairs = select_pairs(read_tsv(panel / pairs))
        observed_pairs = select_pairs(read_tsv(runtime / pairs))
        if not expected_pairs or expected_pairs != observed_pairs:
            raise ValueError(f'Frozen {phage} control component pairs are contradicted')


def worker(root, workspace, resume_from):
    prepared, plan, ledger = verify_preparation(root, workspace)
    runtime = inside(prepared, 'runtime')
    if ROOT.resolve() != runtime:
        raise ValueError('Worker must execute the isolated runtime copy of this driver')
    inputs, cfg = verify_runtime(root, prepared, runtime, plan)
    os.chdir(runtime)
    receipts_path = runtime / 'results/reports/exploratory_stage_receipts.json'
    receipts = json.loads(receipts_path.read_text()) if receipts_path.exists() else {}
    for number in range(1, resume_from):
        receipt = receipts.get(str(number))
        if not receipt or receipt['execution_inputs_sha256'] != sha256(runtime / 'execution_inputs.json'):
            raise ValueError('A required exploratory dependency has no valid receipt')
        for name, expected in receipt['outputs'].items():
            if sha256(inside(runtime, name)) != expected:
                raise ValueError('A completed exploratory stage output changed')
    report = {'analysis': 'exploratory_50_genome_detection_screen', 'status': 'running',
              'mode': 'intermediate', 'started_utc': utcnow(), 'completed_stages': [],
              'public_export_performed': False, 'main_expansion_gate_modified': False,
              'primary_cohort_modified': False, 'h2_analysis_performed': False}
    try:
        for number in range(resume_from, 6):
            verify_preparation(root, workspace)
            verify_runtime(root, prepared, runtime, plan)
            receipts = {k: v for k, v in receipts.items() if int(k) < number}
            write_json(receipts_path, receipts)
            proposed = 600_000_000 if number == 3 else 0
            if number == 1 and not cfg.get('offline', False):
                manifest = runtime / 'data/manifests/intermediate_inputs.tsv'
                targets = [r['path'] for r in read_tsv(manifest)] if manifest.exists() else [
                    'data/raw/genbank/' + r['accession'] + '.gb' for r in ledger]
                proposed = sum(not inside(runtime, name).exists() for name in targets) * 3_000_000
            before = root_guard(cfg, root, proposed)
            budget_check(cfg, proposed)
            began = time.monotonic()
            print(f'Exploratory stage {number}: {STAGES[number][0]}', flush=True)
            module = importlib.import_module('actinolysis.' + STAGES[number][0])
            if number == 1:
                with AcquisitionLedger(module, runtime, ledger):
                    module.run(cfg)
            else:
                getattr(module, STAGES[number][1])(cfg)
            if number == 1:
                check_acquisition(runtime, ledger)
            if number == 5:
                check_controls(root, runtime, plan)
            outputs = {name: sha256(inside(runtime, name)) for name in stage_outputs(runtime, number)}
            after = root_guard(cfg, root)
            budget_check(cfg)
            verify_runtime(root, prepared, runtime, plan)
            receipts[str(number)] = {'status': 'exploratory_stage_complete', 'outputs': outputs,
                'seconds': round(time.monotonic() - began, 3), 'root_budget_before': before,
                'root_budget_after': after, 'proposed_additional_bytes': proposed,
                'execution_inputs_sha256': sha256(runtime / 'execution_inputs.json')}
            write_json(receipts_path, receipts)
        report.update(status='exploratory_complete_awaiting_replay', completed_stages=[1, 2, 3, 4, 5])
    except Exception as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}',
                      completed_stages=sorted(int(k) for k in receipts))
        raise
    finally:
        report['finished_utc'] = utcnow()
        write_json(runtime / 'results/reports/exploratory_execution.json', report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', required=True)
    parser.add_argument('--resume-from', type=int, choices=range(1, 6))
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--source-root', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if args.source_root is None:
            parser.error('Worker requires the original source root')
        worker(args.source_root.resolve(), args.workspace, args.resume_from or 1)
        return
    prepared, plan, ledger = verify_preparation(ROOT, args.workspace)
    runtime = inside(prepared, 'runtime')
    if args.resume_from is None:
        runtime = prepare_runtime(ROOT, prepared, plan, ledger, args.offline)
    verify_runtime(ROOT, prepared, runtime, plan)
    command = [sys.executable, '-B', str(runtime / 'scripts/run_exploratory_screen.py'),
               '--workspace', args.workspace, '--worker', '--source-root', str(ROOT),
               '--resume-from', str(args.resume_from or 1)]
    subprocess.run(command, cwd=runtime, check=True)
    print(json.dumps({'status': 'exploratory_complete_awaiting_replay',
                      'completed_stages': [1, 2, 3, 4, 5], 'public_export_performed': False}))


if __name__ == '__main__':
    main()
