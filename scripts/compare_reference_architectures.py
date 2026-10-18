#!/usr/bin/env python3
"""Export a source-survey replay and compare every frozen scientific artifact."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('reference_survey',ROOT/'scripts/survey_reference_architectures.py')
survey=importlib.util.module_from_spec(spec)
spec.loader.exec_module(survey)

def run(args):
    runtime=ROOT/'.experiments'/args.workspace
    target=ROOT/'.experiments'/args.output
    for path in (runtime,target):
        if not path.resolve().is_relative_to((ROOT/'.experiments').resolve()): raise ValueError('Comparison path is outside .experiments')
    frozen=json.loads((runtime/'survey_inputs.json').read_text())
    if survey.sha256(survey.PROTOCOL)!=frozen['protocol_sha256'] or survey.sha256(ROOT/'scripts/survey_reference_architectures.py')!=frozen['driver_sha256']:
        raise ValueError('Survey protocol or driver changed')
    for name,digest in frozen['snapshot_sha256'].items():
        if survey.sha256(runtime/name)!=digest: raise ValueError('Survey source snapshot changed')
    survey.summarize(runtime,target)
    expected=json.loads((survey.PUBLIC/'artifact_manifest.json').read_text())['artifacts']
    observed=json.loads((target/'artifact_manifest.json').read_text())['artifacts']
    if expected != observed: raise ValueError('Scientific source-survey artifacts differ on replay')
    for name,digest in expected.items():
        if survey.sha256(target/name)!=digest or survey.sha256(survey.PUBLIC/name)!=digest: raise ValueError('Scientific artifact checksum mismatch')
    report={'status':'pass','scientific_artifacts_compared':len(expected),'all_artifacts_match':True,
            'reference_artifact_manifest_sha256':survey.sha256(survey.PUBLIC/'artifact_manifest.json'),
            'artifacts':expected,'comparison_script_sha256':survey.sha256(__file__),
            'scope':'Fresh scientific execution from the same checksum-verified acquired source records; no independent software reinstall or source redownload'}
    report_path=ROOT/args.report if args.report else target/'reproduction_comparison.json'
    survey.write_json(report_path,report)
    print(json.dumps({'status':'pass','scientific_artifacts_compared':len(expected)}))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--report')
    run(parser.parse_args())
