"""Export only predefined, site-compatible TCGA tumor/normal comparisons.

This script prepares expression matrices for R. It does NOT run limma.
The minimum sample count is a configurable planning rule, not a power analysis.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from GDCClassification import BASE, PHENO

# Explicit initial scope: not a Cartesian product of organs and cancers.
# Distinct kidney and lung cancers remain distinct comparisons.
# Broader/ambiguous tissues are deliberately not mapped to compatible controls.
DEFAULT_RULES = [
    ('THCA_thyroid', 'TCGA-THCA', 'thyroid carcinoma', 'thyroid', 'thyroid'),
    ('BRCA_breast', 'TCGA-BRCA', 'breast invasive carcinoma', 'breast', 'breast'),
    ('PRAD_prostate', 'TCGA-PRAD', 'prostate adenocarcinoma', 'prostate', 'prostate'),
    ('KIRC_kidney', 'TCGA-KIRC', 'kidney clear cell carcinoma', 'kidney', 'kidney'),
    ('KIRP_kidney', 'TCGA-KIRP', 'kidney papillary cell carcinoma', 'kidney', 'kidney'),
    ('KICH_kidney', 'TCGA-KICH', 'kidney chromophobe', 'kidney', 'kidney'),
    ('LUAD_lung', 'TCGA-LUAD', 'lung adenocarcinoma', 'lung', 'lung'),
    ('LUSC_lung', 'TCGA-LUSC', 'lung squamous cell carcinoma', 'lung', 'lung'),
    ('LIHC_liver', 'TCGA-LIHC', 'liver hepatocellular carcinoma', 'liver', 'liver'),
    ('COAD_colon', 'TCGA-COAD', 'colon adenocarcinoma', 'colon', 'colon'),
    ('READ_rectum', 'TCGA-READ', 'rectum adenocarcinoma', 'rectum', 'rectum'),
    ('ESCA_esophagus', 'TCGA-ESCA', 'esophageal carcinoma', 'esophagus', 'esophagus'),
    ('STAD_stomach', 'TCGA-STAD', 'stomach adenocarcinoma', 'stomach', 'stomach'),
    ('PAAD_pancreas', 'TCGA-PAAD', 'pancreatic adenocarcinoma', 'pancreas', 'pancreas'),
    ('UCEC_endometrium', 'TCGA-UCEC', 'uterine corpus endometrioid carcinoma', 'endometrium', 'endometrium'),
    ('CESC_cervix', 'TCGA-CESC', 'cervical & endocervical cancer', 'cervix', 'cervix'),
    ('BLCA_bladder', 'TCGA-BLCA', 'bladder urothelial carcinoma', 'bladder_unspecified|bladder_urothelium', 'bladder_urothelium'),
]
RULE_COLUMNS = ['comparison_id', 'project', 'disease', 'tumor_sites', 'normal_sites']


def prepare_comparisons(input_path=BASE / 'resultats_GDC' / 'dataset_enrichi.pkl',
                        genes_path=None, output_dir=BASE / 'resultats_DE_preparation',
                        min_normals=10, min_tumors=10):
    input_path, output_dir = Path(input_path).resolve(), Path(output_dir).resolve()
    genes_path = Path(genes_path) if genes_path else input_path.parent / 'genes.json'
    if min_normals < 3 or min_tumors < 3:
        raise ValueError('This initial workflow requires at least 3 independent patients per group; even 3 is not a guarantee of reliability.')
    # Use a fresh run directory to prevent old exports surviving a stricter rerun.
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f'{output_dir} is nonempty. Use a new --output directory for this run.')
    df = pd.read_pickle(input_path)
    genes = json.loads(genes_path.read_text())
    needed = PHENO + ['gdc_project', 'gdc_site_status', 'gdc_site_category', 'gdc_patient_id']
    if not set(needed).issubset(df.columns) or not set(genes).issubset(df.columns):
        raise ValueError('Input or genes.json does not match the GDC annotation output.')
    if not genes or len(genes) != len(set(genes)) or df['sample'].duplicated().any():
        raise ValueError('Empty/duplicate genes or duplicate samples.')
    if set(genes) & set(needed) or any(str(c).startswith(('gdc_', 'de_')) for c in genes):
        raise ValueError('Metadata must never be included in gene expression.')
    if not all(pd.api.types.is_numeric_dtype(df[g]) for g in genes):
        raise ValueError('Non-numeric expression columns.')
    output_dir.mkdir(parents=True, exist_ok=True)
    df = df.reset_index(drop=True)
    rules = pd.DataFrame(DEFAULT_RULES, columns=RULE_COLUMNS)
    rules['min_normals'] = min_normals
    rules['min_tumors'] = min_tumors
    rules.to_csv(output_dir / 'comparaisons_prevues.csv', index=False)
    meta_cols = [c for c in df.columns if c in PHENO or str(c).startswith('gdc_')]
    audit = df[meta_cols].copy()
    audit['de_group'] = df['sample_type'].map({'Primary Tumor': 'Tumor', 'Solid Tissue Normal': 'Normal'}).fillna('')
    audit['de_comparison_id'] = ''
    audit['de_decision'] = 'excluded'
    audit['de_reason'] = 'project_outside_initial_scope'
    audit.loc[audit['de_group'] == '', 'de_reason'] = 'sample_type_outside_initial_scope'
    plan = []
    included_indices = []

    for rule in rules.to_dict('records'):
        name = rule['comparison_id']
        in_project = df['gdc_project'].eq(rule['project'])
        audit.loc[in_project, 'de_comparison_id'] = name
        audit.loc[in_project, 'de_reason'] = 'site_not_documented_or_compatible'
        same_disease = df['_primary_disease'].eq(rule['disease'])
        audit.loc[in_project & ~same_disease, 'de_reason'] = 'project_disease_mismatch'
        supported_type = audit['de_group'].ne('')
        documented = df['gdc_site_status'].eq('documented')
        base_mask = in_project & same_disease & documented & supported_type
        tumor_mask = base_mask & df['sample_type'].eq('Primary Tumor') & df['gdc_site_category'].isin(rule['tumor_sites'].split('|'))
        normal_mask = base_mask & df['sample_type'].eq('Solid Tissue Normal') & df['gdc_site_category'].isin(rule['normal_sites'].split('|'))
        selected = audit.loc[tumor_mask | normal_mask].copy()
        # One observation per patient and condition. Repeated samples are left
        # for review instead of selected arbitrarily or treated as independent.
        duplicates = selected.duplicated(['gdc_patient_id', 'de_group'], keep=False)
        audit.loc[selected.index[duplicates], 'de_reason'] = 'multiple_samples_same_patient_and_group'
        selected = selected.loc[~duplicates]
        audit.loc[selected.index, 'de_reason'] = 'candidate_for_comparison'
        normals = set(selected.loc[selected['de_group'] == 'Normal', 'gdc_patient_id'])
        tumors = set(selected.loc[selected['de_group'] == 'Tumor', 'gdc_patient_id'])
        record = {**rule, 'n_normal': len(normals), 'n_tumor': len(tumors),
                  'n_paired_patients': len(normals & tumors), 'n_genes_input': len(genes),
                  'n_genes_exported': 0, 'status': '', 'expression_file': '', 'samples_file': ''}
        if len(normals) < min_normals or len(tumors) < min_tumors:
            record['status'] = 'skipped_insufficient_patients'
            audit.loc[selected.index, 'de_reason'] = record['status']
            plan.append(record)
            print(f'{name}: skipped | normal={len(normals)} tumor={len(tumors)}', flush=True)
            continue
        selected = selected.sort_values(['de_group', 'gdc_patient_id', 'sample'])
        x = df.loc[selected.index, genes].to_numpy(dtype=float)
        finite = np.isfinite(x).all(axis=0)
        varying = np.zeros(len(genes), dtype=bool)
        varying[finite] = np.ptp(x[:, finite], axis=0) > 0
        retained = finite & varying
        keep_genes = [g for g, keep in zip(genes, retained) if keep]
        if not keep_genes:
            record['status'] = 'skipped_no_usable_genes'
            audit.loc[selected.index, 'de_reason'] = record['status']
            plan.append(record)
            continue
        directory = output_dir / name
        directory.mkdir()
        removed = pd.DataFrame({'gene': genes, 'reason': np.where(~finite, 'missing_or_nonfinite_within_comparison',
                                                    np.where(~varying, 'constant_within_comparison', 'retained'))})
        removed.loc[~retained].to_csv(directory / 'genes_excluded.csv', index=False)
        # R-friendly orientation: genes in rows, samples in columns. Values are
        # unchanged log expression; no StandardScaler, imputation or averaging.
        expression = pd.DataFrame(x[:, retained].T, index=keep_genes,
                                  columns=selected['sample'].tolist())
        expression.index.name = 'gene'
        expression.to_csv(directory / 'expression.tsv.gz', sep='\t',
                          compression={'method': 'gzip', 'compresslevel': 1})
        selected = selected.rename(columns={'de_group': 'group', 'gdc_patient_id': 'patient_id'})
        selected = selected.drop(columns=['de_decision', 'de_reason'])
        selected['is_paired_patient'] = selected['patient_id'].isin(normals & tumors)
        if selected['sample'].tolist() != expression.columns.tolist():
            raise RuntimeError('Expression/metadata alignment error.')
        selected.to_csv(directory / 'samples.csv', index=False)
        record.update(status='exported_for_expression_qc', n_genes_exported=len(keep_genes),
                      expression_file=f'{name}/expression.tsv.gz', samples_file=f'{name}/samples.csv')
        (directory / 'comparison.json').write_text(json.dumps({**record,
            'contrast': 'Tumor minus Normal',
            'expression_scale_expected': 'original Xena log2(norm_value + 1), not StandardScaler output',
            'missing_gene_policy': 'remove nonfinite genes within this comparison; no imputation',
            'pairing': 'patient_id and is_paired_patient are provided; account for paired samples in R',
            'interpretation': 'site-compatible candidate comparison, still requires expression QC and assessment of confounding',
        }, indent=2))
        audit.loc[selected.index, 'de_decision'] = 'included'
        audit.loc[selected.index, 'de_reason'] = 'exported_for_expression_qc'
        included_indices.extend(selected.index.tolist())
        plan.append(record)
        print(f'{name}: exported | normal={len(normals)} tumor={len(tumors)} genes={len(keep_genes)} paired={len(normals & tumors)}', flush=True)
        del x, expression

    audit.to_csv(output_dir / 'decisions_echantillons.csv', index=False)
    plan_df = pd.DataFrame(plan)
    plan_df.to_csv(output_dir / 'plan_analyses.csv', index=False)
    # Filtered union for Python inspection. R must use the per-comparison exports,
    # which additionally have comparison-specific gene filtering.
    indices = sorted(set(included_indices))
    filtered = df.loc[indices].copy()
    filtered['de_comparison_id'] = audit.loc[indices, 'de_comparison_id'].to_numpy()
    filtered['de_group'] = audit.loc[indices, 'de_group'].to_numpy()
    filtered.to_pickle(output_dir / 'dataset_DE_filtre.pkl')
    (output_dir / 'run_info.json').write_text(json.dumps({
        'created_utc': datetime.now(timezone.utc).isoformat(), 'input': str(input_path),
        'n_input_samples': len(df), 'n_retained_samples': len(filtered),
        'min_normals': min_normals, 'min_tumors': min_tumors,
        'note': 'Sample-count thresholds are planning choices, not a statistical power calculation.'
    }, indent=2))
    print(f'Finished: {len(filtered)}/{len(df)} samples retained. Results: {output_dir}', flush=True)
    return filtered, plan_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=BASE / 'resultats_GDC' / 'dataset_enrichi.pkl')
    parser.add_argument('--genes', type=Path, default=None)
    parser.add_argument('--output', type=Path, default=BASE / 'resultats_DE_preparation')
    parser.add_argument('--min-normals', type=int, default=10)
    parser.add_argument('--min-tumors', type=int, default=10)
    args = parser.parse_args()
    prepare_comparisons(args.input, args.genes, args.output, args.min_normals, args.min_tumors)
