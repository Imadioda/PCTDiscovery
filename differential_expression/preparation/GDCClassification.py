"""Annotate TCGA expression samples with sample-linked GDC anatomy.

This is metadata annotation, not a machine-learning classifier.
Run on a local DataFrame with one row per sample and numeric gene columns.
Only load pickle files created by your own project.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import time
import unicodedata
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

BASE = Path(__file__).resolve().parent
API = 'https://api.gdc.cancer.gov'
PHENO = ['sample', 'sample_type', '_primary_disease']
SAMPLE_RE = re.compile(r'^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}(?:[A-Z])?(?:-.*)?$')
UUID_RE = re.compile(r'^[a-fA-F0-9]{8}(?:-[a-fA-F0-9]{4}){3}-[a-fA-F0-9]{12}$')


def text_key(value):
    """Exact normalized labels; never use a substring to guess an organ."""
    value = unicodedata.normalize('NFKD', str(value or ''))
    value = ''.join(c for c in value if not unicodedata.combining(c))
    return re.sub(r'\s+', ' ', value.strip().casefold())


def clean(value):
    if value is None or pd.isna(value):
        return ''
    value = str(value).strip()
    if value.startswith('[') or text_key(value) in {'', 'unknown', 'not reported', 'not available', 'not applicable', 'none', 'nan'}:
        return ''
    return value


# Explicit anatomical vocabulary. Unlisted labels remain unmapped and are exported
# for review. No label is assigned from the patient's cancer or model predictions.
SITE_GROUPS = {
    'thyroid': ['Thyroid gland', 'Thyroid'],
    'breast': ['Breast'],
    'prostate': ['Prostate', 'Prostate gland'],
    'kidney': ['Kidney'],
    'kidney_mixed': ['Kidney; Mixed'],
    'kidney_cortex': ['Kidney; Cortex'],
    'kidney_medulla': ['Kidney; Medulla'],
    'liver': ['Liver'],
    'lung': ['Lung', 'Lung; Upper Lobe', 'Lung; Lower Lobe', 'Lung; Middle Lobe',
             'Upper lobe of lung', 'Lower lobe of lung', 'Middle lobe of lung', 'lingula'],
    'esophagus': ['Esophagus', 'Esophageal; Distal', 'Esophageal; Middle', 'Esophageal; Proximal'],
    'stomach': ['Stomach', 'gastric', 'Stomach, not otherwise specified', 'Antrum',
                'Fundus of stomach', 'Cardia; Proximal', 'Body of stomach'],
    'colon': ['Colon', 'Sigmoid colon', 'Ascending colon', 'Descending colon',
              'Transverse Colon', 'Cecum', 'Hepatic Flexure', 'Splenic Flexure'],
    'rectum': ['Rectum'],
    'rectosigmoid': ['Rectosigmoid Junction'],
    'pancreas': ['Pancreas'],
    'endometrium': ['Endometrium', 'Endometrium, NOS', 'Endometrium, Lower Segment/Isthmus'],
    'uterus_unspecified': ['Uterus', 'Corpus Uteri', 'Body of uterus'],
    'myometrium': ['Myometrium'],
    'cervix': ['Cervix', 'Cervix uteri', 'Uterine cervix', 'Endocervix', 'Ectocervix'],
    'bladder_unspecified': ['Bladder', 'Urinary bladder'],
    'bladder_urothelium': ['Bladder; macroscopic urothelial tissue'],
    'adrenal': ['Adrenal', 'Adrenal Gland'],
    'skin': ['Skin'],
    'thymus': ['Thymus'],
    'fallopian_tube': ['Fallopian tube'],
    'ovary': ['Ovary', 'Left ovary', 'Right ovary'],
    'seminal_vesicle': ['Seminal Vesicle'],
    'muscle_unspecified': ['Muscle'],
    'skeletal_muscle': ['Skeletal Muscle'],
    'spleen': ['Spleen'],
    'duodenum': ['Duodenum'],
    'lymph_node': ['Lymph node', 'Lymphnode'],
    'gastroesophageal_junction': ['Gastroesophageal Junction'],
}
SITE_LOOKUP = {text_key(label): key for key, labels in SITE_GROUPS.items() for label in labels}


def site_category(label):
    return SITE_LOOKUP.get(text_key(label), '')


class GDCClient:
    """Public API, bounded retries, cached responses, no expression upload."""
    def __init__(self, cache_dir):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)

    def download(self, endpoint, params=None):
        url = API + '/' + endpoint.lstrip('/')
        if params:
            url += '?' + urlencode(params)
        key = hashlib.sha256(url.encode()).hexdigest()
        path = self.cache / (key + '.cache')
        if path.exists():
            return path.read_bytes()
        for attempt in range(4):
            try:
                request = Request(url, headers={'User-Agent': 'TCGA-anatomy-audit/1.0'})
                with urlopen(request, timeout=90) as response:
                    data = response.read()
                tmp = path.with_suffix('.tmp')
                tmp.write_bytes(data)
                tmp.replace(path)
                (self.cache / (key + '.json')).write_text(json.dumps({
                    'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat()
                }, indent=2))
                return data
            except (HTTPError, URLError, TimeoutError) as exc:
                if isinstance(exc, HTTPError) and exc.code not in {429, 500, 502, 503, 504}:
                    raise
                if attempt == 3:
                    raise RuntimeError(f'GDC inaccessible: {endpoint}; cached work preserved.') from exc
                time.sleep(2 ** attempt)

    def json(self, endpoint, params=None):
        return json.loads(self.download(endpoint, params))


def prepare_input(df, gene_columns=None):
    df = df.copy(deep=False)
    if 'sample' not in df.columns and df.index.name == 'sample':
        df = df.reset_index()
    if df.columns.duplicated().any():
        raise ValueError('Duplicate column names: resolve them before annotation.')
    missing = set(PHENO) - set(df.columns)
    if missing:
        raise ValueError(f'Missing metadata columns: {sorted(missing)}')
    if df[PHENO].isna().any().any():
        raise ValueError('Missing sample ID, sample_type or disease.')
    df['sample'] = df['sample'].astype(str).str.strip()
    if df['sample'].duplicated().any():
        raise ValueError('Duplicate sample identifiers: do not silently duplicate observations.')
    if any(str(c).startswith(('gdc_', 'de_')) for c in df.columns):
        raise ValueError('Use the expression input before enrichment, not an already enriched dataset.')
    if gene_columns is None:
        # The project DataFrame contains three phenotype columns and gene columns.
        extra_meta = set(PHENO + ['sample_type_id', 'patient_id'])
        gene_columns = [c for c in df.columns if c not in extra_meta]
    gene_columns = list(gene_columns)
    if len(gene_columns) != len(set(gene_columns)) or not gene_columns:
        raise ValueError('Provide a nonempty, unique list of gene columns.')
    for c in gene_columns:
        if c in PHENO or c not in df or not pd.api.types.is_numeric_dtype(df[c]):
            raise ValueError(f'Invalid gene column: {c}. Pass gene_columns explicitly if needed.')
    return df.reset_index(drop=True), gene_columns


def retrieve_cases(metadata, client, workers=4):
    patients = sorted({s[:12] for s in metadata['sample'] if SAMPLE_RE.fullmatch(s)})
    batches = [patients[i:i + 40] for i in range(0, len(patients), 40)]
    def fetch(ids):
        payload = client.json('cases', {
            'filters': json.dumps({'op': 'in', 'content': {'field': 'submitter_id', 'value': ids}}),
            'expand': 'samples,diagnoses,project', 'size': len(ids)
        })
        if payload['data']['pagination']['total'] > len(ids):
            raise RuntimeError('Unexpected API pagination for case lookup.')
        return payload['data']['hits']
    cases = {}
    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(fetch, ids): ids for ids in batches}
        for i, future in enumerate(as_completed(jobs), 1):
            for case in future.result():
                cases[case['submitter_id']] = case
            print(f'GDC cases: batch {i}/{len(batches)} | {len(cases)} cases found | {(time.monotonic()-start)/60:.1f} min', flush=True)
    return cases


def retrieve_supplements(cases, client, workers=4):
    projects = sorted({c['project']['project_id'] for c in cases.values() if c.get('project')})
    manifest, offset = [], 0
    if not projects:
        return {}, [], []
    filters = {'op': 'and', 'content': [
        {'op': 'in', 'content': {'field': 'cases.project.project_id', 'value': projects}},
        {'op': 'in', 'content': {'field': 'data_format', 'value': ['BCR Biotab']}}
    ]}
    while True:
        result = client.json('files', {'filters': json.dumps(filters),
            'fields': 'file_id,file_name,data_format', 'size': 500, 'from': offset,
            'sort': 'file_id:asc'})['data']
        manifest.extend(h for h in result['hits'] if any(
            part in h['file_name'] for part in ['_ssf_normal_controls_', '_ssf_tumor_samples_']))
        offset += len(result['hits'])
        if offset >= result['pagination']['total']:
            break
        if not result['hits']:
            raise RuntimeError('Incomplete file pagination.')
    evidence = defaultdict(list)
    skipped = []
    def fetch(file):
        data = client.download('data/' + file['file_id'])
        table = pd.read_csv(io.BytesIO(data), sep='\t', skiprows=[1, 2], dtype=str, keep_default_na=False)
        if not {'bcr_sample_uuid', 'bcr_patient_barcode'}.issubset(table.columns):
            raise ValueError(f'Unexpected supplement schema: {file["file_name"]}')
        return file, table
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fetch, file) for file in manifest]
        for i, future in enumerate(as_completed(futures), 1):
            file, table = future.result()
            kind = 'normal' if '_normal_controls_' in file['file_name'] else 'tumor'
            for row in table.to_dict('records'):
                uuid = row['bcr_sample_uuid'].strip().lower()
                # Some Biotab rows aggregate multiple UUIDs with '|'. Do not guess
                # how site lists align to samples: leave them unresolved.
                if not UUID_RE.fullmatch(uuid):
                    skipped.append({'file_id': file['file_id'], 'patient': row['bcr_patient_barcode'],
                                    'raw_sample_uuid': uuid, 'reason': 'non_single_sample_uuid'})
                    continue
                row.update(source_file_id=file['file_id'], source_file_name=file['file_name'], evidence_kind=kind)
                evidence[uuid].append(row)
            print(f'GDC supplements: {i}/{len(manifest)}', flush=True)
    return evidence, manifest, skipped


def resolve_sample(barcode, case):
    samples = case.get('samples', [])
    parts = barcode.split('-')
    if len(parts[3]) == 3:  # Vial letter supplied: exact vial match.
        stem = '-'.join(parts[:4])
        matches = [s for s in samples if s.get('submitter_id') == stem]
    else:
        matches = [s for s in samples if s.get('submitter_id', '')[:15] == barcode[:15]]
    if len(matches) == 1:
        return matches[0], 'unique_sample', len(matches)
    if len(matches) > 1:
        # RNA-seq cannot originate from a diagnostic slide-only vial. Retain only
        # a UNIQUE RNA-bearing candidate; never choose the first vial in a list.
        rna_ids = case.get('submitter_analyte_ids', []) + case.get('submitter_aliquot_ids', [])
        rna = [s for s in matches if any(a.startswith(s['submitter_id'] + '-')
            and len(a.split('-')) >= 5 and a.split('-')[4].endswith('R') for a in rna_ids)]
        if len(rna) == 1:
            return rna[0], 'unique_rna_vial', len(matches)
        return None, 'ambiguous_vials', len(matches)
    return None, 'sample_not_found', 0


def extract_site(row):
    if row['evidence_kind'] == 'normal':
        main = clean(row.get('normal_tissue_anatomic_site'))
        detail = clean(row.get('normal_control_site_other'))
        site = detail if text_key(main) == 'other' else main
    else:
        main = clean(row.get('site_of_disease'))
        detail = clean(row.get('site_of_disease_description'))
        free = clean(row.get('site_of_disease_text'))
        site = free if text_key(main) == 'other' else main
        # Endometrial tumors: an organ-level "Corpus Uteri" label alone does
        # not establish an endometrial sample. Use explicit sample detail.
        if site_category(main) == 'uterus_unspecified' and site_category(detail) in {'endometrium', 'myometrium', 'cervix'}:
            site = detail
    if text_key(site) == 'other' or '|' in site:
        site = ''
    return site, main, detail


def annotate_metadata(metadata, cases, evidence):
    results = []
    for row in metadata[PHENO].to_dict('records'):
        barcode = row['sample']
        rec = {'sample': barcode, 'gdc_patient_id': barcode[:12], 'gdc_project': '',
            'gdc_case_id': '', 'gdc_sample_id': '', 'gdc_sample_barcode': '', 'gdc_sample_type': '',
            'gdc_sample_match': '', 'gdc_candidate_vials': 0, 'gdc_case_primary_site': '',
            'gdc_site_raw': '', 'gdc_site_category': '', 'gdc_site_detail': '',
            'gdc_site_status': 'invalid_barcode', 'gdc_source_urls': '', 'gdc_site_evidence': '',
            'gdc_histology': '', 'gdc_anatomy_main': '', 'gdc_api_sample_site': ''}
        if not SAMPLE_RE.fullmatch(barcode):
            results.append(rec); continue
        case = cases.get(barcode[:12])
        if not case:
            rec['gdc_site_status'] = 'case_not_found'; results.append(rec); continue
        rec.update(gdc_project=case.get('project', {}).get('project_id', ''),
                   gdc_case_id=case['case_id'], gdc_case_primary_site=case.get('primary_site', ''))
        sample, method, count = resolve_sample(barcode, case)
        rec.update(gdc_sample_match=method, gdc_candidate_vials=count)
        if sample is None:
            rec['gdc_site_status'] = method; results.append(rec); continue
        rec.update(gdc_sample_id=sample['sample_id'], gdc_sample_barcode=sample['submitter_id'],
                   gdc_sample_type=sample.get('sample_type', ''),
                   gdc_api_sample_site=clean(sample.get('biospecimen_anatomic_site')))
        if sample.get('sample_type') != row['sample_type']:
            rec['gdc_site_status'] = 'sample_type_conflict'; results.append(rec); continue
        kind = {'Solid Tissue Normal': 'normal', 'Primary Tumor': 'tumor'}.get(row['sample_type'])
        if not kind:
            rec['gdc_site_status'] = 'sample_type_not_in_initial_scope'; results.append(rec); continue
        candidates = [e for e in evidence.get(sample['sample_id'].lower(), [])
                      if e['bcr_patient_barcode'] == barcode[:12] and e['evidence_kind'] == kind]
        rec['gdc_source_urls'] = '; '.join(sorted({API + '/data/' + e['source_file_id'] for e in candidates}))
        labels = []
        raw_evidence = []
        for e in candidates:
            site, main, detail = extract_site(e)
            labels.append(site)
            raw_evidence.append({'site': site, 'main': main, 'detail': detail,
                'site_text': clean(e.get('site_of_disease_text')), 'file_id': e['source_file_id']})
        rec['gdc_site_evidence'] = json.dumps(raw_evidence, ensure_ascii=False)
        rec['gdc_anatomy_main'] = '; '.join(sorted({x['main'] for x in raw_evidence if x['main']}))
        rec['gdc_site_detail'] = '; '.join(sorted({x['detail'] for x in raw_evidence if x['detail']}))
        rec['gdc_histology'] = '; '.join(sorted({clean(e.get('histologic_diagnosis')) for e in candidates
                                                if clean(e.get('histologic_diagnosis'))}))
        labels = sorted({label for label in labels if label})
        # Site evidence only from specimen-linked supplements; diagnosis-level
        # primary_site is retained for auditing, NEVER used as an anatomy fallback.
        keys = {site_category(label) or text_key(label) for label in labels}
        if len(keys) > 1:
            rec['gdc_site_status'] = 'conflicting_site_annotations'
            rec['gdc_site_raw'] = '; '.join(labels)
        elif labels:
            rec['gdc_site_raw'] = labels[0]
            rec['gdc_site_category'] = site_category(labels[0])
            api_site = rec['gdc_api_sample_site']
            # An explicit conflicting sample-level API label triggers review.
            if api_site and site_category(api_site) and rec['gdc_site_category'] and site_category(api_site) != rec['gdc_site_category'] and not (
                site_category(api_site) == 'uterus_unspecified' and rec['gdc_site_category'] in {'endometrium', 'myometrium', 'cervix'}):
                rec['gdc_site_status'] = 'api_supplement_site_conflict'
            else:
                rec['gdc_site_status'] = 'documented' if rec['gdc_site_category'] else 'documented_unmapped_label'
        else:
            rec['gdc_site_status'] = 'site_not_reported' if candidates else 'no_sample_linked_supplement'
        results.append(rec)
    return pd.DataFrame(results)


def classer_gdc(dtf, output_dir=BASE / 'resultats_GDC', gene_columns=None, workers=4):
    """Return the full enriched DataFrame and persist it for the second script."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if workers < 1 or workers > 8:
        raise ValueError('workers must be between 1 and 8.')
    df, genes = prepare_input(dtf, gene_columns)
    client = GDCClient(output_dir / 'cache')
    start = time.monotonic()
    print(f'Input: {len(df)} samples | {len(genes)} genes', flush=True)
    status = client.json('status')
    cases = retrieve_cases(df[PHENO], client, workers)
    evidence, manifest, skipped = retrieve_supplements(cases, client, workers)
    annotations = annotate_metadata(df[PHENO], cases, evidence)
    enriched = df.merge(annotations, on='sample', how='left', validate='one_to_one', sort=False)
    if enriched['sample'].tolist() != df['sample'].tolist():
        raise RuntimeError('Unexpected sample-order change.')
    annotations.merge(df[PHENO], on='sample', validate='one_to_one').to_csv(output_dir / 'annotations_GDC.csv', index=False)
    enriched.to_pickle(output_dir / 'dataset_enrichi.pkl')
    (output_dir / 'genes.json').write_text(json.dumps(genes, indent=2))
    (output_dir / 'manifest_supplements.json').write_text(json.dumps(manifest, indent=2))
    pd.DataFrame(skipped, columns=['file_id', 'patient', 'raw_sample_uuid', 'reason']).to_csv(output_dir / 'supplement_rows_unresolved.csv', index=False)
    audit = enriched.groupby(['gdc_project', 'sample_type', 'gdc_site_status', 'gdc_site_category'], dropna=False).size().reset_index(name='n_samples')
    audit.to_csv(output_dir / 'bilan_annotations.csv', index=False)
    enriched.loc[enriched['gdc_site_status'] == 'documented_unmapped_label',
        ['gdc_project', 'sample_type', 'gdc_site_raw', 'gdc_site_detail']].drop_duplicates().to_csv(output_dir / 'sites_a_revoir.csv', index=False)
    (output_dir / 'run_info.json').write_text(json.dumps({
        'created_utc': datetime.now(timezone.utc).isoformat(), 'gdc_status': status,
        'n_samples': len(df), 'n_genes': len(genes), 'workers': workers,
        'rule': 'sample_uuid + patient_id, no diagnosis-based anatomy imputation',
        'input_expression_unchanged': True,
        'site_vocabulary': SITE_GROUPS,
    }, indent=2))
    print(annotations['gdc_site_status'].value_counts().to_string())
    print(f'Saved: {output_dir / "dataset_enrichi.pkl"} | {(time.monotonic()-start)/60:.1f} min', flush=True)
    return enriched


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=BASE / 'dtf_modele.pkl')
    parser.add_argument('--output', type=Path, default=BASE / 'resultats_GDC')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f'Input not found: {args.input}. Save dtf_modele.to_pickle(...) first.')
    # This pickle is a local intermediate produced by your own trusted project.
    classer_gdc(pd.read_pickle(args.input), args.output, workers=args.workers)
