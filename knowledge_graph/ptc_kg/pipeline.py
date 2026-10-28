"""Independent batch pipeline; imports no v1 training or Streamlit modules."""
import json
import platform
import sys
from pathlib import Path

from . import __version__
from .common import digest, file_hash, read_csv, read_json, utc_now, write_json, write_table
from .pubtator import CollectionError, PubTatorClient, documents
from .graph import KnowledgeGraph, validate_graph

DEFAULT_EXCLUSIONS = Path(__file__).resolve().parent / 'exclusions.json'


def source_inputs(rag_data):
    if not rag_data:
        return {}
    root = Path(rag_data)
    if not (root / 'retained_genes_annotated.csv').is_file():
        raise ValueError('rag_data must contain retained_genes_annotated.csv; omit --rag-data for publication-only mode')
    files = [root / 'literature_all_cancers' / 'articles.json',
             root / 'literature_all_cancers' / 'gene_article_links.csv',
             root / 'retained_genes_annotated.csv',
             root / 'references' / 'hgnc_complete_set.txt']
    return {str(p.relative_to(root)): file_hash(p) for p in files if p.is_file()}


def select_publications(rag_data=None, comparison=None, pmids=None, exclusions=DEFAULT_EXCLUSIONS):
    explicit = set(str(p) for p in (pmids or []))
    if any(not p.isdigit() for p in explicit):
        raise ValueError('Invalid explicit PMID')
    blacklist = read_json(exclusions)
    excluded = {(r['comparison_id'], r['approved_symbol'], r['article_id'])
                for r in blacklist['excluded_links']}
    contexts, articles_by_id = {}, {}
    excluded_publications = set()
    if rag_data:
        root = Path(rag_data) / 'literature_all_cancers'
        articles = read_json(root / 'articles.json')
        ids = [a['article_id'] for a in articles]
        if len(ids) != len(set(ids)):
            raise ValueError('Duplicate article IDs in corpus')
        articles_by_id = {a['article_id']: a for a in articles}
        for article in articles:
            labels = article.get('publication_types', []) + article.get('publication_flags', [])
            if any('retract' in label.lower() or 'expression of concern' in label.lower() for label in labels):
                excluded_publications.add(article['article_id'])
        for row in read_csv(root / 'gene_article_links.csv'):
            if not {'comparison_id', 'approved_symbol', 'article_id'}.issubset(row):
                raise ValueError('Missing corpus-link columns')
            if row['article_id'] not in articles_by_id:
                raise ValueError('Corpus link refers to a missing article')
            if comparison and row['comparison_id'] != comparison:
                continue
            if (row['comparison_id'], row['approved_symbol'], row['article_id']) in excluded:
                continue
            if row['article_id'] in excluded_publications:
                continue
            if not row['article_id'].startswith('MED:') or not row['article_id'][4:].isdigit():
                continue
            pmid = row['article_id'][4:]
            contexts.setdefault(pmid, []).append({
                'comparison_id': row['comparison_id'], 'original_gene_label': row.get('gene', ''),
                'approved_symbol': row['approved_symbol'], 'hgnc_id': row.get('hgnc_id', ''),
                'association_review_status': row.get('review_status', 'unknown'),
                'interpretation': 'corpus_selection_context_only_not_extraction_context',
            })
    blocked_pmids = {x[4:] for x in excluded_publications if x.startswith('MED:')}
    if explicit & blocked_pmids:
        raise ValueError('Explicit PMIDs include a publication flagged in the existing corpus')
    # Explicit seeds override automatic corpus selection, but retain known context.
    chosen = explicit if explicit else set(contexts)
    return sorted(chosen, key=int), contexts, {
        'excluded_flagged_article_ids': sorted(excluded_publications),
        'exclusions_sha256': file_hash(exclusions),
        'explicit_seed_selection': bool(explicit),
    }


def collect(out, rag_data=None, comparison=None, pmids=None, query=None,
            pages=1, limit=10, offline=False, client=None):
    if limit < 1 or limit > 10000:
        raise ValueError('Use a limit between 1 and 10000')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    selected, contexts, policy = select_publications(rag_data, comparison, pmids)
    client = client or PubTatorClient(out / 'cache', offline=offline)
    search_records = []
    if query:
        additional, search_records = client.search(query, pages)
        blocked = {x[4:] for x in policy['excluded_flagged_article_ids'] if x.startswith('MED:')}
        selected = sorted((set(selected) | set(additional)) - blocked, key=int)
    if not selected:
        raise ValueError('No publications selected. Provide rag_data, explicit PMIDs or a query.')
    selected = selected[:limit]
    records = []
    for i, pmid in enumerate(selected, 1):
        print(f'[{i}/{len(selected)}] PubTator PMID {pmid}', flush=True)
        item = {'pmid': pmid, 'corpus_context': contexts.get(pmid, [])}
        try:
            doc, raw, path = client.publication(pmid)
            item.update(status='collected', cache_file=str(path.relative_to(out)),
                        payload_sha256=raw['payload_sha256'], url=raw['url'],
                        retrieved_at_utc=raw['retrieved_at_utc'])
        except (CollectionError, ValueError) as exc:
            item.update(status='failed', error=str(exc))
        records.append(item)
        # A stopped process retains successful cache entries and progress.
        write_json(out / 'collection_progress.json', records)
    manifest = {
        'schema_version': '0.1.0', 'created_at_utc': utc_now(), 'comparison_filter': comparison,
        'query': query, 'search_pages': pages if query else 0, 'search_records': search_records,
        'limit': limit, 'selection_rule': 'unique PMIDs sorted numerically then limited',
        'selection_policy': policy, 'input_sha256': source_inputs(rag_data),
        'publications': records, 'collected': sum(r['status'] == 'collected' for r in records),
        'failed': sum(r['status'] == 'failed' for r in records),
        'offline': offline,
    }
    write_json(out / 'collection_manifest.json', manifest)
    return manifest


def build(collection_dir, output_dir, rag_data=None, allow_partial=False):
    collection_dir, output_dir = Path(collection_dir), Path(output_dir)
    manifest = read_json(collection_dir / 'collection_manifest.json')
    if manifest['failed'] and not allow_partial:
        raise ValueError('Collection is incomplete. Resume collection or explicitly use --allow-partial.')
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError('Graph output directory must be empty; choose a new output path.')
    actual_inputs = source_inputs(rag_data)
    if rag_data and manifest['input_sha256'] != actual_inputs:
        raise ValueError('RAG inputs changed since collection. Use the same snapshot or recollect.')
    graph = KnowledgeGraph()
    for record in manifest['publications']:
        if record['status'] != 'collected':
            continue
        cache_path = (collection_dir / record['cache_file']).resolve()
        if not cache_path.is_relative_to(collection_dir.resolve()):
            raise ValueError('Cache path escapes collection directory')
        raw = read_json(cache_path)
        if digest(raw['payload']) != record['payload_sha256'] or raw['payload_sha256'] != record['payload_sha256']:
            raise ValueError('Cache/manifest hash mismatch')
        matches = [doc for doc in documents(raw['payload']) if str(doc.get('pmid') or doc['id']) == record['pmid']]
        if len(matches) != 1:
            raise ValueError('Cached document PMID does not match selection')
        graph.add_document(matches[0], {
            'url': record['url'], 'retrieved_at_utc': record['retrieved_at_utc'],
            'payload_sha256': record['payload_sha256'], 'cache_file': record['cache_file'],
        }, record['corpus_context'])
    if rag_data:
        root = Path(rag_data)
        results = root / 'retained_genes_annotated.csv'
        reference = root / 'references' / 'hgnc_complete_set.txt'
        if results.is_file():
            graph.add_expression(results, manifest['comparison_filter'],
                                 reference if reference.is_file() else None)
    payload = graph.payload()
    report = validate_graph(payload)
    if not report['valid']:
        raise ValueError('Graph integrity failure: ' + str(report['errors']))
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / 'graph.json', payload)
    fields = {
        'nodes': ['id', 'type', 'namespace', 'labels', 'xrefs', 'attributes'],
        'edges': ['id', 'subject', 'predicate', 'object', 'status', 'directed', 'evidence_ids'],
        'evidence': ['id', 'edge_id', 'publication_id', 'source', 'review_status', 'extraction_score',
                     'source_relation_id', 'role_order', 'original_relation', 'provenance', 'corpus_context',
                     'context_passages', 'supporting_annotation_ids', 'support_level', 'biological_context_status',
                     'original_gene_label', 'hgnc_id', 'effect_log2_scale', 'fdr_bh_treat', 'source_row', 'cohort_details'],
        'mentions': ['pmid', 'annotation_id', 'entity_id', 'text', 'locations', 'passage_index', 'span_verified', 'raw_infons'],
        'issues': ['reason', 'pmid', 'annotation_id', 'relation_id', 'entity_id', 'detail', 'gene', 'hgnc_id', 'comparison_id', 'infons'],
    }
    for name, columns in fields.items():
        write_table(output_dir / f'{name}.csv', payload[name], columns)
    review_rows = []
    evidence_by_id = {x['id']: x for x in payload['evidence']}
    for edge in payload['edges']:
        if edge['status'] != 'automatically_extracted':
            continue
        for evidence_id in edge['evidence_ids']:
            ev = evidence_by_id[evidence_id]
            review_rows.append({
                'edge_id': edge['id'], 'evidence_id': evidence_id,
                'subject': edge['subject'], 'predicate': edge['predicate'], 'object': edge['object'],
                'publication_id': ev['publication_id'], 'extraction_score': ev['extraction_score'],
                'source_relation_id': ev['source_relation_id'], 'role_order': ev['role_order'],
                'human_label': '', 'reviewer': '', 'review_notes': '',
            })
    write_table(output_dir / 'review_queue.csv', review_rows,
                ['edge_id', 'evidence_id', 'subject', 'predicate', 'object', 'publication_id',
                 'extraction_score', 'source_relation_id', 'role_order', 'human_label', 'reviewer', 'review_notes'])
    write_json(output_dir / 'validation_report.json', report)
    code_root = Path(__file__).parent
    write_json(output_dir / 'build_manifest.json', {
        'created_at_utc': utc_now(), 'software': 'ptc_kg', 'version': __version__,
        'python': sys.version, 'platform': platform.platform(),
        'collection_manifest_sha256': file_hash(collection_dir / 'collection_manifest.json'),
        'input_sha256': actual_inputs,
        'code_sha256': {p.name: file_hash(p) for p in sorted(code_root.glob('*.py'))},
        'graph_sha256': file_hash(output_dir / 'graph.json'),
        'partial_collection': bool(manifest['failed']), 'validation': report,
        'limitations': [
            'Graph integrity validation is not scientific validation.',
            'PubTator gene/protein mentions normalize to NCBI Gene; proteins/isoforms are not inferred.',
            'Corpus comparison labels are retrieval context, not biological context of every extracted relation.',
            'No causal direction inferred from BioREx role order or correlation.',
            'Extraction scores are not calibrated biological probabilities.',
            'No co-mention is converted into a biological edge.',
            'No prediction, curated interaction or drug approval is fabricated.',
            'Retraction flags are from the supplied corpus snapshot, not a live comprehensive check.',
            'Graph not ready for GNN training: split design, feature selection and label review still required.',
        ],
    })
    return report
