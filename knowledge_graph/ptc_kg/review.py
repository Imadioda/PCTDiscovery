"""Evidence-level review. Human agreement never changes biological truth status."""
import html
import json
from pathlib import Path
from .common import read_json, read_csv, write_json, write_table, file_hash, utc_now
from .graph import validate_graph

FIELDS = ['graph_sha256', 'edge_id', 'evidence_id', 'publication_id',
          'subject', 'predicate', 'object', 'decision', 'reviewer',
          'support_quote', 'source_locator', 'biological_context', 'notes']
DECISIONS = {'supported', 'unsupported', 'uncertain'}


def load_graph(path):
    graph = read_json(path)
    report = validate_graph(graph)
    if not report['valid']:
        raise ValueError('Invalid graph; run validate first')
    return graph


def rows_for(graph, sha):
    edges = {e['id']: e for e in graph['edges']}
    rows = []
    for evidence in graph['evidence']:
        edge = edges[evidence['edge_id']]
        if edge['status'] != 'automatically_extracted':
            continue
        row = dict.fromkeys(FIELDS, '')
        row.update(graph_sha256=sha, edge_id=edge['id'], evidence_id=evidence['id'],
                   publication_id=evidence['publication_id'], subject=edge['subject'],
                   predicate=edge['predicate'], object=edge['object'])
        rows.append(row)
    return sorted(rows, key=lambda r: (r['edge_id'], r['evidence_id']))


def new_directory(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError('Output directory must be new or empty')
    path.mkdir(parents=True, exist_ok=True)
    return path


def export_review(graph_path, output):
    graph = load_graph(graph_path)
    sha = file_hash(graph_path)
    rows = rows_for(graph, sha)
    evidence = {e['id']: e for e in graph['evidence']}
    labels = {n['id']: ', '.join(n['labels']) for n in graph['nodes']}
    esc = lambda x: html.escape(str(x), quote=True)
    cards = []
    for row in rows:
        ev = evidence[row['evidence_id']]
        pmid = row['publication_id'].removeprefix('PMID:')
        if not pmid.isdigit():
            raise ValueError('Invalid PMID')
        cards.append('<article><h2>' + esc(labels.get(row['subject']) or row['subject'])
                     + ' / ' + esc(row['predicate']) + ' / '
                     + esc(labels.get(row['object']) or row['object']) + '</h2>'
                     + '<p>' + esc(row['edge_id']) + '<br>' + esc(row['evidence_id']) + '</p>'
                     + f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank" rel="noopener noreferrer">PubMed {pmid}</a>'
                     + '<p>Document-level extraction; passages below are context, not verified supporting sentences.</p>'
                     + '<details open><summary>Available context passages</summary><pre>'
                     + esc(json.dumps(ev.get('context_passages', []), ensure_ascii=False, indent=2))
                     + '</pre></details><details><summary>Raw relation and provenance</summary><pre>'
                     + esc(json.dumps(ev, ensure_ascii=False, indent=2)) + '</pre></details></article>')
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>PTCDiscovery evidence review</title>
<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:20px;background:#f5f7fa}article{background:white;padding:24px;margin:20px 0;border:1px solid #ccd4df}pre{white-space:pre-wrap;overflow-wrap:anywhere}h2{font-size:20px}</style>
<h1>Review extracted relations</h1><p>Edit decisions.csv alongside this read-only report. Decisions: supported, unsupported, uncertain; blank means unreviewed. For every decision supply a reviewer and notes. Supported decisions also require a short source quote, source locator and biological context. Agreement with an extraction is not proof of causality or therapeutic efficacy. An unsupported extraction is not a biological negative for GNN training.</p>'''
    output = new_directory(output)
    write_table(output / 'decisions.csv', rows, FIELDS)
    (output / 'review.html').write_text(page + ''.join(cards) + '</html>', encoding='utf-8')
    result = {'valid': True, 'graph_sha256': sha, 'review_items': len(rows)}
    write_json(output / 'review_manifest.json', result)
    return result


def import_review(graph_path, decisions_path, output):
    graph = load_graph(graph_path)
    sha = file_hash(graph_path)
    expected = {r['evidence_id']: r for r in rows_for(graph, sha)}
    rows = read_csv(decisions_path)
    seen, reviews = set(), []
    immutable = FIELDS[:7]
    for row in rows:
        ident = row.get('evidence_id')
        if ident not in expected or ident in seen:
            raise ValueError(f'Unknown or duplicate evidence: {ident}')
        seen.add(ident)
        if any(row.get(k) != expected[ident][k] for k in immutable):
            raise ValueError(f'Graph identity or immutable columns changed: {ident}')
        decision = row.get('decision', '').strip()
        if not decision:
            if any(row.get(k, '').strip() for k in FIELDS[8:]):
                raise ValueError(f'Review fields filled without a decision: {ident}')
            continue
        if decision not in DECISIONS:
            raise ValueError(f'Unknown decision: {decision}')
        required = ['reviewer', 'notes']
        if decision == 'supported':
            required += ['support_quote', 'source_locator', 'biological_context']
        if any(not row.get(k, '').strip() for k in required):
            raise ValueError(f'Missing review justification: {ident}')
        reviews.append({k: row.get(k, '').strip() for k in FIELDS})
    if seen != set(expected):
        raise ValueError('Missing review rows; leave undecided rows blank rather than deleting them')
    summaries = []
    for edge in graph['edges']:
        if edge['status'] != 'automatically_extracted':
            continue
        relevant = [r for r in reviews if r['edge_id'] == edge['id']]
        counts = {k: sum(r['decision'] == k for r in relevant) for k in sorted(DECISIONS)}
        summaries.append({'edge_id': edge['id'], **counts,
                          'unreviewed': len(edge['evidence_ids']) - len(relevant),
                          'mixed_review_outcomes': bool(counts['supported'] and counts['unsupported'])})
    output = new_directory(output)
    report = {'valid': True, 'graph_sha256': sha, 'decisions_sha256': file_hash(decisions_path),
              'imported_at_utc': utc_now(), 'total_evidence': len(expected),
              'reviewed_evidence': len(reviews), 'unreviewed_evidence': len(expected)-len(reviews),
              'training_ready': False,
              'note': 'Human extraction review only. No biological negative labels or automatic promotion to documented knowledge.'}
    write_json(output / 'reviews.json', {'manifest': report, 'reviews': reviews})
    write_table(output / 'edge_review_summary.csv', summaries,
                ['edge_id', 'supported', 'unsupported', 'uncertain', 'unreviewed', 'mixed_review_outcomes'])
    write_json(output / 'review_report.json', report)
    return report
