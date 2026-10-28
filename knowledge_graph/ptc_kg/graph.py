"""Evidence-aware graph builder. Article context never becomes disease causation."""
from collections import Counter
import math
import re
from .common import boolean, digest, read_csv
from .normalize import Entity, normalize, hgnc_crosswalk

SCHEMA_VERSION = '0.1.0'
STATUSES = {'automatically_extracted', 'dataset_observation', 'documented', 'predicted'}
SYMMETRIC = {'association', 'positive_correlation', 'negative_correlation',
             'bind', 'binding', 'interact', 'interaction', 'comparison',
             'cotreatment', 'drug_interaction'}


class KnowledgeGraph:
    def __init__(self):
        self.nodes, self.edges, self.evidence = {}, {}, {}
        self.mentions, self.issues = [], []

    def issue(self, reason, **context):
        self.issues.append({'reason': reason, **context})

    def node(self, entity):
        if entity.id in self.nodes and self.nodes[entity.id]['type'] != entity.type:
            self.issue('identity_type_conflict', entity_id=entity.id,
                       incoming_type=entity.type)
            return False
        node = self.nodes.setdefault(entity.id, {'id': entity.id, 'type': entity.type,
                                                 'namespace': entity.namespace,
                                                 'labels': [], 'xrefs': [], 'attributes': {}})
        if entity.label and entity.label not in node['labels']:
            node['labels'].append(entity.label)
            node['labels'].sort()
        return True

    def assertion(self, subject, predicate, obj, status, evidence, directed=None):
        if subject not in self.nodes or obj not in self.nodes:
            raise ValueError('Dangling assertion')
        if status not in STATUSES:
            raise ValueError('Unknown assertion status')
        if directed is False:
            subject, obj = sorted([subject, obj])
        edge_id = 'edge:' + digest([subject, predicate, obj, status, directed])[:24]
        edge = self.edges.setdefault(edge_id, {
            'id': edge_id, 'subject': subject, 'predicate': predicate,
            'object': obj, 'status': status, 'directed': directed,
            'evidence_ids': [],
        })
        evidence_id = 'evidence:' + digest([edge_id, evidence])[:24]
        self.evidence[evidence_id] = {'id': evidence_id, 'edge_id': edge_id, **evidence}
        if evidence_id not in edge['evidence_ids']:
            edge['evidence_ids'].append(evidence_id)
            edge['evidence_ids'].sort()

    def add_document(self, doc, provenance, corpus_context=None):
        pmid = str(doc.get('pmid') or doc.get('id', ''))
        if not pmid.isdigit():
            raise ValueError('Document needs a PMID')
        paper_id = 'PMID:' + pmid
        title = next((p.get('text', '') for p in doc['passages']
                      if p.get('infons', {}).get('type') == 'title'), pmid)
        self.node(Entity(paper_id, 'Publication', title, 'PMID'))
        annotations = {}
        for passage_index, passage in enumerate(doc['passages']):
            for annotation in passage.get('annotations', []):
                annotation_id = str(annotation.get('id', ''))
                entity, error = normalize(annotation.get('infons', {}), annotation.get('text', ''))
                if error:
                    self.issue(error, pmid=pmid, annotation_id=annotation_id,
                               infons=annotation.get('infons', {}))
                    continue
                if not self.node(entity):
                    continue
                locations = annotation.get('locations', [])
                spans_ok = bool(locations)
                for loc in locations:
                    start = loc.get('offset', -1) - passage.get('offset', 0)
                    length = loc.get('length', 0)
                    text = passage.get('text', '')
                    if start < 0 or length <= 0 or start + length > len(text):
                        spans_ok = False
                    elif len(locations) == 1 and text[start:start + length] != annotation.get('text', ''):
                        spans_ok = False
                if not spans_ok:
                    self.issue('mention_span_mismatch', pmid=pmid, annotation_id=annotation_id)
                annotation_record = {
                    'pmid': pmid, 'annotation_id': annotation_id,
                    'entity_id': entity.id, 'text': annotation.get('text', ''),
                    'locations': locations, 'passage_index': passage_index,
                    'span_verified': spans_ok, 'raw_infons': annotation.get('infons', {}),
                }
                if annotation_id in annotations:
                    self.issue('duplicate_annotation_id', pmid=pmid, annotation_id=annotation_id)
                annotations[annotation_id] = (entity.id, passage_index, spans_ok)
                self.mentions.append(annotation_record)
        # De-duplicate exact same relation exported at several BioC levels.
        seen = set()
        groups = [('document', doc.get('relations', []))]
        for i, passage in enumerate(doc['passages']):
            groups.append((f'passage:{i}', passage.get('relations', [])))
            for j, sentence in enumerate(passage.get('sentences', [])):
                groups.append((f'sentence:{i}:{j}', sentence.get('relations', [])))
        for scope, relations in groups:
            for relation in relations:
                key = digest(relation)
                if key in seen:
                    continue
                seen.add(key)
                info = relation.get('infons', {})
                endpoints = []
                if isinstance(info.get('role1'), dict) and isinstance(info.get('role2'), dict):
                    for role in ('role1', 'role2'):
                        entity, error = normalize(info[role])
                        if error or not self.node(entity):
                            self.issue('unresolved_relation_endpoint', pmid=pmid,
                                       relation_id=relation.get('id'), detail=error)
                            break
                        endpoints.append(entity.id)
                else:
                    # Conventional BioC: nodes refer to annotation identifiers.
                    for node in relation.get('nodes', []):
                        annotation = annotations.get(str(node.get('refid', '')))
                        if annotation:
                            endpoints.append(annotation[0])
                    if len(endpoints) != 2:
                        self.issue('unresolved_relation_endpoint', pmid=pmid,
                                   relation_id=relation.get('id'), detail='Unsupported BioC endpoint layout')
                if len(endpoints) != 2:
                    continue
                raw_type = str(info.get('type', '')).strip()
                if not raw_type:
                    self.issue('missing_relation_type', pmid=pmid, relation_id=relation.get('id'))
                    continue
                predicate = 'pubtator:' + re.sub(r'[^a-z0-9]+', '_', raw_type.lower()).strip('_')
                score = None
                if info.get('score') not in (None, ''):
                    try:
                        score = float(info['score'])
                        if not math.isfinite(score) or not 0 <= score <= 1:
                            raise ValueError()
                    except (ValueError, TypeError):
                        self.issue('invalid_extraction_score', pmid=pmid, relation_id=relation.get('id'))
                        score = None
                # PubTator3 encodes supporting mention IDs in nodes[].role;
                # its nodes[].refid can be the relation ordinal, not a mention.
                support_ids = set()
                for node in relation.get('nodes', []):
                    role = str(node.get('role', ''))
                    support_ids.update(x for x in role.split(',') if x in annotations)
                    refid = str(node.get('refid', ''))
                    if refid in annotations:
                        support_ids.add(refid)
                passage_ids = sorted({annotations[x][1] for x in support_ids})
                context_passages = [
                    {'passage_index': i, 'offset': doc['passages'][i].get('offset', 0),
                     'text': doc['passages'][i].get('text', '')} for i in passage_ids]
                evidence = {
                    'publication_id': paper_id, 'source': 'PubTator3/BioREx',
                    'source_relation_id': relation.get('id', ''),
                    'extraction_scope': scope, 'support_level': 'document',
                    'review_status': 'unreviewed', 'extraction_score': score,
                    'score_interpretation': 'Uncalibrated extraction-model score; not biological confidence',
                    'original_relation': relation, 'role_order': endpoints,
                    'context_passages': context_passages,
                    'supporting_annotation_ids': sorted(support_ids),
                    'corpus_context': corpus_context or [],
                    'biological_context_status': 'not_resolved',
                    'provenance': provenance,
                }
                # No causal direction inferred from role1/role2 or correlation.
                directed = False if predicate.removeprefix('pubtator:') in SYMMETRIC else None
                self.assertion(*[endpoints[0], predicate, endpoints[1]],
                               'automatically_extracted', evidence, directed=directed)

    def add_expression(self, results_file, comparison=None, hgnc_file=None):
        forward, reverse = hgnc_crosswalk(hgnc_file)
        for node_id, hgnc_id in reverse.items():
            if node_id in self.nodes:
                self.nodes[node_id]['xrefs'] = [hgnc_id]
                self.nodes[node_id]['attributes']['taxon_id'] = 'NCBITaxon:9606'
        seen = set()
        for row_index, row in enumerate(read_csv(results_file), 2):
            required = {'gene', 'comparison_id', 'hgnc_id', 'mapping_status',
                        'needs_review', 'selected', 'effect_log2_scale', 'fdr_bh_treat'}
            if not required.issubset(row):
                raise ValueError('Missing expression-table columns: ' + str(sorted(required - row.keys())))
            if comparison and row['comparison_id'] != comparison:
                continue
            pair = (row['comparison_id'], row['gene'])
            if pair in seen:
                raise ValueError(f'Duplicate expression row: {pair}')
            seen.add(pair)
            if not boolean(row['selected']):
                continue
            if row['mapping_status'] != 'mapped' or boolean(row['needs_review']):
                self.issue('expression_mapping_requires_review', comparison_id=row['comparison_id'], gene=row['gene'])
                continue
            hgnc = row['hgnc_id']
            if not re.fullmatch(r'HGNC:\d+', hgnc):
                raise ValueError('Invalid HGNC identifier')
            effect, fdr = float(row['effect_log2_scale']), float(row['fdr_bh_treat'])
            if not math.isfinite(effect) or not math.isfinite(fdr) or not 0 <= fdr <= 1:
                raise ValueError('Non-finite expression effect or invalid FDR')
            gene_id = forward.get(hgnc, hgnc)
            comparison_id = 'PTCComparison:' + row['comparison_id']
            self.node(Entity(gene_id, 'Gene', row.get('approved_symbol') or row['gene'], gene_id.split(':')[0]))
            if gene_id in forward.values():
                self.nodes[gene_id]['xrefs'] = [hgnc]
                self.nodes[gene_id]['attributes']['taxon_id'] = 'NCBITaxon:9606'
            else:
                self.issue('hgnc_entrez_crosswalk_unavailable', gene=row['gene'], hgnc_id=hgnc)
            self.node(Entity(comparison_id, 'StudyComparison', row.get('cancer') or row['comparison_id'], 'PTCComparison'))
            self.assertion(gene_id, 'ptc:differentially_expressed_in', comparison_id,
                           'dataset_observation', {
                               'source': 'PTCDiscovery/limma', 'review_status': 'computed_observation',
                               'original_gene_label': row['gene'], 'hgnc_id': hgnc,
                               'effect_log2_scale': effect, 'fdr_bh_treat': fdr,
                               'expression_scale': 'log2(normalized_value + 1)',
                               'source_row': row_index,
                               'cohort_details': {k: row[k] for k in ['n_normal', 'n_tumor', 'n_paired_patients', 'test', 'source_file'] if k in row},
                           }, directed=True)

    def payload(self):
        return {'schema_version': SCHEMA_VERSION,
                'nodes': sorted(self.nodes.values(), key=lambda x: x['id']),
                'edges': sorted(self.edges.values(), key=lambda x: x['id']),
                'evidence': sorted(self.evidence.values(), key=lambda x: x['id']),
                'mentions': sorted(self.mentions, key=lambda x: (x['pmid'], x['annotation_id'])),
                'issues': sorted(self.issues, key=lambda x: digest(x))}


def validate_graph(data):
    errors = []
    def indexed(name):
        values = data.get(name, [])
        ids = [x['id'] for x in values]
        if len(ids) != len(set(ids)):
            errors.append('duplicate_' + name + '_ids')
        return {x['id']: x for x in values}
    nodes, edges, evidence = indexed('nodes'), indexed('edges'), indexed('evidence')
    for edge in edges.values():
        if edge['subject'] not in nodes or edge['object'] not in nodes:
            errors.append('dangling_edge:' + edge['id'])
        if edge['status'] not in STATUSES:
            errors.append('unknown_status:' + edge['id'])
        if not edge['evidence_ids']:
            errors.append('missing_evidence:' + edge['id'])
        for evidence_id in edge['evidence_ids']:
            if evidence_id not in evidence or evidence[evidence_id]['edge_id'] != edge['id']:
                errors.append('broken_evidence_link:' + evidence_id)
    for item in evidence.values():
        if item['edge_id'] not in edges or item['id'] not in edges[item['edge_id']]['evidence_ids']:
            errors.append('orphan_evidence:' + item['id'])
        if item.get('publication_id') and item['publication_id'] not in nodes:
            errors.append('missing_publication:' + item['id'])
        if edges.get(item['edge_id'], {}).get('status') == 'automatically_extracted':
            if not item.get('publication_id') or not item.get('provenance', {}).get('payload_sha256'):
                errors.append('missing_extraction_provenance:' + item['id'])
    return {'valid': not errors, 'errors': sorted(set(errors)),
            'counts': {k: len(data.get(k, [])) for k in ['nodes', 'edges', 'evidence', 'mentions', 'issues']},
            'node_types': dict(sorted(Counter(n['type'] for n in nodes.values()).items())),
            'edge_statuses': dict(sorted(Counter(e['status'] for e in edges.values()).items())),
            'predicates': dict(sorted(Counter(e['predicate'] for e in edges.values()).items())),
            'issue_reasons': dict(sorted(Counter(i['reason'] for i in data.get('issues', [])).items()))}
