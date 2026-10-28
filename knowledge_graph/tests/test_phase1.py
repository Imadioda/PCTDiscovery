import copy
import io
import json
import tempfile
import unittest
from pathlib import Path
from urllib.error import HTTPError

from ptc_kg.common import digest, write_json
from ptc_kg.graph import KnowledgeGraph, validate_graph
from ptc_kg.normalize import normalize, hgnc_crosswalk
from ptc_kg.pubtator import PubTatorClient, CollectionError, documents
from ptc_kg.pipeline import select_publications, collect, build


def synthetic_document(pmid='1'):
    # Synthetic parser fixture, not biomedical evidence.
    return {'id': pmid, 'passages': [{
        'offset': 0, 'text': 'EGFR cancer', 'infons': {'type': 'abstract'},
        'annotations': [
            {'id': '1', 'text': 'EGFR', 'locations': [{'offset': 0, 'length': 4}],
             'infons': {'type': 'Gene', 'identifier': '1956', 'name': 'EGFR'}},
            {'id': '2', 'text': 'cancer', 'locations': [{'offset': 5, 'length': 6}],
             'infons': {'type': 'Disease', 'identifier': 'MESH:D009369'}},
        ], 'relations': []}], 'relations': [{
            'id': 'R1', 'infons': {'type': 'Association', 'score': '0.8',
                'role1': {'type': 'Gene', 'identifier': '1956', 'name': 'EGFR'},
                'role2': {'type': 'Disease', 'identifier': 'MESH:D009369'}},
            'nodes': [{'refid': '0', 'role': '1,2'}]}]}


def provenance():
    return {'url': 'https://example.org/fixture', 'payload_sha256': 'test-fixture-only'}


class NormalizationTests(unittest.TestCase):
    def test_namespace_not_label(self):
        a, _ = normalize({'type': 'Gene', 'identifier': '1956', 'name': 'EGFR'})
        b, _ = normalize({'type': 'Gene', 'identifier': '13649', 'name': 'EGFR'})
        self.assertNotEqual(a.id, b.id)  # same label != same species/gene

    def test_ambiguous_ids_not_expanded(self):
        for identifier in ['1;2', '1|2', ['1', '2'], '-1']:
            self.assertIsNone(normalize({'type': 'Gene', 'identifier': identifier})[0])

    def test_invalid_annotation_excluded(self):
        self.assertEqual(normalize({'type': 'Gene', 'identifier': '1', 'valid': False})[1], 'invalid_annotation')

    def test_gene_not_fabricated_protein(self):
        entity, _ = normalize({'type': 'Gene', 'identifier': '1956'})
        self.assertEqual((entity.id, entity.type), ('NCBIGene:1956', 'Gene'))

    def test_unknown_variant_preserved_in_issues(self):
        self.assertIsNone(normalize({'type': 'Variant', 'identifier': 'p.V600E'})[0])
        self.assertEqual(normalize({'type': 'Variant', 'identifier': 'rs123'})[0].id, 'dbSNP:rs123')

    def test_native_cellosaurus_id(self):
        node, error = normalize({'type': 'CellLine', 'identifier': 'CVCL:1561', 'normalized': ['1561']})
        self.assertIsNone(error)
        self.assertEqual(node.id, 'Cellosaurus:CVCL_1561')

    def test_conflicting_normalized_id(self):
        node, error = normalize({'type': 'Gene', 'identifier': '1', 'normalized': [2]})
        self.assertIsNone(node)
        self.assertEqual(error, 'conflicting_normalized_identifier')

    def test_crosswalk_ambiguity(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'hgnc.tsv'
            p.write_text('hgnc_id\tentrez_id\tstatus\nHGNC:1\t10\tApproved\nHGNC:2\t10\tApproved\nHGNC:3\t30\tApproved\n')
            self.assertEqual(hgnc_crosswalk(p)[0], {'HGNC:3': 'NCBIGene:30'})


class GraphTests(unittest.TestCase):
    def graph(self, doc=None):
        graph = KnowledgeGraph()
        graph.add_document(doc or synthetic_document(), provenance(), [{'comparison_id': 'LUAD_lung'}])
        return graph

    def test_native_role_layout_and_support(self):
        graph = self.graph()
        data = graph.payload()
        self.assertTrue(validate_graph(data)['valid'])
        self.assertEqual(len(data['edges']), 1)
        self.assertEqual(data['evidence'][0]['supporting_annotation_ids'], ['1', '2'])
        self.assertEqual(data['evidence'][0]['review_status'], 'unreviewed')
        self.assertFalse(data['edges'][0]['directed'])
        self.assertTrue(all(m['span_verified'] for m in data['mentions']))
        self.assertNotIn('PTCComparison:LUAD_lung', graph.nodes)

    def test_conventional_bioc_layout(self):
        doc = synthetic_document()
        doc['relations'][0] = {'id': 'r', 'infons': {'type': 'Association'},
                               'nodes': [{'refid': '1', 'role': 'arg1'}, {'refid': '2', 'role': 'arg2'}]}
        self.assertEqual(len(self.graph(doc).edges), 1)

    def test_comention_not_interaction(self):
        doc = synthetic_document()
        doc['relations'] = []
        self.assertEqual(len(self.graph(doc).edges), 0)

    def test_relation_dedup_within_and_across_papers(self):
        doc = synthetic_document()
        doc['passages'][0]['relations'] = copy.deepcopy(doc['relations'])
        graph = self.graph(doc)
        self.assertEqual(len(graph.evidence), 1)
        graph.add_document(synthetic_document('2'), provenance())
        self.assertEqual(len(graph.edges), 1)
        self.assertEqual(len(graph.evidence), 2)

    def test_no_invented_causal_direction(self):
        doc = synthetic_document()
        doc['relations'][0]['infons']['type'] = 'Inhibition'
        self.assertIsNone(self.graph(doc).payload()['edges'][0]['directed'])

    def test_invalid_scores_not_probabilities(self):
        doc = synthetic_document()
        doc['relations'][0]['infons']['score'] = 'nan'
        data = self.graph(doc).payload()
        self.assertIsNone(data['evidence'][0]['extraction_score'])
        self.assertIn('invalid_extraction_score', [x['reason'] for x in data['issues']])

    def test_unresolved_relation_endpoint_quarantined(self):
        doc = synthetic_document()
        doc['relations'][0]['infons']['role2']['identifier'] = 'D001;D002'
        graph = self.graph(doc)
        self.assertEqual(len(graph.edges), 0)
        self.assertIn('unresolved_relation_endpoint', [x['reason'] for x in graph.issues])

    def test_span_mismatch_reported(self):
        doc = synthetic_document()
        doc['passages'][0]['annotations'][0]['locations'][0]['offset'] = 999
        data = self.graph(doc).payload()
        self.assertIn('mention_span_mismatch', [x['reason'] for x in data['issues']])

    def test_broken_evidence_and_dangling_edge_detected(self):
        data = self.graph().payload()
        data['edges'][0]['object'] = 'missing'
        data['evidence'] = []
        self.assertFalse(validate_graph(data)['valid'])

    def test_expression_is_observation_not_causation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'results.csv'
            p.write_text('gene,comparison_id,hgnc_id,mapping_status,needs_review,selected,effect_log2_scale,fdr_bh_treat\nEGFR,LUAD_lung,HGNC:3236,mapped,False,True,2.0,0.01\n')
            graph = self.graph()
            graph.add_expression(p, 'LUAD_lung')
            data = graph.payload()
            self.assertTrue(validate_graph(data)['valid'])
            self.assertIn('dataset_observation', [e['status'] for e in data['edges']])
            self.assertIn('PTCComparison:LUAD_lung', graph.nodes)
            self.assertEqual(graph.nodes['HGNC:3236']['type'], 'Gene')


class SelectionTests(unittest.TestCase):
    def test_retractions_comparison_and_specific_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); lit = root / 'literature_all_cancers'; lit.mkdir()
            write_json(lit / 'articles.json', [
                {'article_id': 'MED:1', 'publication_flags': ['Retracted Publication']},
                {'article_id': 'MED:42562391'}, {'article_id': 'MED:2'}, {'article_id': 'MED:3'}])
            (lit / 'gene_article_links.csv').write_text(
                'comparison_id,approved_symbol,article_id\n'
                'LUAD_lung,A,MED:1\nLUAD_lung,CLDN18,MED:42562391\n'
                'LUAD_lung,A,MED:2\nBRCA_breast,A,MED:3\n')
            ids, contexts, policy = select_publications(root, 'LUAD_lung')
            self.assertEqual(ids, ['2'])
            self.assertEqual(policy['excluded_flagged_article_ids'], ['MED:1'])
            with self.assertRaises(ValueError):
                select_publications(root, 'LUAD_lung', ['1'])



class ClientTests(unittest.TestCase):
    def test_wrappers(self):
        doc = synthetic_document()
        for data in [doc, [doc], {'PubTator3': [doc]}, {'documents': [doc]}]:
            self.assertEqual(documents(data), [doc])
        with self.assertRaises(ValueError):
            documents({'error': 'no documents'})

    def test_cache_offline_and_tamper(self):
        calls = []
        def opener(request, timeout):
            calls.append(request.full_url)
            return io.BytesIO(json.dumps({'PubTator3': [synthetic_document()]}).encode())
        with tempfile.TemporaryDirectory() as tmp:
            client = PubTatorClient(tmp, opener=opener)
            _, _, p = client.publication('1')
            offline = PubTatorClient(tmp, offline=True)
            offline.publication('1')
            self.assertEqual(len(calls), 1)
            record = json.loads(p.read_text())
            record['payload']['PubTator3'][0]['id'] = 'changed'
            write_json(p, record)
            with self.assertRaises(CollectionError):
                offline.publication('1')
            with self.assertRaises(CollectionError):
                offline.publication('2')

    def test_bounded_http_retry(self):
        calls, waits = [], []
        def opener(request, timeout):
            calls.append(1)
            raise HTTPError(request.full_url, 503, 'offline test', {}, None)
        with tempfile.TemporaryDirectory() as tmp:
            client = PubTatorClient(tmp, attempts=2, opener=opener, sleeper=waits.append)
            with self.assertRaises(CollectionError):
                client.publication('1')
            self.assertEqual(len(calls), 2)
            self.assertFalse(list(Path(tmp).glob('*.json')))

    def test_wrong_pmid_not_cached(self):
        def opener(request, timeout):
            return io.BytesIO(json.dumps({'PubTator3': [synthetic_document('999')]}).encode())
        with tempfile.TemporaryDirectory() as tmp:
            client = PubTatorClient(tmp, opener=opener)
            with self.assertRaises(CollectionError):
                client.publication('1')
            self.assertFalse(list(Path(tmp).glob('*.json')))

    def test_search_contract(self):
        def opener(request, timeout):
            return io.BytesIO(b'{"results":[{"pmid":1},{"pmid":"2"}]}')
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(PubTatorClient(tmp, opener=opener).search('EGFR')[0], ['1', '2'])

    def test_offline_pipeline_and_deterministic_graph(self):
        def opener(request, timeout):
            return io.BytesIO(json.dumps({'PubTator3': [synthetic_document()]}).encode())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            client = PubTatorClient(root / 'collection/cache', opener=opener)
            collect(root / 'collection', pmids=['1'], client=client)
            report = build(root / 'collection', root / 'g1')
            self.assertTrue(report['valid'])
            build(root / 'collection', root / 'g2')
            self.assertEqual((root / 'g1/graph.json').read_bytes(), (root / 'g2/graph.json').read_bytes())
            with self.assertRaises(ValueError):
                build(root / 'collection', root / 'g1')


if __name__ == '__main__':
    unittest.main()
