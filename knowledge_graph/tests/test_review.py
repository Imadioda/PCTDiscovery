import tempfile
import unittest
from pathlib import Path
from test_phase1 import synthetic_document, provenance
from ptc_kg.graph import KnowledgeGraph
from ptc_kg.common import write_json, read_json, read_csv, write_table, file_hash
from ptc_kg.review import export_review, import_review, FIELDS


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        kg = KnowledgeGraph()
        kg.add_document(synthetic_document(), provenance(), [])
        self.graph = self.root / 'graph.json'
        write_json(self.graph, kg.payload())
        self.sha = file_hash(self.graph)
        self.export = self.root / 'export'
        export_review(self.graph, self.export)
        self.csv = self.export / 'decisions.csv'
        self.rows = read_csv(self.csv)

    def run_import(self):
        write_table(self.csv, self.rows, FIELDS)
        return import_review(self.graph, self.csv, self.root / 'imported')

    def supported(self):
        self.rows[0].update(decision='supported', reviewer='Fixture reviewer',
                            notes='Synthetic test only', support_quote='Fixture quote',
                            source_locator='Synthetic paragraph 1', biological_context='Synthetic')

    def test_blank_stays_unreviewed(self):
        result = self.run_import()
        self.assertEqual(result['unreviewed_evidence'], 1)
        self.assertFalse(result['training_ready'])

    def test_supported_preserves_original_graph(self):
        self.supported()
        result = self.run_import()
        self.assertEqual(result['reviewed_evidence'], 1)
        self.assertEqual(file_hash(self.graph), self.sha)
        self.assertEqual(read_json(self.graph)['edges'][0]['status'], 'automatically_extracted')

    def test_missing_quote_rejected(self):
        self.supported(); self.rows[0]['support_quote'] = ''
        with self.assertRaises(ValueError): self.run_import()

    def test_tampered_identity_rejected(self):
        self.rows[0]['subject'] = 'NCBIGene:999'
        with self.assertRaises(ValueError): self.run_import()

    def test_wrong_graph_rejected(self):
        self.rows[0]['graph_sha256'] = 'wrong'
        with self.assertRaises(ValueError): self.run_import()

    def test_duplicate_rejected(self):
        self.rows.append(dict(self.rows[0]))
        with self.assertRaises(ValueError): self.run_import()

    def test_missing_row_rejected(self):
        self.rows.clear()
        with self.assertRaises(ValueError): self.run_import()

    def test_invalid_decision_rejected(self):
        self.rows[0]['decision'] = 'true'
        with self.assertRaises(ValueError): self.run_import()

    def test_no_overwrite(self):
        with self.assertRaises(ValueError): export_review(self.graph, self.export)

    def test_html_escapes_source_text(self):
        g = read_json(self.graph)
        g['nodes'][0]['labels'] = ['<script>alert(1)</script>']
        write_json(self.graph, g)
        target = self.root / 'escaped'
        export_review(self.graph, target)
        page = (target / 'review.html').read_text()
        self.assertNotIn('<script>', page)
        self.assertIn('&lt;script&gt;', page)
