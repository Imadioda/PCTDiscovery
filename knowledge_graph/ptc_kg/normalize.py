"""ID-based normalization only: never merge nodes merely by their names."""
import re
from dataclasses import dataclass
from .common import read_csv


@dataclass(frozen=True)
class Entity:
    id: str
    type: str
    label: str
    namespace: str


def normalize(info, mention=''):
    if info.get('valid') is False or str(info.get('valid')).lower() == 'false':
        return None, 'invalid_annotation'
    kind = str(info.get('type', '')).casefold().replace('_', '')
    specs = {
        'gene': ('Gene', 'NCBIGene', r'\d+', ('NCBIGene:', 'Gene:')),
        'protein': ('Gene', 'NCBIGene', r'\d+', ('NCBIGene:', 'Gene:')),
        'disease': ('Disease', 'MESH', r'[DC]\d+', ('MESH:', 'MeSH:')),
        'chemical': ('Chemical', 'MESH', r'[DC]\d+', ('MESH:', 'MeSH:')),
        'species': ('Species', 'NCBITaxon', r'\d+', ('NCBITaxon:', 'Taxonomy:')),
        'cellline': ('CellLine', 'Cellosaurus', r'CVCL_[A-Za-z0-9]+', ('Cellosaurus:', 'CVCL:')),
        'mutation': ('Variant', 'dbSNP', r'rs\d+', ('dbSNP:',)),
        'variant': ('Variant', 'dbSNP', r'rs\d+', ('dbSNP:',)),
    }
    if kind not in specs:
        return None, 'unsupported_entity_type'
    node_type, namespace, pattern, prefixes = specs[kind]
    raw = info.get('identifier', info.get('Identifier', info.get('normalized_id', '')))
    if raw is None or str(raw).strip() in {'', '-', '-1'}:
        return None, 'missing_identifier'
    # Multiple candidate IDs are not expanded into asserted identities.
    if isinstance(raw, list) or any(s in str(raw) for s in (';', '|', ',')):
        return None, 'ambiguous_identifier'
    value = str(raw).strip()
    for prefix in prefixes:
        if value.startswith(prefix):
            value = value[len(prefix):]
            break
    if namespace == 'Cellosaurus' and not value.startswith('CVCL_'):
        value = 'CVCL_' + value
    if not re.fullmatch(pattern, value):
        return None, 'unsupported_identifier'
    normalized = info.get('normalized')
    if isinstance(normalized, list) and len({str(x) for x in normalized}) > 1:
        return None, 'ambiguous_identifier'
    if isinstance(normalized, list) and len(normalized) == 1:
        alternate = str(normalized[0])
        for prefix in prefixes:
            if alternate.startswith(prefix):
                alternate = alternate[len(prefix):]
                break
        if namespace == 'Cellosaurus' and not alternate.startswith('CVCL_'):
            alternate = 'CVCL_' + alternate
        if namespace in {'NCBIGene', 'NCBITaxon'} and alternate.isdigit():
            alternate = str(int(alternate))
            value = str(int(value))
        if alternate != value:
            return None, 'conflicting_normalized_identifier'
    if namespace in {'NCBIGene', 'NCBITaxon'} and int(value) < 1:
        return None, 'invalid_identifier'
    if namespace in {'NCBIGene', 'NCBITaxon'}:
        value = str(int(value))
    return Entity(namespace + ':' + value, node_type,
                  str(info.get('name') or mention or value), namespace), None


def hgnc_crosswalk(path):
    if not path:
        return {}, {}
    rows = read_csv(path, delimiter='\t')
    by_hgnc, by_entrez = {}, {}
    for row in rows:
        if row.get('status') != 'Approved':
            continue
        hgnc, entrez = row.get('hgnc_id', ''), row.get('entrez_id', '')
        if not re.fullmatch(r'HGNC:\d+', hgnc) or not re.fullmatch(r'\d+', entrez):
            continue
        by_hgnc.setdefault(hgnc, set()).add('NCBIGene:' + str(int(entrez)))
        by_entrez.setdefault('NCBIGene:' + str(int(entrez)), set()).add(hgnc)
    # Disallow both forward ambiguity and many-to-one crosswalks.
    accepted = {h: next(iter(ids)) for h, ids in by_hgnc.items() if len(ids) == 1
                and len(by_entrez[next(iter(ids))]) == 1}
    reverse = {value: key for key, value in accepted.items()}
    return accepted, reverse
