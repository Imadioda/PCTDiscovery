import argparse
import json
from pathlib import Path
from .common import read_json
from .graph import validate_graph
from .pipeline import collect, build
from .pubtator import CollectionError
from .review import export_review, import_review


def main():
    parser = argparse.ArgumentParser(description='PTCDiscovery knowledge graph — phase 1')
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ['collect', 'run']:
        sub = commands.add_parser(name)
        sub.add_argument('--rag-data', type=Path)
        sub.add_argument('--comparison', default='LUAD_lung')
        sub.add_argument('--pmids', nargs='+')
        sub.add_argument('--query', help='Optional PubTator search enrichment')
        sub.add_argument('--pages', type=int, default=1)
        sub.add_argument('--limit', type=int, default=10)
        sub.add_argument('--collection', type=Path, required=True)
        sub.add_argument('--offline', action='store_true')
        if name == 'run':
            sub.add_argument('--output', type=Path, required=True)
            sub.add_argument('--allow-partial', action='store_true')
    sub = commands.add_parser('build')
    sub.add_argument('--collection', type=Path, required=True)
    sub.add_argument('--output', type=Path, required=True)
    sub.add_argument('--rag-data', type=Path)
    sub.add_argument('--allow-partial', action='store_true')
    sub = commands.add_parser('validate')
    sub.add_argument('graph', type=Path)
    for name in ['review-export', 'review-import']:
        sub = commands.add_parser(name)
        sub.add_argument('--graph', type=Path, required=True)
        sub.add_argument('--output', type=Path, required=True)
        if name == 'review-import':
            sub.add_argument('--decisions', type=Path, required=True)
    args = parser.parse_args()
    if getattr(args, 'comparison', None) == 'all':
        args.comparison = None
    try:
        if args.command in {'collect', 'run'}:
            result = collect(args.collection, args.rag_data, args.comparison,
                             args.pmids, args.query, args.pages, args.limit, args.offline)
            print(json.dumps({k: result[k] for k in ['collected', 'failed']}, indent=2))
            if args.command == 'collect':
                return 1 if result['failed'] else 0
        if args.command in {'build', 'run'}:
            result = build(args.collection, args.output, args.rag_data, args.allow_partial)
        elif args.command == 'review-export':
            result = export_review(args.graph, args.output)
        elif args.command == 'review-import':
            result = import_review(args.graph, args.decisions, args.output)
        elif args.command == 'validate':
            result = validate_graph(read_json(args.graph))
        print(json.dumps(result, indent=2))
        return 0 if result['valid'] else 1
    except (ValueError, OSError, KeyError, CollectionError) as exc:
        parser.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    raise SystemExit(main())
