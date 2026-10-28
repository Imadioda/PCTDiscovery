import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def digest(value):
    if not isinstance(value, bytes):
        value = json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()
    return hashlib.sha256(value).hexdigest()


def file_hash(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                     sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def read_csv(path, delimiter=','):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def boolean(value):
    if str(value).strip().lower() not in {'true', 'false'}:
        raise ValueError(f'Expected true/false, got {value!r}')
    return str(value).strip().lower() == 'true'


def write_table(path, rows, fields):
    """Nested values use JSON, not Python repr, inside CSV cells."""
    fields = list(fields) + sorted({k for row in rows for k in row} - set(fields))
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, sort_keys=True)
                             if isinstance(v, (dict, list)) else v
                             for k, v in row.items() if k in fields})
