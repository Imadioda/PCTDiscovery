"""PubTator 3 REST client, snapshots, bounded retries and offline replay."""
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .common import digest, read_json, utc_now, write_json

BASE = 'https://www.ncbi.nlm.nih.gov/research/pubtator3-api'


class CollectionError(RuntimeError):
    pass


def documents(payload):
    """Current PubTator3 wrapper and conventional BioC collection layouts."""
    if isinstance(payload, list):
        result = []
        for entry in payload:
            result.extend(documents(entry))
        return result
    if not isinstance(payload, dict):
        raise ValueError('BioC payload must be an object or list')
    if 'PubTator3' in payload:
        return documents(payload['PubTator3'])
    if 'documents' in payload:
        return documents(payload['documents'])
    if 'id' in payload and isinstance(payload.get('passages'), list):
        return [payload]
    raise ValueError('Unrecognized BioC structure')


class PubTatorClient:
    def __init__(self, cache, offline=False, interval=1.0, attempts=3,
                 timeout=45, opener=urlopen, sleeper=time.sleep):
        if interval < 0.5 or not 1 <= attempts <= 5:
            raise ValueError('Use interval >= 0.5 s and 1–5 attempts')
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline, self.interval, self.attempts = offline, interval, attempts
        self.timeout, self.opener, self.sleeper = timeout, opener, sleeper
        self.last_request = 0.0

    def get(self, endpoint, parameters):
        url = BASE + '/' + endpoint.lstrip('/') + '?' + urlencode(parameters)
        path = self.cache / (digest(url) + '.json')
        if path.exists():
            record = read_json(path)
            if record.get('url') != url or digest(record.get('payload')) != record.get('payload_sha256'):
                raise CollectionError(f'Corrupt cache: {path}')
            return record, path
        if self.offline:
            raise CollectionError(f'No cached response in offline mode: {endpoint} {parameters}')
        for attempt in range(self.attempts):
            wait = self.interval - (time.monotonic() - self.last_request)
            if wait > 0:
                self.sleeper(wait)
            self.last_request = time.monotonic()
            try:
                request = Request(url, headers={'User-Agent': 'PTCDiscovery-KG/0.1',
                                               'Accept': 'application/json'})
                with self.opener(request, timeout=self.timeout) as response:
                    raw = response.read()
                try:
                    payload = json.loads(raw)
                except (ValueError, UnicodeError) as exc:
                    raise CollectionError('Non-JSON PubTator response') from exc
                if endpoint.startswith('publications/'):
                    matches = [d for d in documents(payload) if str(d.get('pmid') or d['id']) == parameters['pmids']]
                    if len(matches) != 1:
                        raise CollectionError('Export response missing or duplicating the requested PMID')
                elif endpoint == 'search/':
                    if not isinstance(payload, dict) or not isinstance(payload.get('results'), list):
                        raise CollectionError('Unrecognized search response')
                record = {'url': url, 'retrieved_at_utc': utc_now(),
                          'payload_sha256': digest(payload),
                          'response_bytes_sha256': digest(raw), 'payload': payload}
                # Only valid successful JSON responses enter the cache.
                write_json(path, record)
                return record, path
            except HTTPError as exc:
                retryable = exc.code in {429, 500, 502, 503, 504}
                if not retryable or attempt == self.attempts - 1:
                    raise CollectionError(f'HTTP {exc.code} from PubTator') from exc
                delay = min(30, 2 ** (attempt + 1))
                try:
                    delay = min(30, max(delay, float(exc.headers.get('Retry-After', 0))))
                except (ValueError, TypeError, AttributeError):
                    pass
                self.sleeper(delay)
            except (URLError, TimeoutError, OSError) as exc:
                if attempt == self.attempts - 1:
                    raise CollectionError(f'PubTator transport error: {type(exc).__name__}') from exc
                self.sleeper(min(30, 2 ** (attempt + 1)))
        raise CollectionError('Request failed')

    def publication(self, pmid):
        if not str(pmid).isdigit():
            raise ValueError('PMID must contain digits only')
        record, path = self.get('publications/export/biocjson', {'pmids': str(pmid)})
        docs = documents(record['payload'])
        matching = [doc for doc in docs if str(doc.get('pmid') or doc['id']) == str(pmid)]
        if len(matching) != 1:
            raise CollectionError(f'Expected one BioC document for PMID {pmid}; got {len(matching)}')
        return matching[0], record, path

    def search(self, query, pages=1):
        if not query.strip() or not 1 <= pages <= 5:
            raise ValueError('Nonempty query and 1–5 pages required')
        pmids, records = set(), []
        for page in range(1, pages + 1):
            record, path = self.get('search/', {'text': query, 'page': page})
            data = record['payload']
            if not isinstance(data, dict) or not isinstance(data.get('results'), list):
                raise CollectionError('Unrecognized PubTator search response')
            records.append({'path': str(path), 'sha256': record['payload_sha256'],
                            'url': record['url']})
            for item in data['results']:
                value = str(item.get('pmid', ''))
                if not value.isdigit():
                    raise CollectionError('Search result has no valid PMID')
                pmids.add(value)
            if not data['results']:
                break
        return sorted(pmids, key=int), records
