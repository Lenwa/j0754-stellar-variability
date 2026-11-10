"""Bounded, reproducible public-archive GET snapshots for J0754 history."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests


def snapshot(out, name, url, timeout=60):
    out.mkdir(parents=True, exist_ok=True)
    path = out / name
    record = {"file": name, "request_url": url,
              "retrieved_utc": datetime.now(timezone.utc).isoformat()}
    try:
        response = requests.get(url, timeout=timeout)
        path.write_bytes(response.content)
        record.update(final_url=response.url, status=response.status_code,
                      content_type=response.headers.get("Content-Type"),
                      bytes=len(response.content),
                      sha256=hashlib.sha256(response.content).hexdigest())
    except requests.RequestException as exc:
        record["error"] = str(exc)
    receipts = out / "query_receipts.jsonl"
    with receipts.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record) + "\n")
    print(json.dumps({key: value for key, value in record.items()
                      if key not in ('request_url', 'final_url')}), flush=True)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--name", required=True)
    parser.add_argument("--url", required=True)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--tap-query", help="Read-only ADQL appended to a public TAP sync URL")
    args = parser.parse_args()
    if args.tap_query:
        args.url = requests.Request('GET', args.url, params=dict(
            REQUEST='doQuery', LANG='ADQL', FORMAT='votable', QUERY=args.tap_query)).prepare().url
    snapshot(args.out, args.name, args.url, args.timeout)
