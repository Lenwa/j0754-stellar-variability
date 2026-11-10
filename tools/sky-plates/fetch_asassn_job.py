"""Read an already-submitted public Sky Patrol V1 job; never submit a new job."""
import argparse
import datetime
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import time
from urllib.parse import urljoin

import requests


def bounded_minutes(value):
    minutes = float(value)
    if not math.isfinite(minutes) or not 0 <= minutes <= 60:
        raise argparse.ArgumentTypeError("wait-minutes must be finite and between 0 and 60")
    return minutes


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.csv_links = []

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == "a" and "href" in values:
            href = values["href"]
            if any(s in href.lower() for s in ("csv", "json", ".dat", "download")):
                print("Result link:", href, flush=True)
            if href.endswith(".csv"):
                self.csv_links.append(href)


def run(job_id, out, wait_minutes):
    wait_minutes = bounded_minutes(wait_minutes)
    out.mkdir(parents=True, exist_ok=True)
    base = "https://asas-sn.osu.edu"
    job = f"{base}/sky-patrol/computation/{job_id}"
    deadline = time.monotonic()+wait_minutes*60
    while True:
        response = requests.get(job+"/status", timeout=35)
        response.raise_for_status()
        state = response.json()
        stamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%dT%H%M%SZ")
        (out / f"asassn_v1_status_{stamp}.json").write_bytes(response.content)
        print(stamp, json.dumps(state), flush=True)
        if state["status"] == "complete":
            result_url = f"{base}/sky-patrol/coordinate/{state['result']}"
            receipt = out / "asassn_v1_result_request.json"
            if receipt.exists() and json.loads(receipt.read_text())["job_url"] != job:
                raise ValueError("This output contains a different job; use a fresh directory")
            page = out / "asassn_v1_result.html"
            if not page.exists():
                response = requests.get(result_url, timeout=35)
                response.raise_for_status()
                page.write_bytes(response.content)
                receipt.write_text(json.dumps(dict(job_url=job,result_url=response.url,retrieved_utc=stamp),indent=2)+"\n")
            print("Result URL:", result_url, flush=True)
            links = Links()
            links.feed(page.read_text(encoding="utf-8"))
            if len(links.csv_links) != 1:
                raise ValueError("Expected exactly one CSV export link on the saved result page")
            csv_url = urljoin(result_url,links.csv_links[0])
            csv_file = out / "asassn_v1_lightcurve.csv"
            if not csv_file.exists():
                response = requests.get(csv_url,timeout=45)
                response.raise_for_status()
                if not response.content.startswith(b"HJD,UT Date,Camera,"):
                    raise ValueError("CSV export did not return the expected native table")
                csv_file.write_bytes(response.content)
                (out / "asassn_v1_csv_request.json").write_text(json.dumps(dict(url=response.url,status=response.status_code,
                    bytes=len(response.content),retrieved_utc=stamp),indent=2)+"\n")
            print("Native CSV bytes:",csv_file.stat().st_size,flush=True)
            return
        if state["status"] == "error":
            raise RuntimeError(state)
        if time.monotonic() >= deadline:
            print("Still pending; resume the same job URL. No new job was submitted.", flush=True)
            return
        time.sleep(min(60, max(0, deadline-time.monotonic())))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--job-id", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--wait-minutes", type=bounded_minutes, default=0)
    a = p.parse_args()
    run(a.job_id, a.out, a.wait_minutes)
