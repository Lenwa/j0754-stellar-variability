"""Snapshot public Sky Patrol parquet blocks using the official client's protocol."""
import argparse
import base64
import datetime
import io
import json
from pathlib import Path

import pyarrow.parquet as pq
import requests


def save_response(out, name, raw):
    """Preserve bad responses separately so they cannot poison a retry's cache."""
    if name.endswith(".parquet"):
        try:
            pq.read_table(io.BytesIO(raw))
        except (ValueError, OSError):
            rejected = out / (name+".rejected.bin")
            number = 1
            while rejected.exists():
                rejected = out / (name+f".rejected-{number}.bin")
                number += 1
            rejected.write_bytes(raw)
            raise
    (out / name).write_bytes(raw)


def run(previous, out):
    out.mkdir(parents=True, exist_ok=True)
    old = json.loads((previous / "isotonic_pixel_results.json").read_text())
    ids = [int(x) for x in old["results"][-1]["selected_controls"]]
    ids.append(667401913089980544)
    cols = ["asas_sn_id", "ra_deg", "dec_deg", "gaia_id", "tic_id"]
    payload = dict(tar_ids=ids, catalog="stellar_main", id_col="gaia_id",
                   cols=cols, format="arrow", download=True)
    base = "http://asassn-lb01.ifa.hawaii.edu:9006"
    receipt_path = out / "asassn_requests.json"
    receipt = (json.loads(receipt_path.read_text()) if receipt_path.exists() else
               dict(created_utc=datetime.datetime.now(datetime.UTC).isoformat(),
                    payload=payload, client_reference="https://github.com/asas-sn/skypatrol/blob/master/pyasassn/client.py",
                    requests=[]))

    def get_file(method, url, name, **kwargs):
        if (out / name).exists():
            return (out / name).read_bytes()
        r = requests.request(method, url, timeout=60, **kwargs)
        receipt["requests"].append(dict(method=method, url=r.url, status=r.status_code,
                                         file=name, bytes=len(r.content)))
        receipt_path.write_text(json.dumps(receipt, indent=2)+"\n")
        r.raise_for_status()
        save_response(out, name, r.content)
        return r.content

    # A fresh lookup registers the download job; a cached index alone cannot do so.
    if (out / "asassn_index.parquet").exists() and not (out / "asassn_lightcurves.parquet").exists():
        raise RuntimeError("Partial cached query: use a fresh output directory to register its download job")
    raw = get_file("POST", base+"/lookup_targets/catalog_list", "asassn_index.parquet", json=payload)
    index = pq.read_table(io.BytesIO(raw)).to_pylist()
    print("Catalogue:", json.dumps(index), flush=True)
    if not index:
        return
    servers = json.loads(get_file("GET", base+"/get_block_servers", "asassn_block_servers.json"))
    query_id = (f"listlen-{len(ids)}_listfirst-{ids[0]}_listend-{ids[-1]}"
                "|catalog-stellar_main|id_col-gaia_id|cols-"+"/".join(cols))
    query_hash = base64.encodebytes(query_id.encode()).decode()
    # There are fewer than 1000 sources, hence one native block.
    failures = []
    for host in servers:
        url = f"http://{host}:9006/get_block/query_hash-{query_hash}-block_idx-0-catalog-stellar_main"
        try:
            raw = get_file("GET", url, "asassn_lightcurves.parquet")
            table = pq.read_table(io.BytesIO(raw))
            print("Light curves:", table.num_rows, "rows; columns", table.column_names, flush=True)
            break
        except (requests.RequestException, ValueError, OSError) as e:
            failures.append(str(e))
    else:
        raise RuntimeError("All public block hosts failed: "+repr(failures))
    receipt["block_failures"] = failures
    receipt["catalogue_rows"] = len(index)
    receipt["lightcurve_rows"] = table.num_rows
    receipt_path.write_text(json.dumps(receipt, indent=2)+"\n")
    # Public variable-catalogue identity was obtained by its coordinate-search form.
    variable = "https://asas-sn.osu.edu/variables/05a882b3-3477-5b94-b81b-17a186f6b4c1"
    for suffix, name in [("", "asassn_variable_target.html"), (".json", "asassn_variable_target.json"),
                         (".csv", "asassn_variable_target.csv")]:
        get_file("GET", variable+suffix, name)
    # Verify the absence through the single-target and master-cone query routes.
    single = dict(payload, tar_ids=[667401913089980544])
    single_query = ("listlen-1_listfirst-667401913089980544_listend-667401913089980544"
                    "|catalog-stellar_main|id_col-gaia_id|cols-"+"/".join(cols))
    radius, ra, dec = .01, 118.69117410685, 16.694621231895
    cone_cols = ["asas_sn_id", "ra_deg", "dec_deg", "catalog_sources"]
    cone = dict(catalog="master_list", cols=cone_cols, format="arrow", download=True)
    cone_query = f"conectr-{radius}_conera-{ra}_conedec-{dec}|catalog-master_list|cols-"+"/".join(cone_cols)
    for stem, body, lookup, query_id, catalog, block_name in [
        ("asassn_target", single, base+"/lookup_targets/catalog_list", single_query,
         "stellar_main", "asassn_target_lightcurve.parquet"),
        ("asassn_master_cone", cone, base+f"/lookup_cone/radius{radius}_ra{ra}_dec{dec}", cone_query,
         "master_list", "asassn_master_cone_lightcurves.parquet")]:
        # A block that is not cached needs a fresh lookup even if its index exists.
        if not (out / block_name).exists() and (out / (stem+"_index.parquet")).exists():
            raise RuntimeError("Partial cached query: use a fresh output directory to register its download job")
        get_file("POST", lookup, stem+"_index.parquet", json=body)
        encoded = base64.encodebytes(query_id.encode()).decode()
        get_file("GET", f"http://{servers[0]}:9006/get_block/query_hash-{encoded}-block_idx-0-catalog-{catalog}", block_name)
    receipt["supplementary_query_definitions"] = dict(single=single, cone=cone, radius=radius, ra=ra, dec=dec)
    receipt_path.write_text(json.dumps(receipt, indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--previous", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    run(a.previous, a.out)
