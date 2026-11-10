"""Download bounded row bands, retaining only calibrated-WCS native scan crops."""
import argparse
import io
import json
import math
import pathlib
import requests
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.wcs import WCS
from plate_pixels import decode_band, crop_origin


def byte_range(session, url, first, last):
    with session.get(url, headers={"Range":f"bytes={first}-{last}",
                                  "Accept-Encoding":"identity"}, stream=True, timeout=55) as r:
        r.raise_for_status()
        if r.status_code != 206 or r.headers.get("Content-Range", "").split("/")[0] != f"bytes {first}-{last}":
            raise ValueError("Server did not honour bounded byte range; full download refused")
        raw = r.content
        if len(raw) != last-first+1:
            raise ValueError("Incomplete range response")
        return raw, {k:v for k,v in r.headers.items() if k.lower() in
                     ["content-range", "content-length", "etag", "last-modified"]}


def run(old, out, size, orientation="x", maximum_band_mb=90, edge_clamp=False):
    out.mkdir(parents=True, exist_ok=True)
    scans = parse_single_table(old/"applause_scans.xml").to_table()
    # Orthogonal scans are sensitivity checks of one physical exposure.
    selected = [r for r in scans if int(r["plate_id"]) in [413,446,4191,4413,4599]
                and str(r["filename_scan"]).endswith(f"_{orientation}.fits")]
    session = requests.Session()
    receipts = []
    for row in selected:
        plate = int(row["plate_id"])
        output = out/f"plate_{plate}_{orientation}.fits"
        receipt_path = out/f"plate_{plate}_{orientation}_request.json"
        legacy_receipt = out/f"plate_{plate}_request.json"
        if orientation == "x" and not receipt_path.exists() and legacy_receipt.exists():
            receipt_path = legacy_receipt
        if output.exists() and receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt["crop"][2:] != [size,size]:
                raise ValueError("Cached crop has a different size; use a separate output directory")
            receipts.append(receipt)
            print(f"Cached plate {plate}", flush=True)
            continue
        url = "https://www.plate-archive.org/files/"+str(row["filename_scan"])
        prefix, metadata = byte_range(session, url, 0, 57599)
        stream = io.BytesIO(prefix)
        header = fits.Header.fromfile(stream)
        offset = stream.tell()
        (out/f"plate_{plate}_{orientation}_original_header.txt").write_text(repr(header), encoding="utf-8")
        x, y = map(float, WCS(header).world_to_pixel_values(118.69125, 16.694722222))
        nx, ny = header["NAXIS1"], header["NAXIS2"]
        x0, y0 = crop_origin(x,y,nx,ny,size,clamp=edge_clamp)
        first, last = offset+2*nx*y0, offset+2*nx*(y0+size)-1
        if last-first+1 > maximum_band_mb*1_000_000:
            raise ValueError("Row-band transfer exceeds configured cap")
        raw, band_metadata = byte_range(session, url, first, last)
        crop, h = decode_band(raw, header, x0, y0, size, size)
        fits.PrimaryHDU(crop.astype(np.float32), header=h).writeto(output, checksum=True)
        receipt = dict(plate_id=plate, scan_id=int(row["scan_id"]), url=url,
                       first_byte=first, last_byte=last, transferred_bytes=len(prefix)+len(raw),
                       original_header_bytes=offset, original_shape=[ny,nx],
                       original_target_pixel=[float(x),float(y)], crop=[x0,y0,size,size],
                       prefix_metadata=metadata, band_metadata=band_metadata)
        receipt_path.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
        receipts.append(receipt)
        print(f"Saved plate {plate}: {size} px crop from {len(raw)/1e6:.1f} MB row band", flush=True)
    name = "scan_downloads.json" if orientation == "x" else "scan_downloads_y.json"
    (out/name).write_text(json.dumps(receipts,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--old",type=pathlib.Path,required=True)
    p.add_argument("--out",type=pathlib.Path,required=True)
    p.add_argument("--size",type=int,default=1300)
    p.add_argument("--orientation",choices=["x","y"],default="x")
    p.add_argument("--maximum-band-mb",type=float,default=90)
    p.add_argument("--edge-clamp",action="store_true")
    a=p.parse_args()
    run(a.old,a.out,a.size,a.orientation,a.maximum_band_mb,a.edge_clamp)
