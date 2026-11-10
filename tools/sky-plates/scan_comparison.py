"""Record raw orthogonal-scan diagnostics without inferring stellar disappearance."""
import argparse
import json
import pathlib
import warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS,FITSFixedWarning
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plate_pixels import aperture


def run(out):
    warnings.simplefilter("ignore",FITSFixedWarning)
    gaia=Table.read(out/"gaia_field.csv",format="ascii.csv")
    ids=[667401913089980544,667024467067954560,667023509291889792]
    fig,axes=plt.subplots(3,2,figsize=(7,10),layout="constrained")
    records=[]
    for column,ori in enumerate(["x","y"]):
        name="isotonic_pixel_results.json" if ori=="x" else "isotonic_y_pixel_results.json"
        record=json.loads((out/name).read_text())["results"][-1]
        receipts=json.loads((out/("scan_downloads.json" if ori=="x" else "scan_downloads_y.json")).read_text())
        receipt=next(r for r in receipts if r["plate_id"]==4599)
        with fits.open(out/f"plate_4599_{ori}.fits") as h:image=h[0].data.astype(float);w=WCS(h[0].header)
        for row,sid in enumerate(ids):
            q=gaia[gaia["source_id"]==sid][0];year=record["year"]
            ra=float(q["ra"])+(year-2016)*float(q["pmra"])/(3.6e6*np.cos(np.radians(float(q["dec"]))))
            de=float(q["dec"])+(year-2016)*float(q["pmdec"])/3.6e6
            x,y=map(float,w.world_to_pixel_values(ra,de))
            x+=record["position_shift_pixels"][0];y+=record["position_shift_pixels"][1]
            measure=aperture(image,x,y,6,record["polarity"])
            original_x,original_y=x+receipt["crop"][0],y+receipt["crop"][1]
            ny,nx=receipt["original_shape"]
            records.append(dict(source_id=str(sid),orientation=ori,plate_id=4599,
                                original_pixel=[original_x,original_y],
                                nearest_scan_boundary_pixels=min(original_x,original_y,nx-original_x,ny-original_y),
                                signal=measure["signal"],blank_sigma=record["apertures"]["6"]["blank_sigma"],
                                calibration_gate=record["apertures"]["6"]["calibration_gate"]))
            ix,iy=round(x),round(y);patch=image[iy-60:iy+61,ix-60:ix+61]
            bg=np.median(patch);noise=1.4826*np.median(abs(patch-bg))
            ax=axes[row,column];ax.imshow(bg-patch,origin="lower",cmap="gray",vmin=-noise*3,vmax=noise*20)
            ax.add_patch(plt.Circle((60+x-ix,60+y-iy),6,fill=False,color="orange"))
            ax.set_xticks([]);ax.set_yticks([])
            ax.set_title(f"{'J0754' if row==0 else 'Comparison â€¦'+str(sid)[-7:]} | {ori.upper()} scan",fontsize=10)
    fig.suptitle("Same 1937 photographic plate, two scan orientations\nPhotographic signal shown positive; scales set separately",fontsize=12)
    fig.savefig(out/"orthogonal_scan_check.png",dpi=150);plt.close(fig)
    output=dict(measurements=records,conclusion="Both troublesome positions have weak signal in both orientations; scanner-only error is less plausible, but historical variability and shared plate/coverage effects remain unresolved.",
                limits="The scans are the same physical plate. Photometric fits have different control membership and are not an isolated scanner-response experiment. No calibrated stellar upper limit inferred.")
    (out/"scan_comparison_summary.json").write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(output,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True)
    run(p.parse_args().out)
