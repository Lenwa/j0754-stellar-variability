"""Bounded ASAS retrieval for the sixteen previously selected J0754 controls."""
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from fetch_historical_extension import snapshot


def designation(ra, dec):
    seconds = round(ra / 15 * 3600)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    dec_tenths = round(abs(dec) * 600)
    degrees, minute_tenths = divmod(dec_tenths, 600)
    return f'{hours:02d}{minutes:02d}{seconds:02d}{"+" if dec >= 0 else "-"}{degrees:02d}{minute_tenths/10:04.1f}'


if __name__ == '__main__':
    base = Path('results/sky-plates-j0754-followup-20261108')
    out = Path('results/sky-plates-j0754-models-20261108/history/controls')
    ids = json.loads((base / 'isotonic_pixel_results.json').read_text())['results'][-1]['selected_controls']
    with (base / 'gaia_field.csv').open(newline='') as stream:
        field = {r['source_id']: r for r in csv.DictReader(stream)}
    selection = []
    for source_id in ids:
        row = dict(field[source_id])
        row['asas_query_designation'] = designation(float(row['ra']), float(row['dec']))
        row['url'] = 'https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?' + row['asas_query_designation'] + ',asas3'
        selection.append(row)
    (out / 'selected_control_queries.json').write_text(json.dumps(selection, indent=2))
    def fetch(row):
        return snapshot(out, row['source_id'] + '.txt', row['url'], timeout=25)
    with ThreadPoolExecutor(max_workers=4) as executor:
        list(executor.map(fetch, selection))
