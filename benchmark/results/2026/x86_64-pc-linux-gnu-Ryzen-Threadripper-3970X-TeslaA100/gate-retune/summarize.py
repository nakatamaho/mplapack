#!/usr/bin/env python3
"""Rebuild summary.csv from preserved, pre-retune completed calls only."""
import csv
import json
from pathlib import Path
import statistics

base = Path(__file__).resolve().parent
groups = {}
for path in sorted((base / 'raw').glob('*-*.jsonl')):
    for line in path.read_text().splitlines():
        r = json.loads(line)
        if 'backend' not in r or r.get('returncode') != 0:
            continue
        key = tuple(r[k] for k in ('backend', 'routine', 'n', 'nrhs', 'uplo'))
        groups.setdefault(key, {}).setdefault(r['mode'], []).append(r)
fields = ['backend', 'routine', 'n', 'nrhs', 'uplo', 'complete']
fields += [m + '_' + f for m in ('cpu', 'default', 'forced')
           for f in ('rounds', 'cold', 'warm', 'warm_min', 'warm_max', 'entry_traces')]
fields += ['forced_over_cpu', 'default_over_cpu']
with (base / 'summary.csv').open('w') as out:
    writer = csv.DictWriter(out, fields, lineterminator='\n')
    writer.writeheader()
    for key, modes in sorted(groups.items()):
        expected = 2 if key[1][1:] in ('getri', 'potri', 'trtri', 'lauum') and key[2] >= 2048 else 3
        row = dict(zip(fields[:5], key))
        row['complete'] = all(len(modes.get(m, [])) >= expected for m in ('cpu', 'default', 'forced'))
        medians = {}
        for mode, records in modes.items():
            warm = [statistics.median(c['seconds'] for c in r['calls'][1:]) for r in records]
            medians[mode] = statistics.median(warm)
            row.update({mode+'_rounds': len(records),
                        mode+'_cold': statistics.median(r['calls'][0]['seconds'] for r in records),
                        mode+'_warm': medians[mode], mode+'_warm_min': min(warm),
                        mode+'_warm_max': max(warm),
                        mode+'_entry_traces': sum(r['entry_trace_count'] for r in records)})
        for mode in ('forced', 'default'):
            if mode in medians and 'cpu' in medians:
                row[mode+'_over_cpu'] = medians[mode] / medians['cpu']
        writer.writerow(row)
print(f'{len(groups)} points, {sum(len(rs) for modes in groups.values() for rs in modes.values())} completed process records')
