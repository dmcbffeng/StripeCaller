import json
import sys

payload = {}
presets = ["HiC_5000","MiC_5000"]

prekeys = ['reference_genome', 'nstrata_blank', 'chrs', 'norm', 'threshold', 'max_range', 'resolution', 'min_length', 'min_distance', 'merge', 'window_size', 'step', 'sigma', 'rel_height']

prevals = [
['hg38', 1, ['chr1','chr5','chr8'], "balanced", 0.15, 2000000, 5000, 300000, 2, 8, 35, 1, 0.3],
['hg38', 50, ['chr1','chr5','chr8'], "balanced", 0.15, 2000000, 5000, 300000, 2, 8, 600, 12, 0.3],
] 

for j in range(len(presets)):
    payload[presets[j]] = {k:v for k,v in zip(prekeys,prevals[j])}

with open("presets.json", "w") as out:
    json.dump(payload, out)
