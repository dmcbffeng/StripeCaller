import os
import sys
sys.path.append(os.path.abspath("../../.."))
for p in sys.path:
    print(p)
# from caller.StripeCaller import stripe_caller_all as stripe_caller
from StripeCaller import stripe_caller


stripe_caller(
    # hic_file='../../GM12878.hic',
    hic_file='/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/GM12878/GM12878.hic',
    reference_genome='hg38',
    output_file='GM12878_chr2.bedpe',
    chromosomes=['chr2'],
    norm='balanced',
    threshold=0.15,
    resolution=5000,
    max_range=2000000,
    min_length=300000,
    min_distance=300000,
    merge=3,
    window_size=10,
    centromere_file=None,
    N_threads=1,
    nstrata_blank=1,
    step=80,
    sigma=2,
    rel_height=0.3
)


