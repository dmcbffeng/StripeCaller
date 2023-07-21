import os
import sys
sys.path.append(os.path.abspath("../"))
for p in sys.path:
    print(p)
# from caller.StripeCaller import stripe_caller_all as stripe_caller
from Quagga import stripe_caller


stripe_caller(
    hic_file='../../O3.hic',
    # hic_file='/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/GM12878/GM12878.hic',
    reference_genome='hg19',
    output_file='O3_chr1.bedpe',
    chromosomes=['chr1'],
    norm='balanced',
    threshold=0.15,
    resolution=25000,
    max_range=2000000,
    min_length=300000,
    min_distance=300000,
    window_size=5,
    centromere_file='removed_regions.bed',
    N_threads=1,
    nstrata_blank=1,
    sigma=2,
    rel_height=0.2,
    log_path='O3_example_log.txt'
)


