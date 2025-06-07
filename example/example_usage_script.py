from Quagga import stripe_caller


stripe_caller(
    hic_file='/Documents/data/hic/GM12878.hic',
    reference_genome='hg38',
    output_file='GM12878_chr7.bedpe',
    chromosomes=['chr7'],
    norm='balanced',
    threshold=0.15,
    resolution=10000,
    max_range=2000000,
    min_length=100000,
    min_distance=100000,
    window_size=5,
    centromere_file='removed_regions.bed',
    N_threads=1,
    nstrata_blank=1,
    sigma=2,
    rel_height=0.2,
)


