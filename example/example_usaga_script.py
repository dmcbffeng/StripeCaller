from StripeCaller import stripe_caller


stripe_caller(
    hic_file='../../GM12878.hic',
    reference_genome='hg38',
    output_file='GM12878_chr1_chr2.bedpe',
    chromosomes=['chr1', 'chr2'],
    norm='balanced',
    threshold=0.15,
    resolution=5000,
    max_range=200000,
    min_length=300000,
    min_distance=300000,
    merge=3,
    window_size=10,
    centromere_file='removed_regions.bed',
    N_threads=1,
    nstrata_blank=1,
    step=80,
    sigma=2,
    rel_height=0.3
)


