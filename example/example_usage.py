from StripeCaller import stripe_caller


if __name__ == '__main__':
    stripe_caller(
        hic_file='../../GM12878.hic',
        reference_genome='hg38',
        chromosomes=['chr1'],
        output_file='GM12878_chr1.bedpe',
        norm='balanced',
        threshold=0.15,
        max_range=2000000, resolution=5000,
        min_length=300000, min_distance=500000,
        merge=2, window_size=8,
        centromere_file='removed_regions.bed',
        N_threads=1,

        nstrata_blank=0, step=400, sigma=12., rel_height=0.3
    )


