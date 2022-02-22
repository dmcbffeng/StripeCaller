import sys
sys.path.append("../StripeCaller/")
from caller.StripeCaller import stripe_caller_all # stripe_caller
from caller.buildarg import stripe_parser

if __name__ == '__main__':
    args = stripe_parser().parse_args()
    print(vars(args))
    
    stripe_caller_all(
        hic_file=args.hic,
        reference_genome=args.reference_genome,
        chromosomes=args.chrs,
        output_file=args.output_file,
        norm=args.norm,
        threshold=args.threshold,
        max_range=args.max_range, resolution=args.resolution,
        min_length=args.min_length, min_distance=args.min_distance,
        merge=args.merge, window_size=args.window_size,
        centromere_file=args.centromere_file,
        N_threads=args.N_threads,

        nstrata_blank=args.nstrata_blank,
        step=args.step,
        sigma=args.sigma,
        rel_height=args.rel_height
    )
