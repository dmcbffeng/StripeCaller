import sys
sys.path.append("../StripeCaller/")
from caller.StripeCaller import stripe_caller_all # stripe_caller
from caller.buildarg import stripe_parser


if __name__ == "__main__":
    args = stripe_parser().parse_args()
    if args.presets != "":
        jfile = open('presets.json')
        data = json.load(jfile)[args.presets]
        for key in data:
            args.__dict__[key] = data[key]

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
    

#python example_usage.py --hic '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/GM12878/GM12878.hic' --output "TEST.txt" --chr "chr1" --rg 'hg38' --max_range 2000000 --resolution 5000 --min_length 300000 --min_distance 500000 --merge 2 --window_size 8 --N_threads 26 --step 35 --sigma 1 --rel_height 0.3 --norm balanced --thr 0.15
