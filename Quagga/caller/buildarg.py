import os
import argparse
import textwrap
import sys


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def stripe_parser():
    parser = MyParser(
            description='Call stripes',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--hic',
        type=str,
        help='Hi-C file path')

    parser.add_argument(
        '--rg', dest="reference_genome",
        type=str,
        help='reference_genome name (e.g., hg38)')

    parser.add_argument(
        '--chrs',
        nargs="+",
        default=["chr1"],
        help='which chromosomes to calculate, seperated by comma')

    parser.add_argument(
        '--output', dest="output_file",
        type=str,
        help='output bedpe path')

    parser.add_argument(
        '--norm',
        type=str,
        help='Hi-C normalization. Recommend: "balanced", can also be "none"')

    parser.add_argument(
        '--thr', dest='threshold',
        type=float,
        help='P value threshold')

    parser.add_argument(
        '--max_range',
        type=int,
        help='max distance off the diagonal to be calculated')

    parser.add_argument(
        '--resolution',
        type=int,
        help='resolution')

    parser.add_argument(
        '--min_length',
        type=int,
        help='minimum length of stripes')

    parser.add_argument(
        '--min_distance',
        type=int,
        help='threshold for removing stripes too far away from the diagonal')

    parser.add_argument(
        '--max_width',
        type=int,
        default=50000,
        help='maximum stripe width')

    parser.add_argument(
        '--window_size',
        type=int,
        help='size of the window for calculating enrichment score')

    parser.add_argument(
        '--centromere_file',
        type=str,
        default=None,
        help='file for centromere region exclusion')

    parser.add_argument(
        '--nstrata_blank',
        type=int,
        default=0,
        help='bins from main diagonal that are ignored for detection',
        )

    parser.add_argument(
        '--sigma',
        type=float,
        default=12,
        help='Sigma for determining relative peaks in 1D edge detection'
    )

    parser.add_argument(
        '--rel_height',
        type=float,
        default=0.3,
        help='height relative to 1D peak width for stripe width detection'
        )

    parser.add_argument(
        '--presets',
        type=str,
        default="",
        choices=['HiC_5000', 'MiC_5000'],
        help='set preset to \"hiC\" for Hi-C or \"miC\" for micro-C for preset parameter args'
        )

    parser.add_argument(
        '--N_cores',
        type=int,
        default=1,
        help='choose number of CPU cores'
        )

    parser.add_argument(
        '--log_path',
        type=str,
        default='test.txt',
        help='choose number of CPU cores'
    )

    args, unknown = parser.parse_known_args()

    if args.log_path is not None and len(args.log_path) > 0:
        if os.path.exists(args.log_path):
            f = open(args.log_path, 'a')
        else:
            f = open(args.log_path, 'w')
        f.write('Program started with command line...\n')
        f.write('===============\n')
        f.write('#Args passed from command line:\n')
        f.write('#Arg_name Arg_value Arg_PythonReference\n')
        f.write(f'hic {args.hic} {id(args.hic)}\n')
        f.write(f'reference_genome {args.reference_genome} {id(args.reference_genome)}\n')
        f.write(f'chrs {args.chrs} {id(args.chrs)}\n')
        f.write(f'output_file {args.output_file} {id(args.output_file)}\n')
        f.write(f'norm {args.norm} {id(args.norm)}\n')
        f.write(f'threshold {args.threshold} {id(args.threshold)}\n')
        f.write(f'max_range {args.max_range} {id(args.max_range)}\n')
        f.write(f'resolution {args.resolution} {id(args.resolution)}\n')
        f.write(f'min_length {args.min_length} {id(args.min_length)}\n')
        f.write(f'min_distance {args.min_distance} {id(args.min_distance)}\n')
        f.write(f'max_width {args.max_width} {id(args.max_width)}\n')
        f.write(f'window_size {args.window_size} {id(args.window_size)}\n')
        f.write(f'N_cores {args.N_cores} {id(args.N_cores)}\n')
        f.write(f'nstrata_blank {args.nstrata_blank} {id(args.nstrata_blank)}\n')
        f.write(f'sigma {args.sigma} {id(args.sigma)}\n')
        f.write(f'rel_height {args.rel_height} {id(args.rel_height)}\n')
        f.write(f'centromere_file {args.centromere_file} {id(args.centromere_file)}\n')
        if len(unknown):
            f.write(f'The following args are not recognized:{unknown}')
        f.write('===============\n\n')
        f.close()

    if len(unknown):
        print(' The following args are not recognized:', unknown)
    return args


if __name__ == '__main__':
    stripe_parser()



