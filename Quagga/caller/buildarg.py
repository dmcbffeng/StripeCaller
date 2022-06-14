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
        default = ["chr1"],
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
        default=10000000000,
        help='merge stripes which are close to each other (# of bins)')

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
    args, unknown = parser.parse_known_args()

    print(' The following args are not recognized:', unknown)
    return args




