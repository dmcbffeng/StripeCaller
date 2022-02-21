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
        type=str,
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
        '--merge',
        type=int,
        default=1,
        help='merge stripes which are close to each other (# of bins)')

    parser.add_argument(
        '--window_size',
        type=int,
        default=1,
        help='size of the window for calculating enrichment score')

    return parser

