# -*- coding: utf-8 -*-
"""
This file is for loading a single HiC contact map from different type of file.
"""

import numpy as np
from gzip import open as gopen
import os.path
from .hic_straw import straw
from .cool import dump

my_path = os.path.abspath(os.path.dirname(__file__))


def get_chromosome_lengths(ref, res=1):
    """
    Get lengths for all chromosomes in a given resolution according to the reference genome.

    Args:
        ref (str or dict): name of reference genome, eg.'mm10', 'hg38'
        res (int): resolution

    Returns:
        chromosomes (set):
        lengths (dict): eg. {'chr1': 395, 'chr2': 390, ...}
    """

    def _res_change(length, res):
        if length % res == 0:
            return length // res
        else:
            return length // res + 1

    if isinstance(ref, str):
        try:
            path = os.path.join(my_path, 'reference_genome/' + ref)
            rg = open(path)
        except FileNotFoundError:
            try:
                path = os.path.join(my_path, 'reference_genome\\' + ref)
                rg = open(path)
            except FileNotFoundError:
                raise FileNotFoundError('Reference genome {0} not supported!'.format(ref))
        rg = [line.strip().split() for line in rg]
        lengths = {lst[0]: _res_change(int(lst[1]), res) for lst in rg}
    elif isinstance(ref, dict):
        lengths = {elm: _res_change(ref[elm], res) for elm in ref}
    else:
        raise ValueError('Unsupported reference genome!')

    chromosomes = set(lengths.keys())
    return chromosomes, lengths


def file_line_generator(file, format=None, chrom=None, header=0, resolution=1,
                        resolution_adjust=True, mapping_filter=0., gzip=False):
    """
    For formats other than .hic and .mcool

    Args:
        file (str): File path.

        format (int or list): Format
            1) [chr1 pos1 chr2 pos2 mapq1 mapq2];
            2) [chr1 pos1 chr2 pos2 score];
            3) [chr1 pos1 chr2 pos2].

        chrom (str): Chromosome to extract.

        header (None or int): Whether to skip header (1 line).

        mapping_filter (float): The threshold to filter some reads by map quality.


    Returns:
        No return value.
        Yield a line each time in the format of (position_1, position_2, contact_reads).
    """

    f = gopen(file) if gzip else open(file)
    if header:
        for _ in range(header):
            next(f)
    for line in f:
        _line = line.decode('utf-8') if gzip else line
        if _line.startswith('#'):
            continue
        lst = _line.strip().split()
        if len(format) not in [4, 5, 6]:
            raise ValueError('Wrong custom format!')

        if format[0] != 0 and format[2] != 0:
            # chr1 chr2
            c1, c2 = lst[format[0] - 1], lst[format[2] - 1]
            if (c1 != chrom and 'chr' + c1 != chrom) or (c2 != chrom and 'chr' + c2 != chrom):
                continue

        if len(format) == 6:  # [chr1 pos1 chr2 pos2 mapq1 mapq2]
            # mapq1 mapq2
            q1, q2 = float(lst[format[4] - 1]), float(lst[format[5] - 1])
            if q1 < mapping_filter or q2 < mapping_filter:
                continue

        # pos1 pos2
        p1, p2 = int(lst[format[1] - 1]), int(lst[format[3] - 1])
        if resolution_adjust:
            p1 = p1 // resolution  # * resolution
            p2 = p2 // resolution  # * resolution

        if len(format) == 4 or len(format) == 6:  # [chr1 pos1 chr2 pos2]
            v = 1.0
        elif len(format) == 5:
            v = float(lst[format[4] - 1])
        else:
            raise ValueError('Wrong custom format!')

        yield p1, p2, v
    f.close()


def load_HiC(file, ref_genome, format=None,
             chromosome=None, resolution=10000, norm='KR',
             max_distance=5000000, **kwargs):
    """
    Load HiC contact map into a matrix
    Args:
        file (str): File path.
        ref_genome (str): reference genome
        format (str): Now support .txt, .hic, and .mcool file.
        chromosome (str): Specify the chromosome.
        resolution (int): Resolution.
        max_distance (int): max contact distance to keep.

    Return:
        Numpy.array: loaded contact map
    """
    chroms, sizes = get_chromosome_lengths(ref_genome)
    size = sizes[chromosome]
    format = format.lower()
    norm = norm.lower()

    if format in ['hic', '.hic']:
        if norm in ['kr', 'balanced', 'balance']:
            gen = straw('KR', file, chromosome, chromosome, 'BP', resolution)
        elif norm in ['none']:
            gen = straw('NONE', file, chromosome, chromosome, 'BP', resolution)
        else:
            raise ValueError('Unrecognized norm: ' + norm)
    elif format in ['mcool', 'cool', '.cool', '.mcool']:
        gen = dump(file, resolution=resolution, range=chromosome, range2=chromosome, header=header > 0)

    elif format in ['pairs', 'pair', '.pair', '.pairs']:
        if norm in ['kr', 'balanced', 'balance']:
            raise ValueError('Do not provide: ' + norm + ' for .pairs format! ' +
                             'Please generate .hic or .mcool first.')
        elif norm in ['none']:
            gen = file_line_generator(
                file, format=[2, 3, 4, 5], chrom=chromosome, header=0, gzip=file.endswith('gz'),
                resolution=resolution, mapping_filter=0
            )
        else:
            raise ValueError('Unrecognized norm: ' + norm)
    else:
        raise ValueError('Unrecognized format: ' + format)

    n_strata = int(np.ceil(max_distance / resolution))
    length = int(np.ceil(size / resolution))
    strata = [np.zeros((length - i,)) for i in range(n_strata)]
    norm_factors = np.ones((length,))
    recorded = set()

    for p1, p2, c, nx, ny in gen:
        # print(p1, p2, c, nx, ny)
        if abs(p1 - p2) < n_strata:
            if nx * ny > 0:
                val = c / nx / ny
            else:
                val = 0
            # if not 0 <= val < 10000:
            #     print(p1, p2, c, nx, ny)
            strata[abs(p1 - p2)][min(p1, p2)] += val
            if p1 not in recorded:
                norm_factors[p1] = nx
                recorded.add(p1)
                # print(p1, nx)
            if p2 not in recorded:
                norm_factors[p2] = ny
                recorded.add(p2)
                # print(p2, ny)

    return strata, norm_factors

