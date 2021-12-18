import os
from hic_straw import straw


if __name__ == '__main__':
    hic_path = '../../../GM12878.hic'
    for a, b, c in straw(norm='KR', infile=hic_path, chr1loc='chr1:10000:100000', chr2loc='chr1:10000:100000',
                         unit='BP', binsize=5000, is_synapse=False):
        print(a, b, c)

    for a, b, c in straw(norm='KR', infile=hic_path, chr1loc='chr1:10000:100000', chr2loc='chr1:10000:100000',
                         unit='BP', binsize=5000, is_synapse=False):
        print(a, b, c)
