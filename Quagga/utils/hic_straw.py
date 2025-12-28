# Quagga/utils/hic_straw_adapter.py
"""
Adapter for modern hicstraw library (v1.3.1+).

The bundled hic_straw.py cannot parse .hic files created by newer versions
of Juicer. This adapter wraps the modern hicstraw library to provide the
same interface that load_HiC.py expects.

Modern hicstraw returns pre-normalized contact values, so normalization
factors are returned as 1.0 (identity) since division is already applied.

Note: Modern hicstraw returns genomic coordinates (base pairs), while the
legacy straw returned bin indices. This adapter converts coordinates back
to bin indices for compatibility with load_HiC.py.
"""

import hicstraw


def straw(norm, infile, chr1loc, chr2loc, unit, binsize, is_synapse=False):
    """
    Extract contact records from a .hic file.
    
    Wraps modern hicstraw.straw() to match the interface expected by
    Quagga's load_HiC module.
    
    Args:
        norm (str): Normalization type - 'NONE', 'VC', 'VC_SQRT', or 'KR'
        infile (str): Path to .hic file
        chr1loc (str): Chromosome with optional range, e.g. 'chr1' or 'chr1:0:1000000'
        chr2loc (str): Chromosome with optional range
        unit (str): 'BP' (base pairs) or 'FRAG' (fragments)
        binsize (int): Resolution in base pairs
        is_synapse (bool): Legacy parameter, ignored
    
    Yields:
        tuple: (bin1, bin2, count, norm1, norm2)
            bin1 (int): Bin index for first position
            bin2 (int): Bin index for second position  
            count (float): Contact count (normalized if norm != 'NONE')
            norm1 (float): 1.0 (normalization pre-applied by modern library)
            norm2 (float): 1.0 (normalization pre-applied by modern library)
    """
    # Modern hicstraw signature differs from legacy:
    #   Legacy:  straw(norm, file, chr1, chr2, unit, binsize)
    #   Modern:  straw(data_type, norm, file, chr1, chr2, unit, binsize)
    records = hicstraw.straw('observed', norm, infile, chr1loc, chr2loc, unit, binsize)
    
    for record in records:
        # Modern hicstraw returns genomic coordinates (bp), not bin indices.
        # Convert to bin indices for compatibility with load_HiC.py.
        bin1 = record.binX // binsize
        bin2 = record.binY // binsize
        yield bin1, bin2, record.counts, 1.0, 1.0