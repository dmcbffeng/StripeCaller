from ..utils.load_HiC import load_HiC
from .functions import pick_max_positions, enrichment_score, find_max_slice, merge_positions, strata2horizontal, strata2vertical
import numpy as np


def _stripe_caller(mat, max_range=150000, resolution=1000,
                   min_length=30000, closeness=50000,
                   stripe_width=1, merge=1, window_size=8, threshold=0.01):
    assert max_range % resolution == 0
    assert min_length % resolution == 0

    # Step 2: for different distance ranges pick the "local maximum" positions
    print(' Finding local maximum for different contact distances...')
    positions = {}
    # Split the max range into small distance ranges
    for dis in range(0, max_range, min_length):
        _min = dis
        if dis + 2 * min_length > max_range:
            _max = max_range
        else:
            _max = dis + min_length
        print(f'  {_min}-{_max}', end=' ')
        distance_range = (_min // resolution, _max // resolution)
        pos_h = pick_max_positions(mat, distance_range=distance_range, line_width=stripe_width, window_size=window_size)
        print(len(pos_h))
        for p in pos_h:
            if p not in positions:
                positions[p] = []
            positions[p].append(distance_range)
    print('  Total:', len(positions))

    # Step 3: find the accurate range of stripe
    print(' Finding the spanning range for each stripe...')
    all_positions = []
    lst = sorted(positions.keys())
    for i, idx in enumerate(lst):
        # print(i, idx)
        if idx <= window_size or idx >= mat.shape[0] - window_size:
            continue
        arr = enrichment_score(mat, idx, line_width=stripe_width, distance_range=(0, max_range // resolution),
                               window_size=window_size)
        arr = arr + np.log10(threshold)
        head, tail, _max = find_max_slice(arr)
        all_positions.append((idx, head, tail, _max))

    # Step 4: Merging
    print(' Merging...')
    if not all_positions:
        raise ValueError("No statistically significant candidate stripes found(enrichment_score()). Try different args: stripe_width, max_range, resolution, window_size")
    all_positions = merge_positions(all_positions, merge)
    print(len(all_positions))

    print(' Filtering by distance and length ...')
    new_positions = []
    for elm in all_positions:
        # print(elm, end=' ')
        if (elm[3] - elm[2]) * resolution >= min_length and elm[2] * resolution <= closeness:
            # print(True)
            new_positions.append(elm)
        else:
            # print(False)
            pass
    print(len(new_positions))

    # Step 5: Statistical test
    results = []
    print(' Statistical Tests...')
    for elm in new_positions:
        [st, ed, head, tail, score] = elm
        # p = stat_test(mat, st, ed, stripe_width, head, tail, window_size)
        # print(idx * resolution, p)
        if score > threshold:
            results.append((st, (ed + 1), head, tail, score))
    print(len(results))
    return results


def stripe_caller_all(
        hic_file, reference_genome,
        chromosomes,
        output_file,
        norm='balanced',
        threshold=0.01,
        max_range=150000, resolution=1000,
        min_length=30000, min_distance=50000,
        stripe_width=1, merge=1, window_size=8
):
    """
    The main function for calling stripes
    Args:
        hic_file (str): file path
        reference_genome (str): reference genome
        chromosomes (list): which chromosomes to calculate
        output_file (str): output bedpe path
        norm (str): recommend: "balanced", can also be "none"
        threshold (float): p value threshold
        max_range (int): max distance off the diagonal to be calculated
        resolution (int): resolution
        min_length (int): minimum length of stripes
        min_distance (int): threshold for removing stripes too far away from the diagonal
        stripe_width (int): stripe width (# of bins at the given resolution)
        merge (int): merge stripes which are close to each other (# of bins)
        window_size (int): size of the window for calculating enrichment score

    """
    if hic_file.lower().endswith('hic'):
        _format = 'hic'
    elif hic_file.lower().endswith('cool'):
        _format = 'cool'
    elif hic_file.lower().endswith('pairs') or hic_file.lower().endswith('pairs.gz'):
        _format = 'pairs'
    else:
        raise ValueError('Unrecognized format for: ' + hic_file)

    f = open(output_file, 'w')
    f.write('#chr1\tx1\tx2\tchr2\ty1\ty2\tenrichment\n')

    for ch in chromosomes:
        print(f'Calling for {ch}...')
        strata, norm_factors = load_HiC(
            file=hic_file, ref_genome=reference_genome, format=_format,
            chromosome=ch, resolution=resolution,
            norm=norm, max_distance=max_range + min_length)

        # horizontal
        mat = strata2horizontal(strata)
        results = _stripe_caller(mat, threshold=threshold,
                                 max_range=max_range, resolution=resolution,
                                 min_length=min_length, closeness=min_distance,
                                 stripe_width=stripe_width, merge=merge, window_size=window_size)
        for (st, ed, hd, tl, sc) in results:
            f.write(f'{ch}\t{st*resolution}\t{ed*resolution}\t{ch}\t{max((st+hd), ed)*resolution}\t{(ed+tl)*resolution}\t{sc}\n')

        # vertical
        mat = strata2vertical(strata)
        results = _stripe_caller(mat, threshold=threshold,
                                 max_range=max_range, resolution=resolution,
                                 min_length=min_length, closeness=min_distance,
                                 stripe_width=stripe_width, merge=merge, window_size=window_size)
        for (st, ed, hd, tl, sc) in results:
            f.write(f'{ch}\t{(st-tl)*resolution}\t{min((ed-hd), st)*resolution}\t{ch}\t{st*resolution}\t{ed*resolution}\t{sc}\n')

    f.close()


