from ..utils.load_HiC import load_HiC
from .functions import enrichment_score2, find_max_slice, phased_max_slice_arr, merge_positions, get_stripe_and_widths
from .mat_ops import strata2vertical, strata2horizontal, strata2triu, blank_diagonal

import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import repeat


def _stripe_caller(
        mat, positions,
        max_range=150000, resolution=1000,
        min_length=30000, closeness=50000,
        merge=1, window_size=8, threshold=0.01,
        N=1
):
    assert max_range % resolution == 0
    assert min_length % resolution == 0

    def pack_tuple(*args):
        return (*args,)

    all_positions = []
    #
    # f = open(f'h_arr_chr1.txt', 'w')

    targeted_range = pack_tuple(0, max_range // resolution)
    if N > 1:  # parallel if #CPUs set
        lst = [idx for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]
        wtd = [int(positions[idx]) for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]

        with Pool(N) as pool:
            # arr = pool.starmap(enrichment_score2, zip(lst, wtd, len(lst)*[targeted_range], len(lst)*[window_size]))
            arr = pool.starmap(enrichment_score2,
                               zip(repeat(mat), lst, wtd, repeat(targeted_range), repeat(window_size)))
        arr += np.log10(threshold)

        with Pool(N) as pool:
            all_positions = (pool.starmap(phased_max_slice_arr, zip(lst, arr)))

    else:
        lst = [idx for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]
        wtd = [int(positions[idx]) for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]
        #         print(lst,wtd)
        for i, idx in enumerate(lst):
            if idx <= window_size or idx >= mat.shape[0] - window_size:
                continue
            # print(idx, lst[i], wtd[i])
            # targeted_range = pack_tuple(0, max_range // resolution)
            arr = enrichment_score2(mat, idx, int(wtd[i]),
                                    distance_range=targeted_range,
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
        merge=1, window_size=8,
        centromere_file=None,
        nstrata_blank=0, step=36, sigma=12., rel_height=0.3
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

    centro = {}
    if centromere_file is not None:
        for line in open(centromere_file):
            [ch, st, ed] = line.strip().split()[:3]
            st, ed = int(st), int(ed)
            assert ch.startswith('chr')
            if ch not in centro:
                centro[ch] = []
            centro[ch].append((st, ed))

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
        # TO DO: max distance!!!

        # full mat for calling candidate stripes
        mat = strata2triu(strata)
        mat = blank_diagonal(mat, nstrata_blank)
        h_Peaks, v_Peaks = get_stripe_and_widths(
            mat, step=step, sigma=sigma, rel_height=rel_height
        )

        # horizontal
        mat = strata2horizontal(strata)
        results = _stripe_caller(mat, positions=h_Peaks, threshold=threshold,
                                 max_range=max_range, resolution=resolution,
                                 min_length=min_length, closeness=min_distance,
                                 merge=merge, window_size=window_size)
        for (st, ed, hd, tl, sc) in results:
            in_centro = False
            if ch in centro:
                for (centro_st, centro_ed) in centro[ch]:
                    if centro_st <= st * resolution <= centro_ed or centro_st <= ed * resolution <= centro_ed:
                        in_centro = True
            if not in_centro:
                f.write(f'{ch}\t{st*resolution}\t{ed*resolution}\t{ch}\t{max((st+hd), ed)*resolution}\t{(ed+tl)*resolution}\t{sc}\n')

        # vertical
        mat = strata2vertical(strata)
        results = _stripe_caller(mat, positions=v_Peaks, threshold=threshold,
                                 max_range=max_range, resolution=resolution,
                                 min_length=min_length, closeness=min_distance,
                                 merge=merge, window_size=window_size)
        for (st, ed, hd, tl, sc) in results:
            in_centro = False
            if ch in centro:
                for (centro_st, centro_ed) in centro[ch]:
                    if centro_st <= st * resolution <= centro_ed or centro_st <= ed * resolution <= centro_ed:
                        in_centro = True
            if not in_centro:
                f.write(f'{ch}\t{(st-tl)*resolution}\t{min((ed-hd), st)*resolution}\t{ch}\t{st*resolution}\t{ed*resolution}\t{sc}\n')
    f.close()


