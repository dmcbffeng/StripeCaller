import sys
sys.path.append("..")
from ..utils.load_HiC import *
from .functions import enrichment_score2, find_max_slice, phased_max_slice_arr, merge_positions, get_stripe_and_widths, get_stripe_and_widths_new
from .mat_ops import strata2vertical, strata2horizontal, blank_diagonal_sparse_from_strata

import numpy as np
from multiprocessing import Pool, cpu_count
# from functools import partial
from itertools import repeat


__version__ = '0.0.1'


def _stripe_caller(
        mat, positions,
        max_range=150000, resolution=1000,
        min_length=30000, closeness=50000,
        window_size=8, threshold=0.01,
        N=1,
        norm_factors=None, stats_test_log=({}, {})
):
    assert max_range % resolution == 0
    assert min_length % resolution == 0

    if norm_factors is None:
        norm_factors = np.ones((len(mat),))

    def pack_tuple(*args):
        return (*args,)

    # Step 5: Statistical test
    print(' Statistical Tests...')
    all_positions = []
    targeted_range = pack_tuple(0, max_range // resolution)

    if N > 1:  # parallel if #CPUs set
        lst = [idx for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]
        wtd = [max(int(positions[idx]), 1) for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]

        with Pool(N) as pool:
            # arr = pool.starmap(enrichment_score2, zip(lst, wtd, len(lst)*[targeted_range], len(lst)*[window_size]))
            arr = pool.starmap(enrichment_score2,
                               zip(repeat(mat), lst, wtd, repeat(norm_factors), repeat(targeted_range), repeat(window_size)))
        arr += np.log10(threshold)

        with Pool(N) as pool:
            all_positions = (pool.starmap(phased_max_slice_arr, zip(lst, arr, wtd)))

    else:
        lst = [idx for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]
        wtd = [max(int(positions[idx]), 1) for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]

        for i, idx in enumerate(lst):
            if idx <= window_size or idx >= mat.shape[0] - window_size:
                continue
            arr = enrichment_score2(mat, idx, int(wtd[i]),
                                    distance_range=targeted_range,
                                    window_size=window_size,
                                    norm_factors=norm_factors, stats_test_log=stats_test_log
                                    )

            arr = arr + np.log10(threshold)
            head, tail, _max = find_max_slice(arr)
            all_positions.append((idx, head, tail, _max, wtd[i]))

    if not all_positions:
        raise ValueError("No statistically significant candidate stripes found(enrichment_score()). "
                         "Try different args: stripe_width, max_range, resolution, window_size")
    all_positions = merge_positions(all_positions)
    print(len(all_positions))  # [st, ed, head, tail, score]

    # Step 5: filter by length
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
    return new_positions


def stripe_caller_all(
        hic_file, reference_genome,
        chromosomes,
        output_file,
        norm='balanced',
        threshold=0.01,
        max_range=150000, resolution=1000,
        min_length=30000, min_distance=50000,
        max_width=25000, window_size=8,
        centromere_file=None,
        N_threads=1,
        nstrata_blank=0,
        sigma=12.,
        rel_height=0.3,
        gabor_freq=0.1
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
        max_width (int): maximum width of stripes
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

    # Stats test record
    _calculated_values = {}
    _poisson_stats = {}

    for ch in chromosomes:
        print(f'Calling for {ch}...')
        print(' Loading contact matrix...')
        strata, norm_factors = load_HiC(
            file=hic_file, ref_genome=reference_genome, format=_format,
            chromosome=ch, resolution=resolution, norm=norm,
            max_distance=max_range + min_length
        )
        # np.savetxt('norm_factors.txt', norm_factors)
        print(' Finish loading contact matrix...')

        # check whether the chromosome is empty
        _sum = sum([np.sum(elm) for elm in strata])
        if _sum < 1e-5:
            print(f' Warning: The matrix for {ch} with {norm} normalization is empty. Skipping the chromosome.')
            continue

        # full mat for calling candidate stripes
        # print(' Finding candidate peaks:')
        # mat = blank_diagonal_sparse_from_strata(strata, nstrata_blank)
        # h_Peaks, v_Peaks = get_stripe_and_widths(
        #     mat, step=step, sigma=sigma, rel_height=rel_height
        # )
        # print('  H:', len(h_Peaks), ', V:', len(v_Peaks))
        
        # horizontal
        print(' Horizontal:')
        mat = strata2horizontal(strata)
        # np.save('hmat.npy', mat)

        print(' Finding candidate peaks:')
        h_Peaks = get_stripe_and_widths_new(
            mat, (max_range + min_length) // resolution, nstrata_blank,
            sigma=sigma, rel_height=rel_height, max_width=max_width // resolution,
            gabor_freq=gabor_freq, gabor_theta=0
        )
        print(f' {len(h_Peaks)} identified')

        if h_Peaks:
            results = _stripe_caller(mat, positions=h_Peaks, threshold=threshold,
                                     max_range=max_range, resolution=resolution,
                                     min_length=min_length, closeness=min_distance,
                                     window_size=window_size, N=N_threads,
                                     norm_factors=norm_factors, stats_test_log=(_calculated_values, _poisson_stats)
                                     )
        else:
            results = []
        for (st, ed, hd, tl, sc) in results:
            in_centro = False
            if ch in centro:
                for (centro_st, centro_ed) in centro[ch]:
                    if centro_st <= st * resolution <= centro_ed or centro_st <= ed * resolution <= centro_ed:
                        in_centro = True
            if not in_centro:
                f.write(f'{ch}\t{st*resolution}\t{ed*resolution}\t{ch}\t{max((st+hd), ed)*resolution}\t{(ed+tl)*resolution}\t{sc}\n')

        # vertical
        print(' Vertical:')
        mat = strata2vertical(strata)
        # np.save('vmat.npy', mat)

        print(' Finding candidate peaks:')
        v_Peaks = get_stripe_and_widths_new(
            mat, (max_range + min_length) // resolution, nstrata_blank,
            sigma=sigma, rel_height=rel_height, max_width=max_width // resolution,
            gabor_freq=gabor_freq, gabor_theta=0
        )
        print(f' {len(v_Peaks)} identified')

        if v_Peaks:
            results = _stripe_caller(mat, positions=v_Peaks, threshold=threshold,
                                     max_range=max_range, resolution=resolution,
                                     min_length=min_length, closeness=min_distance,
                                     window_size=window_size, N=N_threads,
                                     norm_factors=norm_factors, stats_test_log=(_calculated_values, _poisson_stats)
                                     )
        else:
            results = []

        # f2 = open(f'peaks_{ch}.txt', 'w')
        # f2.write('H\n')
        # for h in h_Peaks:
        #     f2.write(f'{h * resolution}\t{h_Peaks[h]}\n')
        # f2.write('V\n')
        # for v in v_Peaks:
        #     f2.write(f'{v * resolution}\t{v_Peaks[v]}\n')
        # f2.close()

        # f2 = open(f'peaks_{ch}.bedpe', 'w')
        # # f2.write('H\n')
        # for h in h_Peaks:
        #     width = max(int(h_Peaks[h]), 1)
        #     f2.write(f'{ch}\t{(h - width // 2) * resolution}\t{(h + width // 2 + 1) * resolution}\t{ch}\t{(h + 10) * resolution}\t{(h + 1000) * resolution}\t{1.0}\n')
        # # f2.write('V\n')
        # for v in v_Peaks:
        #     width = max(int(v_Peaks[v]), 1)
        #     f2.write(f'{ch}\t{(v - 1000) * resolution}\t{(v - 10) * resolution}\t{ch}\t{(v - width // 2) * resolution}\t{(v + width // 2 + 1) * resolution}\t{0.0}\n')
        # f2.close()

        for (st, ed, hd, tl, sc) in results:
            in_centro = False
            if ch in centro:
                for (centro_st, centro_ed) in centro[ch]:
                    if centro_st <= st * resolution <= centro_ed or centro_st <= ed * resolution <= centro_ed:
                        in_centro = True
            if not in_centro:
                f.write(f'{ch}\t{(st-tl)*resolution}\t{min((ed-hd), st)*resolution}\t{ch}\t{st*resolution}\t{ed*resolution}\t{sc}\n')
    f.close()

