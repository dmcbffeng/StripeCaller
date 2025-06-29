import sys
sys.path.append("..")
from ..utils.load_HiC import *
from .functions import enrichment_score2, phased_enrichment_score2, find_max_slice, phased_max_slice_arr, merge_positions, get_stripe_and_widths, get_stripe_and_widths_new
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
        norm_factors=None, stats_test_log=({}, {}),
        log_path=None, nbinom_p=None
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
            phased_outputs = pool.starmap(phased_enrichment_score2,
                               zip(repeat(mat), lst, wtd,
                                   repeat(norm_factors),
                                   repeat(targeted_range),
                                   repeat(window_size),
                                   repeat(nbinom_p)))

        arr, all_exp, all_obs, loggerhash = [], [], [], {}
        ctr=0
        for out in phased_outputs:
            arr0, all_exp0, all_obs0, stat_log0 = out
            arr += [arr0]
            all_exp += [all_exp0]
            all_obs += [all_obs0]
            stats_test_log[0].update(stat_log0[0])
            stats_test_log[1].update(stat_log0[1])
            loggerhash[lst[ctr]]={"arr": arr0,"exp": all_exp0,"obs": all_obs0}
            ctr+=1            
        arr += np.log10(threshold)

        with Pool(N) as pool:
            all_positions = (pool.starmap(phased_max_slice_arr, zip(lst, arr, wtd)))

        if log_path is not None and len(log_path) > 0:
            f = open(log_path, 'a')
            for idx, head, tail, _max, width in all_positions:
                f.write(f'Anchor: {idx * resolution}, Width: {int(width)}\n')
                f.write(f'Observed: {list(loggerhash[idx]["obs"])}\n')
                f.write(f'Expected: {list(loggerhash[idx]["exp"])}\n')
                f.write(f'-logPval: {list(loggerhash[idx]["arr"])}\n')
                f.write(f'Head and tail: {head} {tail}\n')
                f.write('===============\n')
            f.close()
                
    else:
        lst = [idx for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]
        wtd = [max(int(positions[idx]), 1) for idx in list(sorted(positions.keys())) if
               not idx <= window_size or not idx >= mat.shape[0] - window_size]

        for i, idx in enumerate(lst):
            if idx <= window_size or idx >= mat.shape[0] - window_size:
                continue

            arr, all_exp, all_obs = enrichment_score2(mat, idx, int(wtd[i]),
                                    distance_range=targeted_range,
                                    window_size=window_size, nbinom_p=nbinom_p,
                                    norm_factors=norm_factors, stats_test_log=stats_test_log
                                    )

            arr = arr + np.log10(threshold)
            head, tail, _max = find_max_slice(arr)
            # if tail - head > 0:
            #     _max = _max / (tail - head) - np.log10(threshold)
            #     all_positions.append((idx, head, tail, _max, wtd[i]))
            all_positions.append((idx, head, tail, _max, wtd[i]))

            if log_path is not None and len(log_path) > 0:
                f = open(log_path, 'a')
                f.write(f'Anchor: {idx * resolution}, Width: {int(wtd[i])}\n')
                f.write(f'Observed: {list(all_obs)}\n')
                f.write(f'Expected: {list(all_exp)}\n')
                f.write(f'-logPval: {list(arr)}\n')
                f.write(f'Head and tail: {head} {tail}\n')
                f.write('===============\n')
                f.close()

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
        gabor_freq=0.1,
        log_path=None,
        nbinom_p=None
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

    if log_path is not None and len(log_path) > 0:
        if os.path.exists(log_path):
            f2 = open(log_path, 'a')
        else:
            f2 = open(log_path, 'w')
        f2.write('Program started with the stripe_caller_all() function...\n')
        f2.write('===============\n')
        f2.write('#Args received by the function:\n')
        f2.write('#Arg_name Arg_value Arg_PythonReference\n')
        f2.write(f'hic {hic_file} {id(hic_file)}\n')
        f2.write(f'reference_genome {reference_genome} {id(reference_genome)}\n')
        f2.write(f'chrs {chromosomes} {id(chromosomes)}\n')
        f2.write(f'output_file {output_file} {id(output_file)}\n')
        f2.write(f'norm {norm} {id(norm)}\n')
        f2.write(f'threshold {threshold} {id(threshold)}\n')
        f2.write(f'max_range {max_range} {id(max_range)}\n')
        f2.write(f'resolution {resolution} {id(resolution)}\n')
        f2.write(f'min_length {min_length} {id(min_length)}\n')
        f2.write(f'min_distance {min_distance} {id(min_distance)}\n')
        f2.write(f'max_width {max_width} {id(max_width)}\n')
        f2.write(f'window_size {window_size} {id(window_size)}\n')
        f2.write(f'N_cores {N_threads} {id(N_threads)}\n')
        f2.write(f'nstrata_blank {nstrata_blank} {id(nstrata_blank)}\n')
        f2.write(f'sigma {sigma} {id(sigma)}\n')
        f2.write(f'rel_height {rel_height} {id(rel_height)}\n')
        f2.write(f'centromere_file {centromere_file} {id(centromere_file)}\n')
        f2.write('===============\n\n')
        f2.close()

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
    f.write('#chr1\tx1\tx2\tchr2\ty1\ty2\tP_val\n')

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

        if log_path is not None and len(log_path) > 0:
            f2 = open(log_path, 'a')
            f2.write('===============\n')
            f2.write(f'Loaded the chromosome: {ch}, sum of contacts: {_sum}\n')
            f2.write('===============\n')
            f2.close()

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
            gabor_freq=gabor_freq, gabor_theta=1
        )
        print(f' {len(h_Peaks)} identified')

        if log_path is not None and len(log_path) > 0:
            f2 = open(log_path, 'a')
            f2.write('===============\n')
            f2.write('Horizontal:')
            f2.write(f'{len(h_Peaks)} candidate stripes identified with Gabor filtering and peak calling\n')
            f2.write('Starting the statistical test for each candidate stripe')
            f2.write('===============\n')
            f2.close()

        if h_Peaks:
            results = _stripe_caller(mat, positions=h_Peaks, threshold=threshold,
                                     max_range=max_range, resolution=resolution,
                                     min_length=min_length, closeness=min_distance,
                                     window_size=window_size, N=N_threads,
                                     norm_factors=norm_factors, stats_test_log=(_calculated_values, _poisson_stats),
                                     log_path=log_path, nbinom_p=nbinom_p
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
                f.write(f'{ch}\t{st*resolution}\t{ed*resolution}\t{ch}\t{max((st+hd), ed)*resolution}\t{(ed+tl)*resolution}\t{10 ** (- sc)}\n')

        # vertical
        print(' Vertical:')
        mat = strata2vertical(strata)
        # np.save('vmat.npy', mat)

        print(' Finding candidate peaks:')
        v_Peaks = get_stripe_and_widths_new(
            mat, (max_range + min_length) // resolution, nstrata_blank,
            sigma=sigma, rel_height=rel_height, max_width=max_width // resolution,
            gabor_freq=gabor_freq, gabor_theta=1
        )
        print(f' {len(v_Peaks)} identified')

        if log_path is not None and len(log_path) > 0:
            f2 = open(log_path, 'a')
            f2.write('===============\n')
            f2.write('Vertical:')
            f2.write(f'{len(h_Peaks)} candidate stripes identified with Gabor filtering and peak calling\n')
            f2.write('Starting the statistical test for each candidate stripe')
            f2.write('===============\n')
            f2.close()

        if v_Peaks:
            results = _stripe_caller(mat, positions=v_Peaks, threshold=threshold,
                                     max_range=max_range, resolution=resolution,
                                     min_length=min_length, closeness=min_distance,
                                     window_size=window_size, N=N_threads,
                                     norm_factors=norm_factors, stats_test_log=(_calculated_values, _poisson_stats),
                                     log_path=log_path, nbinom_p=nbinom_p
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
                f.write(f'{ch}\t{(st-tl)*resolution}\t{min((ed-hd), st)*resolution}\t{ch}\t{st*resolution}\t{ed*resolution}\t{10 ** (- sc)}\n')
    f.close()

