import os
import re
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import kruskal, poisson
import matplotlib.pyplot as plt
from AVL_tree import AVLTree


stats_log = []
_calculated_values = {}
_poisson_stats = {}


def hic2txt(hic_file, ch, resolution=1000, output='temp.txt'):
    """
    Dump .hic file into contact lists
    :param hic_file: (str) .hic file path
    :param ch: (str) chromosome
    :param resolution: (int) resolution to use
    :param output: (str) temporary output path
    """
    juicer = '/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar'
    # cmd = f'java -jar {juicer} dump observed KR {hic_file} {ch} {ch} BP {resolution} {output}'
    cmd = f'java -jar {juicer} dump observed KR {hic_file} {ch} {ch} BP {resolution} {output}'
    os.system(cmd)


def load_chrom_sizes(reference_genome):
    """
    Load chromosome sizes for a reference genome
    """
    my_path = os.path.abspath(os.path.dirname(__file__))
    # rg_path = f'{my_path}/reference_genome/{reference_genome}'
    rg_path = f'{my_path}/{reference_genome}'
    f = open(rg_path)
    lengths = {}
    for line in f:
        [ch, l] = line.strip().split()
        lengths[ch] = int(l)
    return lengths


def txt2horizontal(txt, length, max_range, resolution=1000):
    """
        :param txt: str, path of input .txt file
        :param length: chromosome length
        :param max_range: int, max distance
        :param resolution: int, default: 25000
        """
    assert max_range % resolution == 0
    f = open(txt)
    n_bins = length // resolution + 1
    rg = max_range // resolution
    mat = np.zeros((n_bins, rg))
    cnt = 0
    for line in f:
        if cnt % 5000000 == 0:
            print('  ', cnt)
        cnt += 1
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)
        if max(p1, p2) >= n_bins * resolution:
            continue
        if p1 > p2:
            p1, p2 = p2, p1
        p1, p2 = p1 // resolution, p2 // resolution
        if p2 - p1 >= rg:
            continue
        mat[p1, p2 - p1] += v
    return mat


def txt2vertical(txt, length, max_range, resolution=1000):
    """
        :param txt: str, path of input .txt file
        :param length: chromosome length
        :param max_range: int, max distance
        :param resolution: int, default: 25000
        """
    assert max_range % resolution == 0
    f = open(txt)
    n_bins = length // resolution + 1
    rg = max_range // resolution
    mat = np.zeros((n_bins, rg))
    cnt = 0
    for line in f:
        if cnt % 5000000 == 0:
            print('  ', cnt)
        cnt += 1
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)
        if max(p1, p2) >= n_bins * resolution:
            continue
        if p1 > p2:
            p1, p2 = p2, p1
        p1, p2 = p1 // resolution, p2 // resolution
        if p2 - p1 >= rg:
            continue
        mat[p2, p2 - p1] += v
    return mat


def pick_max_positions2(mat, distance_range=(10, 160), line_width=1, window_size=10):
    st, ed = distance_range
    stats = np.sum(mat[:, st:ed], axis=1)
    all_pos = []

    all_peaks, _ = find_peaks(stats, distance=window_size * 2)
    for idx in all_peaks:
        check = enrichment_score(mat, idx, line_width, (st, ed), window_size)
        if np.sum(check) > 0:
            all_pos.append(idx)
    return all_pos


def enrichment_score(mat, idx, line_width=1, distance_range=(20, 40), window_size=10):
    # st, ed = max(distance_range[0], window_size), min(distance_range[1], mat.shape[1] - window_size)
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):
        if j < window_size + half or j >= mat.shape[1] - window_size - half:
            continue
        y = j - distance_range[0]
        line_min = min(np.mean(mat[x1:x2, j-window_size//2-half:j-half]),
                       np.mean(mat[x1:x2, j+1+half:j+window_size//2+half+1]))
        neighbor_mean = max(np.mean(mat[idx-window_size:x1, j-window_size//2-half:j+window_size//2+half+1]),
                            np.mean(mat[x2+1:idx+window_size+1, j-window_size//2-half:j+window_size//2+half+1]))
        new_mat[y] = line_min - neighbor_mean
    return new_mat


def enrichment_score2(mat, idx, line_width=1, distance_range=(5, 100), window_size=10):
    # st, ed = max(distance_range[0], window_size), min(distance_range[1], mat.shape[1] - window_size)
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):
        # print(j)
        if j < window_size + half or j >= mat.shape[1] - window_size - half:
            continue
        y = j - distance_range[0]
        line_min = np.median(mat[x1:x2, j - window_size//2 - half:j + window_size//2 + half + 1])
        neighbor_mean = max(np.mean(mat[idx-window_size:x1, j-window_size//2-half:j+window_size//2+half+1]),
                            np.mean(mat[x2+1:idx+window_size+1, j-window_size//2-half:j+window_size//2+half+1]))
        # print(line_min, neighbor_mean)
        # new_mat[y] = line_min / (neighbor_mean + 1e-9) - 1
        lower_b = 1  # This should be (1 / KR_norm_factors)
        _exp = max(neighbor_mean, lower_b)
        _obs = int(line_min)

        global _calculated_values, _poisson_stats, stats_log
        tolerance = 0.01

        # check if obs is calculated before
        if _obs in _calculated_values:
            (_upper, _lower) = _calculated_values[_obs].search(_exp)
            if _upper is not None and (_upper.key - _exp) < tolerance * _exp:
                _exp = _upper.key
                _exp_idx = _upper.val
                mlog_p_val = _poisson_stats[(_exp_idx, _obs)]
            else:
                _exp_idx = _calculated_values[_obs].insert(_exp)
                Poiss = poisson(_exp)
                p_val = 1 - Poiss.cdf(_obs)
                if 0 < p_val < 1:
                    mlog_p_val = - np.log10(p_val)
                else:  # Some p values are too small, -log(0) will return an error, so we use -1 to temporarily replace
                    mlog_p_val = -1
                _poisson_stats[(_exp_idx, _obs)] = mlog_p_val
                stats_log.append([_exp, _obs, mlog_p_val])
        else:
            _calculated_values[_obs] = AVLTree()
            _exp_idx = _calculated_values[_obs].insert(_exp)
            Poiss = poisson(_exp)
            p_val = 1 - Poiss.cdf(_obs)
            if 0 < p_val < 1:
                mlog_p_val = - np.log10(p_val)
            else:  # Some p values are too small, -log(0) will return an error, so we use -1 to temporarily replace
                mlog_p_val = -1
            _poisson_stats[(_exp_idx, _obs)] = mlog_p_val
            stats_log.append([_exp, _obs, mlog_p_val])

        new_mat[y] = mlog_p_val
    new_mat[new_mat < 0] = np.max(new_mat)  # Replace all "-1"s with the largest -log(p)
    return new_mat


def find_max_slice(arr):
    _max, head, tail = 0, 0, 0
    _max_ending, h, t = 0, 0, 0
    i = 0
    while i < len(arr):
        _max_ending = _max_ending + arr[i]
        if _max_ending < 0:
            h, t = i + 1, i + 1
            _max_ending = 0
        else:
            t = i + 1
        if _max_ending > _max:
            head, tail, _max = h, t, _max_ending
        i += 1
    return head, tail, _max


def merge_positions(lst, merge_range):
    def _merge(small_lst):
        st = min([elm[0] for elm in small_lst])
        ed = max([elm[1] for elm in small_lst])
        head = min([elm[2] for elm in small_lst])
        tail = max([elm[3] for elm in small_lst])
        score = max([elm[4] for elm in small_lst])
        return [st, ed, head, tail, score]

    new_lst = []
    temp = []
    for i, (idx, head, tail, score) in enumerate(lst):
        if i == 0:
            temp.append([idx, idx, head, tail, score])
        elif idx - temp[-1][1] <= merge_range:
            temp.append([idx, idx, head, tail, score])
        else:
            new_lst.append(_merge(temp))
            temp = [[idx, idx, head, tail, score]]
    new_lst.append(_merge(temp))
    return new_lst


def _stripe_caller(mat, max_range=150000, resolution=1000,
                   min_length=30000, closeness=50000,
                   stripe_width=1, merge=1, window_size=8, threshold=0.01):
    assert max_range % resolution == 0
    assert min_length % resolution == 0

    # Step 2: for different distance ranges pick the "local maximum" positions
    print(' Finding local maximum for different contact distances...')
    positions = {}
    # Split the max range into small distance ranges
    for dis in range(0, max_range - min_length + 1, min_length):
        _min = dis
        if dis + 2 * min_length > max_range:
            _max = max_range
        else:
            _max = dis + min_length
        print(f'  {_min}-{_max}', end=' ')
        distance_range = (_min // resolution, _max // resolution)
        pos_h = pick_max_positions2(mat, distance_range=distance_range, line_width=stripe_width,
                                    window_size=window_size)
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
        arr = enrichment_score2(mat, idx, line_width=stripe_width,
                                distance_range=(0, max_range // resolution),
                                window_size=window_size)
        arr = arr + np.log10(threshold)
        head, tail, _max = find_max_slice(arr)
        all_positions.append((idx, head, tail, _max))

    # Step 4: Merging
    print(' Merging...')
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
        # if score > threshold:
        results.append((st, (ed + 1), head, tail, score / (tail - head) - np.log10(threshold)))
    print(len(results))
    return results


def stripe_caller_all(
        hic_file,
        chromosomes,
        output_file,
        threshold=50.,
        max_range=150000, resolution=1000,
        min_length=30000, closeness=50000,
        stripe_width=1, merge=1, window_size=8,
        centromere_file=None
):
    ch_sizes = load_chrom_sizes('hg38')

    f = open(output_file, 'w')
    f.write('#chr1\tx1\tx2\tchr2\ty1\ty2\tenrichment\n')

    centro = {}
    if centromere_file is not None:
        for line in open(centromere_file):
            [ch, st, ed] = line.strip().split()[:3]
            st, ed = int(st), int(ed)
            assert ch.startswith('chr')
            if ch not in centro:
                centro[ch] = []
            centro[ch].append((st, ed))

    for ch in chromosomes:
        # if ch == 'chr14':
        #     continue

        print(f'Calling for {ch}...')
        hic2txt(hic_file, ch, resolution=resolution, output='temp.txt')

        # horizontal
        mat = txt2horizontal('temp.txt', length=ch_sizes[ch], max_range=max_range + min_length, resolution=resolution)
        if np.sum(mat) == 0:
            continue
        results = _stripe_caller(mat, threshold=threshold,
                                 max_range=max_range, resolution=resolution,
                                 min_length=min_length, closeness=closeness,
                                 stripe_width=stripe_width, merge=merge, window_size=window_size)

        for (st, ed, hd, tl, sc) in results:
            in_centro = False
            if ch in centro:
                for (centro_st, centro_ed) in centro[ch]:
                    if centro_st <= st * resolution <= centro_ed or centro_st <= ed * resolution <= centro_ed:
                        in_centro = True
            if not in_centro:
                f.write(f'{ch}\t{st * resolution}\t{ed * resolution}\t{ch}\t{max((st + hd), ed) * resolution}\t{(ed + tl) * resolution}\t{np.power(10, -sc)}\n')

        # vertical
        mat = txt2vertical('temp.txt', length=ch_sizes[ch], max_range=max_range + min_length, resolution=resolution)
        results = _stripe_caller(mat, threshold=threshold,
                                 max_range=max_range, resolution=resolution,
                                 min_length=min_length, closeness=closeness,
                                 stripe_width=stripe_width, merge=merge, window_size=window_size)
        for (st, ed, hd, tl, sc) in results:
            in_centro = False
            if ch in centro:
                for (centro_st, centro_ed) in centro[ch]:
                    if centro_st <= st * resolution <= centro_ed or centro_st <= ed * resolution <= centro_ed:
                        in_centro = True
            if not in_centro:
                f.write(f'{ch}\t{max(0, (st - tl)) * resolution}\t{min((ed - hd), st) * resolution}\t{ch}\t{st * resolution}\t{ed * resolution}\t{np.power(10, -sc)}\n')

        os.remove('temp.txt')

    f.close()


def visualize_stats_log(stats_log, name):
    plt.figure(figsize=(12, 12))
    stats_log = np.array(stats_log)
    plt.scatter(stats_log[:, 0], stats_log[:, 1], c=stats_log[:, 1], s=1)
    plt.xlabel('Expected')
    plt.ylabel('Observed')
    plt.colorbar()
    plt.title(f'{len(stats_log)} Poisson Tests')
    plt.savefig(name)
    plt.close()


if __name__ == '__main__':
    import time

    chromosomes = [f'chr{i}' for i in list(range(1, 23)) + ['X']]
    # chromosomes = ['chr1', 'chr2']

    start_time = time.time()
    hic_file = '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/GM12878/GM12878.hic'
    thr = 0.15
    stripe_caller_all(
        hic_file=hic_file,
        chromosomes=chromosomes,
        output_file='GM12878_HiC_stripes_all.bedpe',
        threshold=thr,
        max_range=2000000, resolution=5000,
        min_length=300000, closeness=500000,
        stripe_width=3, merge=3, window_size=10,
        centromere_file='removed_regions.bed'
    )
    print("bulkHiC/GM12878/GM12878.hic took:\n--- %s seconds ---" % (time.time() - start_time))
    # 87 seconds
    visualize_stats_log(stats_log, 'GM_all_stats_log.png')

    start_time = time.time()
    hic_file = '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/H1-hESC/raw/4DNFIQYQWPF5.hic'
    thr = 0.15
    stripe_caller_all(
        hic_file=hic_file,
        chromosomes=chromosomes,
        output_file='H1_HiC_stripes_all.bedpe',
        threshold=thr,
        max_range=2000000, resolution=5000,
        min_length=300000, closeness=500000,
        stripe_width=3, merge=3, window_size=10,
        centromere_file='removed_regions.bed'
    )
    print("bulkHiC/H1/H1.hic took:\n--- %s seconds ---" % (time.time() - start_time))
    # 87 seconds
    visualize_stats_log(stats_log, 'H1_all_stats_log.png')


