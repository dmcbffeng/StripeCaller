import numpy as np
from scipy.signal import find_peaks
from scipy.stats import poisson


def hic_loader(
        file_name, chr, file_format
):
    """
    Load HiC contact maps from .hic/.mcool or other formats and store into xxx
    :return:
    """
    pass


def strata2horizontal(strata):
    hmat = np.zeros((len(strata[0]), len(strata)))
    for i in range(len(strata)):
        hmat[:len(strata[i]), i] = strata[i]
    return hmat


def strata2vertical(strata):
    vmat = np.zeros((len(strata[0]), len(strata)))
    for i in range(len(strata)):
        vmat[i:, i] = strata[i]
    return vmat


def pick_max_positions(mat, distance_range=(10, 160), line_width=1, window_size=10):
    st, ed = distance_range
    stats = np.sum(mat[:, st:ed], axis=1)
    all_pos = []

    all_peaks, _ = find_peaks(stats, distance=window_size * 2)
    for idx in all_peaks:
        check = enrichment_score_strict(mat, idx, line_width, (st, ed), window_size)
        if np.sum(check) > 0:
            all_pos.append(idx)
    return all_pos


def enrichment_score_strict(mat, idx, line_width=1, distance_range=(20, 40), window_size=10):
    # st, ed = max(distance_range[0], window_size), min(distance_range[1], mat.shape[1] - window_size)
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):
        if j < window_size + half or j >= mat.shape[1] - window_size - half:
            continue
        y = j - distance_range[0]
        line_min = min(np.mean(mat[x1:x2, j-window_size-half:j-half]),
                       np.mean(mat[x1:x2, j+1+half:j+window_size+half+1]))
        neighbor_mean = max(np.mean(mat[idx-window_size:x1, j-window_size-half:j+window_size+half+1]),
                            np.mean(mat[x2+1:idx+window_size+1, j-window_size-half:j+window_size+half+1]))
        new_mat[y] = line_min - neighbor_mean
    return new_mat


def enrichment_score(mat, idx, line_width=1, distance_range=(5, 100), window_size=10):
    # st, ed = max(distance_range[0], window_size), min(distance_range[1], mat.shape[1] - window_size)
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):
        # print(j)
        if j < window_size + half or j >= mat.shape[1] - window_size - half:
            continue
        y = j - distance_range[0]
        # line_min = np.median(np.concatenate(
        #     [mat[x1:x2, j-window_size-half:j-half], mat[x1:x2, j+1+half:j+window_size+half+1]]
        # ))
        line_min = np.median(
            [mat[x1:x2, j - window_size - half:j + window_size + half + 1]]
        )
        neighbor_mean = max(np.mean(mat[idx-window_size:x1, j-window_size-half:j+window_size+half+1]),
                            np.mean(mat[x2+1:idx+window_size+1, j-window_size-half:j+window_size+half+1]))

        ###################
        # NEED A LOWER BOUND FOR EXPECTED VAL!
        ###################
        lower_b = 0
        upper_mlogp = 10
        _exp = max(neighbor_mean, lower_b)

        Poiss = poisson(_exp)
        p_val = 1 - Poiss.cdf(line_min)
        new_mat[y] = min(- np.log10(p_val), upper_mlogp) if p_val > 0 else upper_mlogp
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



