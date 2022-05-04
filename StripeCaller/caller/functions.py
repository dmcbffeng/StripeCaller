import sys
sys.path.append("..")
import numpy as np
import scipy.sparse as sp
from scipy.stats import poisson
from scipy import signal
from scipy.ndimage.filters import gaussian_filter1d
from .mat_ops import subsetNpMatrix
from utils.AVL_tree import AVLTree
import math

# _calculated_values = {}
# _poisson_stats = {}


def get_stripe_and_widths(mat, step=1800, sigma=12., rel_height=0.3):
    """
    From the sparse contact map, generate candidate vertical / horizontal stripes.

    Parameters
    ----------
    mat: sp.csr_matrix
        Contact matrix

    step: int
        Step length, i.e., the size of a sub-matrix for finding stripes

    sigma: float


    rel_height: float


    Returns
    ----------
    h_Peaks, v_Peaks: dict
        The dictionaries for horizontal / vertical stripes {location: width}

    """
    v_Peaks = {}
    h_Peaks = {}

    ## peak finder
    for ind in range(step, mat.shape[1] + step, step):
        upper = ind
        lower = ind - step
        # print(lower, upper)
        mat_slice = mat[lower:upper, lower:upper]
        hM, hW, vM, vW = getPeakAndWidths(
            mat_slice, step // 12, sigma=sigma, rel_height=rel_height
        )
        hM += lower
        vM += lower
        for i in range(len(hM)):
            if hM[i] not in h_Peaks.keys():
                h_Peaks[hM[i]] = hW[i]
        for i in range(len(vM)):
            if vM[i] not in v_Peaks.keys():
                v_Peaks[vM[i]] = vW[i]

    for ind in range(step + step // 2, mat.shape[1] + step - step // 2, step):
        upper = ind
        lower = ind - step
        mat_slice = mat[lower:upper, lower:upper]
        # print(lower, upper)
        hM, hW, vM, vW = getPeakAndWidths(mat_slice, step // 12, sigma=sigma, rel_height=rel_height)
        hM += lower
        vM += lower
        for i in range(len(hM)):
            if hM[i] not in h_Peaks.keys():
                h_Peaks[hM[i]] = hW[i]
        for i in range(len(vM)):
            if vM[i] not in v_Peaks.keys():
                v_Peaks[vM[i]] = vW[i]

    return h_Peaks, v_Peaks


def getPeakAndWidths(matrix_in, gap=600, sigma=12., rel_height=0.3):
    HMAT = [matrix_in[i,i:i+gap].mean() for i in range(0, matrix_in.shape[0])]

    hFiltered = gaussian_filter1d(HMAT, sigma=sigma)
    hMax = signal.argrelmax(hFiltered)[0]
    # hMin = signal.argrelmin(hFiltered)[0]
    hwidths = signal.peak_widths(hFiltered, peaks=hMax, rel_height=rel_height)[0]

    VMAT = [matrix_in[:,i].mean() for i in range(0,gap)]
    VMAT += [matrix_in[i-gap:i,i].mean() for i in range(gap,matrix_in.shape[0]-gap)]
    VMAT += [matrix_in[:,i].mean() for i in range(matrix_in.shape[0]-gap,matrix_in.shape[0])]

    vFiltered = gaussian_filter1d(VMAT, sigma=sigma)
    vMax = signal.argrelmax(vFiltered)[0]
    # vMin = signal.argrelmin(vFiltered)[0]
    vwidths = signal.peak_widths(vFiltered, peaks=vMax, rel_height=rel_height)[0]

    return hMax, hwidths, vMax, vwidths


def enrichment_score2(mat, idx, line_width, norm_factors, distance_range=(20, 40), window_size=10,
                      stats_test_log=({}, {})):
    """
    Calculate the enrichment score of a stripe given its location, width and the contact matrix

    Parameters:
    ----------
    mat: np.array (2D)
        Contact matrix generated with strata2horizontal() or strata2vertical()

    idx: int
        The location (index) of the candidate stripe

    line_width: int
        Stripe width (# of bins)

    norm_factors: np.array (1D)
        The vector of normalization factors of the contact map.

    distance_range: tuple
        The distance range (# of bins) for the diagonal for calculating the scores

    window_size: int
        Window size (# of bins)

    stats_test_log: tuple of dict
        Previous log for accelerating statistical tests


    Returns
    ----------
    new_mat: np.array (1D)
        The enrichment score of each pixel along the candidate stripe

    """
    _calculated_values, _poisson_stats = stats_test_log

    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width
    if x1 == x2:
        x2 += 1

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):
        y = j - distance_range[0]
        _min_temp = subsetNpMatrix(mat, (x1, x2), (j - window_size - half, j + window_size + half + 1))
        # if _min_temp == "Empty":
        #     continue
        line_min = np.median([_min_temp])
        # print(_min_temp, line_min)
        _inner_neighbor = subsetNpMatrix(mat, (idx - half - window_size, x1),
                                         (j - window_size - half, j + window_size + half + 1))
        _outer_neighbor = subsetNpMatrix(mat, (x2 + 1, idx + half + window_size + 1),
                                         (j - window_size - half, j + window_size + half + 1))

        if _outer_neighbor.size == 0 or _inner_neighbor.size == 0:
            continue
            
        neighbor_mean = max(np.mean(_inner_neighbor), np.mean(_outer_neighbor))

        # There should be a lower bound for the expected value,
        # otherwise situations like (exp=0.01 and obs=0.02) would also be significant
        # Currently we can set this to 0 until KR norm factors can be loaded
        lower_b = 1 / norm_factors[idx]  # This should be (1 / KR_norm_factors) if we refer to JuiceTools HICCUPS
        _exp = max(neighbor_mean, lower_b)
        _obs = int(line_min)  # the same as floor function when line_min > 0

        # _calculated_values: store all calculated exp-obs pairs in dictionary, in which keys are obs since
        #     they are always integers. Each _calculated_values[obs] is a binary tree for quick searching,
        #     and each tree leaf is a exp value corresponding to the obs value. Since exp values are float,
        #     there is also an integer index attached for searching the exp-obs in dictionary _poisson_stats
        #     (float cannot be dict keys).
        # _poisson_stats: record all calculated result in a dict. It should be
        #     _poisson_stats[(_exp, _obs)] = -log10(p). But _exp is a float and cannot be a dict key, we give
        #     each _exp a unique index and use the index.
        # stats_log: record all p value calculation. Just for benchmarking. Delete this when publishing.

        # global _calculated_values, _poisson_stats  # , stats_log
        tolerance = 0.02

        # check if obs is a value calculated before
        if _obs in _calculated_values:
            # Find the nearest _exp values which were calculated before
            # One larger, one smaller
            (_upper, _lower) = _calculated_values[_obs].search(_exp)
            # If _upper is close enough to _exp, directly use the p value from (_upper-_obs) pair
            if _upper is not None and (_upper.key - _exp) < tolerance * _exp:
                _exp = _upper.key
                _exp_idx = _upper.val  # The integer index for _upper (float cannot be dict keys!)
                mlog_p_val = _poisson_stats[(_exp_idx, _obs)]
            else:
                # Else, calculate p value for _obs-_exp pair and store them in _calculated_values and _poisson_stats
                _exp_idx = _calculated_values[_obs].insert(_exp)  # insert to the binary tree and return an index
                Poiss = poisson(_exp)
                p_val = 1 - Poiss.cdf(_obs)
                if 0 < p_val < 1:
                    mlog_p_val = - np.log10(p_val)
                else:  # Some p values are too small, -log(0) will return an error, so we use -1 to temporarily replace
                    mlog_p_val = -1
                _poisson_stats[(_exp_idx, _obs)] = mlog_p_val
                # stats_log.append([_exp, _obs, mlog_p_val])
        else:  # If _obs is not used before, generate a new binary tree _calculated_values[_obs]
            _calculated_values[_obs] = AVLTree()
            _exp_idx = _calculated_values[_obs].insert(_exp)
            # calculate p value for _obs-_exp pair and store them in _calculated_values and _poisson_stats
            Poiss = poisson(_exp)
            p_val = 1 - Poiss.cdf(_obs)
            if 0 < p_val < 1:
                mlog_p_val = - np.log10(p_val)
            else:  # Some p values are too small, -log(0) will return an error, so we use -1 to temporarily replace
                mlog_p_val = -1
            _poisson_stats[(_exp_idx, _obs)] = mlog_p_val
            # stats_log.append([_exp, _obs, mlog_p_val])

        # Store enrichment score in new_mat
        new_mat[y] = mlog_p_val
    new_mat[new_mat < 0] = np.max(new_mat)  # Replace all "-1"s with the largest -log(p)

    return new_mat


def find_max_slice(arr):
    """
    Given the enrichment score of each pixel along the candidate stripe,
    find the slice with the largest sum as a called stripe.

    Parameters
    ----------
    arr: np.array (1D)
        Enrichment score

    Returns
    ----------
    head: int
        The start location of the slice

    tail: int
        The end location of the slice

    _max: float
        The sum of scores of the slice

    """
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


def phased_max_slice_arr(idx, arr_parallel, width):
    head, tail, _max = find_max_slice(arr_parallel)
    return (idx, head, tail, _max, width)


def merge_positions(lst):#, merge_range):
    """
    Merge stripes that are too close to each other

    Parameters
    ----------
    lst: list
        The stripe list

    Returns
    ----------
    new_lst: list
        The stripe list after merging
    """
    def _merge(small_lst):
        st = min([elm[0] for elm in small_lst])
        ed = max([elm[1] for elm in small_lst])
        head = min([elm[2] for elm in small_lst])
        tail = max([elm[3] for elm in small_lst])
        score = max([elm[4] for elm in small_lst])
        return [st, ed, head, tail, score]

    new_lst = []
    temp = []
    for i, (idx, head, tail, score, width) in enumerate(lst):
        if i == 0:
            temp.append([idx, idx, head, tail, score])
        elif idx - temp[-1][1] <= math.ceil(width):
            temp.append([idx, idx, head, tail, score])
        else:
            new_lst.append(_merge(temp))
            temp = [[idx, idx, head, tail, score]]
    new_lst.append(_merge(temp))
    return new_lst



