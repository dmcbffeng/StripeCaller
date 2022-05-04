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


def blank_diagonal2(matr, strata = False):
    """
    in: edgelist, strata (n entries off main diagonal to zero)
    out:matrix with n blanked diagonal
    """
#     int_shape = (int(max(matr[:,1])+1), int(max(matr[:,1])+1))
#     coo = sp.coo_matrix((matr[:, 2], (matr[:, 0], matr[:, 1])), shape=int_shape, dtype=matr.dtype)
#     csr = coo.tocsr()
#     csr_org = csr.copy()
    coo = sp.coo_matrix(matr)
    csr = coo.tocsr()
    lil = csr.tolil()
    for i in range(lil.shape[0]):
        lil[i:i+strata+1,i:i+strata+1] = 0
    csr = lil.tocsr()
    return csr


def blank_diagonal(mat, nstrata=0):
    return np.triu(mat, nstrata) if nstrata > 0 else mat


def strata2triu(strata):
    mat = np.zeros((len(strata[0]), len(strata[0])))
    for i in range(len(strata)):
        for j in range(len(strata[i])):
            mat[j, j + i] = strata[i][j]
    return mat


