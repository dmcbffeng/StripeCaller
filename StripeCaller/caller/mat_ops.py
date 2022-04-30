import numpy as np
import scipy.sparse as sp


def subsetNpMatrix(matrix, row_bounds, column_bounds):
    rows = np.array([x for x in range(row_bounds[0], row_bounds[1]) if 0 <= int(x) < matrix.shape[0]])
    cols = np.array([y for y in range(column_bounds[0], column_bounds[1]) if 0 <= int(y) < matrix.shape[1]])
    if len(rows)==0 or len(cols)==0:
        return np.empty(0)
    subset = (matrix.ravel()[(cols + (rows * matrix.shape[1]).reshape((-1, 1))).ravel()]).reshape(rows.size, cols.size)
    return subset


def strata2triu(strata):
    mat = np.zeros((len(strata[0]), len(strata[0])))
    for i in range(len(strata)):
        for j in range(len(strata[i])):
            mat[j, j + i] = strata[i][j]
    return mat


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


def blank_diagonal_sparse_from_strata(strata, nstrata=0):
    """
    # >>> strata = [np.array([1, 2, 3, 4]), np.array([5, 6, 7]), np.array([8, 9])]
    # >>> blank_diagonal_sparse_from_strata(strata, nstrata=1)
    """
    size = len(strata[0])
    padded_strata = [
        np.concatenate([np.zeros(shape=(i,)), strata[i]]) for i in range(len(strata))
    ]
    for i in range(nstrata):
        padded_strata[i] = np.zeros((size,))
    diags = np.arange(len(strata))
    return sp.spdiags(padded_strata, diags, size, size, format='csr')


if __name__ == '__main__':
    strata = [np.array([1, 2, 3, 4]), np.array([5, 6, 7]), np.array([8, 9])]
    mat = blank_diagonal_sparse_from_strata(strata, nstrata=1)
    print(mat.toarray())
    # [[0. 5. 8. 0.]
    #  [0. 0. 6. 9.]
    #  [0. 0. 0. 7.]
    #  [0. 0. 0. 0.]]


