import numpy as np
import scipy.sparse as sp
import scipy.ndimage as ndi

def subsetNpMatrix(matrix, row_bounds, column_bounds):
    """
    Obtain the correct slice of a matrix considering x and y boundaries

    Parameters
    ----------
    matrix: np.array
        The original matrix

    row_bounds: tuple
        The upper and lower bounds for row indices

    column_bounds: tuple
        The upper and lower bounds for column indices

    Returns
    ----------
    subset: np.array
        The sub-matrix after slicing
    """
    rows = np.array([x for x in range(row_bounds[0], row_bounds[1]) if 0 <= int(x) < matrix.shape[0]])
    cols = np.array([y for y in range(column_bounds[0], column_bounds[1]) if 0 <= int(y) < matrix.shape[1]])
    if len(rows)==0 or len(cols)==0:
        return np.empty(0)
    subset = (matrix.ravel()[(cols + (rows * matrix.shape[1]).reshape((-1, 1))).ravel()]).reshape(rows.size, cols.size)
    return subset


def strata2horizontal(strata):
    """
    From strata loaded from load_HiC(), construct the matrix for finding horizontal stripes.
    Assuming the full contact map is N * N and the first p strata are loaded in memory,
    for the k-th row, we pick up the p-pixel slice start from the matrix diagonal (i.e., the k-th element of the row),
    which is row[k: k + p].
    The result is a N * p matrix.

    Parameters
    ----------
    strata: list of np.array (1D)
        Loaded Hi-C strata

    Returns
    ----------
    hmat: np.array (2D)
        The N * p matrix for finding horizontal stripes

    """
    hmat = np.zeros((len(strata[0]), len(strata)))
    for i in range(len(strata)):
        hmat[:len(strata[i]), i] = strata[i]
    return hmat


def strata2vertical(strata):
    """
    From strata loaded from load_HiC(), construct the matrix for finding vertical stripes.
    Assuming the full contact map is N * N and the first p strata are loaded in memory,
    for the k-th column, we pick up the p-pixel slice start from the matrix diagonal (i.e., the k-th element of the column),
    which is column[k: k + p].
    The result is a N * p matrix.

    Parameters
    ----------
    strata: list of np.array (1D)
        Loaded Hi-C strata

    Returns
    ----------
    hmat: np.array (2D)
        The N * p matrix for finding vertical stripes

    """
    vmat = np.zeros((len(strata[0]), len(strata)))
    for i in range(len(strata)):
        vmat[i:, i] = strata[i]
    return vmat


def blank_diagonal_sparse_from_strata(strata, nstrata=0):
    """
    Constructing a full contact map from strata and remove the first few strata near the diagonal

    Parameters
    ----------
    strata: list of np.array (1D)
        Loaded Hi-C strata

    nstrata: int
        Number of strata to remove

    Returns
    ----------
    res: sp.csr_matrix
        The sparse contact map

    Example
    ----------
    # >>> strata = [np.array([1, 2, 3, 4]), np.array([5, 6, 7]), np.array([8, 9])]
    # >>> mat = blank_diagonal_sparse_from_strata(strata, nstrata=1)
    # >>> mat.toarray()
    # [[0. 5. 8. 0.]
    #  [0. 0. 6. 9.]
    #  [0. 0. 0. 7.]
    #  [0. 0. 0. 0.]]
    """
    size = len(strata[0])
    padded_strata = [
        np.concatenate([np.zeros(shape=(i,)), strata[i]]) for i in range(len(strata))
    ]
    for i in range(nstrata):
        padded_strata[i] = np.zeros((size,))
    diags = np.arange(len(strata))
    res = sp.spdiags(padded_strata, diags, size, size, format='csr')
    return res

def power(image, kernel):
    # Normalize images for better comparison.
    image = (image - image.mean()) / image.std()
    return np.sqrt(ndi.convolve(image, np.real(kernel), mode='wrap')**2 +
                   ndi.convolve(image, np.imag(kernel), mode='wrap')**2)

if __name__ == '__main__':
    strata = [np.array([1, 2, 3, 4]), np.array([5, 6, 7]), np.array([8, 9])]
    mat = blank_diagonal_sparse_from_strata(strata, nstrata=1)
    print(mat.toarray())
    # [[0. 5. 8. 0.]
    #  [0. 0. 6. 9.]
    #  [0. 0. 0. 7.]
    #  [0. 0. 0. 0.]]


