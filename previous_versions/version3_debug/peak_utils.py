import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy import signal
from scipy.ndimage.filters import gaussian_filter1d
from multiprocessing import Process, Manager

def blank_diagonal2(matr, strata = 50):
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
    lil=csr.tolil()
    for i in range(lil.shape[0]):
        lil[i:i+strata+1,i:i+strata+1] = 0
    csr = lil.tocsr()
    return csr

def getPeakAndWidths(matrix_in, gap=600, sigma=12, rel_height=0.3):
#    
    HMAT = [matrix_in[i,i:i+gap].mean() for i in range(0,matrix_in.shape[0])]
#
    hFiltered = gaussian_filter1d(HMAT, sigma=sigma)
    hMax = signal.argrelmax(hFiltered)[0]
    hMin = signal.argrelmin(hFiltered)[0]
    hwidths = signal.peak_widths(hFiltered, peaks=hMax, rel_height=rel_height)[0]
#
    VMAT = [matrix_in[:,i].mean() for i in range(0,gap)]
    VMAT += [matrix_in[i-gap:i,i].mean() for i in range(gap,matrix_in.shape[0]-gap)]
    VMAT += [matrix_in[:,i].mean() for i in range(matrix_in.shape[0]-gap,matrix_in.shape[0])]
#
    vFiltered = gaussian_filter1d(VMAT, sigma=sigma)
    vMax = signal.argrelmax(vFiltered)[0]
    vMin = signal.argrelmin(vFiltered)[0]
    vwidths = signal.peak_widths(vFiltered, peaks=vMax, rel_height=rel_height)[0]
#    
    return hMax, hwidths, vMax, vwidths


def matrix2peaks(mat,stride=2500,step=600):
    hM_res = list()
    hW_res = list()
    vM_res = list()
    vW_res = list()
    for i in range(1,mat.shape[0],stride):
        lower = i-1
        upper = i+stride-1
        if upper > mat.shape[0]:
            upper = mat.shape[0]
        Z_tester = mat[lower:upper,lower:upper]
        hM,hW,vM,vW = getPeakAndWidths(Z_tester, step)
        hM += lower
        vM += lower
        hM_res.extend(hM)
        vM_res.extend(vM)
        hW_res.extend(hW)
        vW_res.extend(vW)
    return hM_res, hW_res, vM_res, vW_res

