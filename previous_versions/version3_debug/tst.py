import numpy as np


def subsetNpMatrix(matrix, row_bounds, column_bounds):
    rows = np.array([x for x in range(row_bounds[0], row_bounds[1]) if 0 <= int(x) < matrix.shape[0]])
    cols = np.array([y for y in range(column_bounds[0], column_bounds[1]) if 0 <= int(y) < matrix.shape[1]])
    #     if rows.size==0 or cols.size==0:
    #         print("--------------")
    #         print(f"rows: {rows}")
    #         print(f"row_bounds: {row_bounds}")
    #         print(f"cols: {cols}")
    #         print(f"column_bounds: {column_bounds}")
    #         return "Empty"
    #     if row_bounds[0]>row_bounds[1]:
    #         print(f"rowbounds0 > rowbounds1, {row_bounds[0]}, {row_bounds[1]}")
    subset = (matrix.ravel()[(cols + (rows * matrix.shape[1]).reshape((-1, 1))).ravel()]).reshape(rows.size, cols.size)
    return subset


x = np.random.normal(0, 1, (100, 100))


sub = x[-1:5, 96:101]
print(sub.shape)

sub = subsetNpMatrix(x, (-1, 5), (96, 101))
print(sub.shape)

