import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


current_idx = 0
folder = 'temp'
while os.path.exists(f'{folder}/original_{current_idx}.npy'):
    if current_idx >= 500:
        break

    print(current_idx)

    mat3 = np.load(f'{folder}/head_tail_{current_idx}.npy')
    st, ed = mat3[0], mat3[1]
    if ed - st < 50:
        current_idx += 1
        continue

    plt.figure(figsize=(10, 6))
    plt.subplot(211)
    mat1 = np.load(f'{folder}/original_{current_idx}.npy')
    mat1 = np.log10(mat1 + 1)
    sns.heatmap(mat1, square=False, cmap='coolwarm',
                cbar_kws=dict(use_gridspec=False, location="top"))

    plt.subplot(212)
    mat2 = np.load(f'{folder}/enrichment_{current_idx}.npy')  # .reshape([1, -1])
    # sns.heatmap(mat2, square=False, cmap='coolwarm',
    #             cbar_kws=dict(use_gridspec=False, location="top"))
    plt.plot(np.arange(len(mat2)), mat2)
    plt.xlim([0, len(mat2)])
    thr = 0.15
    plt.ylim([0, max(max(mat2), -np.log10(thr)) + 0.6])
    plt.axhline(y=-np.log10(thr), color='gray')
    plt.ylabel('-log10(P)')
    plt.axhline(y=-np.log10(thr), color='red', lw=3, xmin=st / len(mat2), xmax=ed / len(mat2))
    print(st, ed, max(-np.log10(thr), max(mat2)) + 0.5)

    # plt.tight_layout()
    plt.savefig(f'temp_fig/fig_{current_idx}.png')
    plt.close()

    current_idx += 1





