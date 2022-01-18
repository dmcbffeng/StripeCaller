import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


current_idx = 0
folder = 'temp'
while os.path.exists(f'{folder}/original_{current_idx}.npy'):
    print(current_idx)
    plt.figure(figsize=(10, 4))
    plt.subplot(121)
    mat1 = np.load(f'{folder}/original_{current_idx}.npy')
    sns.heatmap(mat1, square=False, cmap='coolwarm',
                cbar_kws=dict(use_gridspec=False, location="top"))

    plt.subplot(122)
    mat2 = np.load(f'{folder}/enrichment_{current_idx}.npy').reshape([1, -1])
    sns.heatmap(mat2, square=False, cmap='coolwarm',
                cbar_kws=dict(use_gridspec=False, location="top"))

    # plt.tight_layout()
    plt.savefig(f'{folder}/fig_{current_idx}.png')
    plt.close()

    current_idx += 1



