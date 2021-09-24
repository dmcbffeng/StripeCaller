# StripeCaller

### Summary



### Installation
  **Required Python Packages**
  - Python (version >= 3.6)
  - numpy (version >= 1.15.4)
  - scipy (version >= 1.0)
  - simplejson
  - six
  - h5py

  **Install from GitHub**

  You can install the package with following command:

  ```console
    $ git clone https://github.com/liu-bioinfo-lab/(to be updated)
    $ cd StripeCaller
    $ python setup.py install
  ```

  **Install from PyPI**

  ```console
    $ pip install (to be updated)
  ```


### Usage
   **Supported Formats**
  - .hic / .mcool / .pairs / .pairs.gz

  **Import Package**
  ```console
  >>>import StripeCaller
  ```

  **Call Stripes**

  ```console
  >>>(to be updated)
  ```
  - hic_file (str): file path
  - reference_genome (str): reference genome
  - chromosomes (list): which chromosomes to calculate
  - output_file (str): output bedpe path
  - norm (str): recommend: "balanced", can also be "none"
  - threshold (float): p value threshold
  - max_range (int): max distance off the diagonal to be calculated
  - resolution (int): resolution
  - min_length (int): minimum length of stripes
  - min_distance (int): threshold for removing stripes too far away from the diagonal
  - stripe_width (int): stripe width (# of bins at the given resolution)
  - merge (int): merge stripes which are close to each other (# of bins)
  - window_size (int): size of the window for calculating enrichment score

  **Recommended Parameter settings**
  - HiC:
  
  - Micro-C:



### Citation
(to be updated)

### References

Durand, Neva C., et al. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell systems 3.1 (2016): 99-101.

Abdennur, Nezar, and Leonid A. Mirny. "Cooler: scalable storage for Hi-C data and other genomically labeled arrays." Bioinformatics 36.1 (2020): 311-316.

Vian, Laura, et al. "The energetics and physiological impact of cohesin extrusion." Cell 173.5 (2018): 1165-1178.

