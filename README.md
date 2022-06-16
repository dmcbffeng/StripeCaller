# StripeCaller - Quagga

### Summary
The python package for calling stripes from Hi-C/Micro-C contact maps.

Stripes are architectural features of the genome that represent asymmetric extrusions of CTCF and cohesion, thought to play regulatory roles in a cell's development or state.
Few established methods for automatic stripe identification exist, and outstanding methods cannot do so in an unbiased way.
To this end, we developed Quagga, a tool for detection and statistical verification of genomic architectural stripes from Hi-C or Micro-C chromatin contact maps.


### Installation
  **Required Python Packages**
  - Python (version >= 3.6)
  - numpy (version >= 1.15.4)
  - scipy (version >= 1.0)
  - scikit-image
  - pandas (version >= 1.0.0)
  - simplejson
  - six
  - h5py
  - requests

  **Install from GitHub**

  You can install the package with following command:

  ```console
    $ git clone https://github.com/dmcbffeng/StripeCaller.git
    $ cd StripeCaller
    $ python setup.py install
  ```


### Usage
  **Supported Formats**
  - .hic / .mcool / .pairs / .pairs.gz

  **Import Package**
  ```console
  >>>from Quagga import stripe_caller
  ```

  **Call Stripes**

  ```console
  >>>stripe_caller(
    hic_file='GM12878.hic',
    reference_genome='hg38',
    chromosomes=['chr2'],
    output_file='GM12878_chr2.bedpe',
    norm='balanced',
    threshold=0.15,
    resolution=5000,
    max_range=2000000,
    min_length=200000,
    min_distance=200000,
    max_width=25000,
    window_size=35,
    centromere_file='removed_regions.bed',
    N_threads=1,
    nstrata_blank=10,
    sigma=3,
    rel_height=0.3,
    gabor_freq=0.1
  )
  ```
   **Call stripes through command line interface**
   ```console
   stripe_caller --hic /path/to/microC/file/HFFc6.hic --output HFFc6_chr1.bedpe --chr 'chr1' --rg 'hg38' --max_range 2000000\
   --resolution 5000 --min_length 200000 --min_distance 500000 --window_size 35 --sigma 12 --N_cores 26  --rel_height 0.3\
   --norm balanced --thr 0.15 --centromere_file removed_regions.bed --nstrata_blank 10
   
   ```
  
  
  - hic_file (str): file path
  - reference_genome (str): reference genome
  - chromosomes (list): which chromosomes to calculate
  - output_file (str): output bedpe path
  - norm (str): Hi-C normalization method. Recommend: "balanced", can also be "none"
  - threshold (float): p value threshold
  - resolution (int): resolution
  - max_range (int): max distance off the diagonal to be calculated
  - min_length (int): minimum length of stripes
  - min_distance (int): threshold for removing stripes too far away from the diagonal
  - stripe_width (int): stripe width (# of bins at the given resolution)
  - max_width (int): maximum width of stripes
  - window_size (int): size of the window for calculating enrichment score
  - centromere_file (str): the bed file for centromere regions in which will be removed during calculation (can be None)
  - N_threads (int): the number of threads
  - nstrata_blank (int): cells from main diagonal of contact matrix to set to zero
  - sigma (float): Threshold of peak detection of file-summation spectrogram
  - rel_height (float): Relative proportional height from peak of file-summation spectrogram to measure peak's width
  - gabor_freq (float): Frequency for the gabor kernel applied to matrix image, directionally filters for lines  

  **Recommended Parameter settings**
  - HiC: nstrata_blank=1, norm="balanced", threshold=0.15, max_range=2000000, resolution=5000, min_length=300000, min_distance=2, max_width=25000, window_size=35, sigma=1, rel_height=0.3
  - Micro-C: nstrata_blank=50, norm="balanced", threshold=0.15, max_range=2000000, resolution=1000, min_length=300000, min_distance=2, max_width=5000, window_size=600, sigma=12, rel_height=0.3


### Output

The output is a bedpe file in which the first columns 1 to 3 correspond to x axis and columns 2 to 4 correspond to y axis.
They are both in the order of "chromosome - start coordinate - end coordinate".
These two ranges define a rectangle region.
We do not cluster them into "horizontal" or "vertical" in the output and the direction can be inferred by comparing the lengths of two ranges.

The bedpe file can be loaded by JuiceBox (https://aidenlab.gitbook.io/juicebox/desktop).
On the top left corner, by choosing "View - Show annotation panel - 2D annotations - Add local",
users can import the called stripes to a pre-loaded Hi-C contact map.
Here is an example region:

![GitHub Logo](/example/JuiceBox_example_region.png)


### Citation
(to be updated)

### References

Durand, Neva C., et al. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell systems 3.1 (2016): 99-101.

Abdennur, Nezar, and Leonid A. Mirny. "Cooler: scalable storage for Hi-C data and other genomically labeled arrays." Bioinformatics 36.1 (2020): 311-316.

Vian, Laura, et al. "The energetics and physiological impact of cohesin extrusion." Cell 173.5 (2018): 1165-1178.

