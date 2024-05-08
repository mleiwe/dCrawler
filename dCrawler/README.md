## dCrawler 

This is the Python3 implementation of dCrawler is a fully independent clustering algorithm that only requires a distance threshold (Th(d)) to perform clustering.

Featured in the preprint:  **Automated neuronal reconstruction with super-multicolour fluorescence imaging** , Lewie et,al (2022) [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.20.512984v1)

## Implementation details

A few things are different in the Python version to make it a bit more efficient, the results are the same as the MATLAB code.
- Vectorized operations are used wherever possible to speed up computations.
- The `cKDTree` from the `scipy.spatial` module is used for efficient nearest neighbor search. Instead of calculating distances to all centroids for each point, the k-d tree is used to find the nearest centroid quickly.


![clustering_process](https://github.com/Elsword016/dCrawler/assets/29883365/2f7e6394-50e5-452a-b398-4e3022bf2ce1)

## Installation

```bash
pip install dCrawler
```

## Build from source - local development
Recommended to build a separate environment to prevent any possible errors
- Clone the repository
- Build the package with the command `python setup.py sdist`
- Then `pip install .`

## Usage

```python
from dCrawler import dCrawler
# Initialize the dCrawler object
crawler = dCrawler(threshold=1.0)
crawler.fit(data)
centroids,clusters = crawler.centroids,crawler.clusters
```




## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Cite
```bash
@article {Leiwe2022.10.20.512984,
	author = {Marcus N. Leiwe and Satoshi Fujimoto and Toshikazu Baba and Daichi Moriyasu and Biswanath Saha and Richi Sakaguchi and Shigenori Inagaki and Takeshi Imai},
	title = {Automated neuronal reconstruction with super-multicolour fluorescence imaging},
	elocation-id = {2022.10.20.512984},
	year = {2022},
	doi = {10.1101/2022.10.20.512984},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Fluorescence imaging is widely used for the mesoscopic mapping of neuronal connectivity. However, neurite reconstruction is challenging, especially when neurons are densely labelled. Here we report a strategy for the fully automated reconstruction of densely labelled neuronal circuits. Firstly, we established stochastic {\textquotedblleft}super-multicolour{\textquotedblright} labelling with up to seven different fluorescent proteins using the Tetbow method. With this method, each neuron was labelled with a unique combination of fluorescent proteins, which were then imaged and separated by linear unmixing. We also established an automated neurite reconstruction pipeline based on the quantitative analysis of multiple dyes (QDyeFinder). To classify colour combinations, we used a newly developed unsupervised clustering algorithm, dCrawler, in which data points in multi-dimensional space were clustered based on a given threshold distance. Our new strategy allows for the reconstruction of neurites for up to hundreds of neurons at a millimetre scale without manual tracing.Competing Interest StatementTI, MNL, and SF have filed a patent application for QDyeFinder.},
	URL = {https://www.biorxiv.org/content/early/2022/10/20/2022.10.20.512984},
	eprint = {https://www.biorxiv.org/content/early/2022/10/20/2022.10.20.512984.full.pdf},
	journal = {bioRxiv}
}
```
