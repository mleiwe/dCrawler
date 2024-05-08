## dCrawler 

This is the Python3 implementation of dCrawler is a fully independent clustering algorithm that only requires a distance threshold (Th(d)) to perform clustering.

## Installation

```bash
pip install dCrawler
```

## Usage

```python
from dCrawler import dCrawler
# Initialize the dCrawler object
crawler = dCrawler(threshold=1.0)
crawler.fit(data)
centroids,clusters = crawler.centroids,crawler.clusters
```

## Animation


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