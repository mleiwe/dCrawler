# dCrawler
dCrawler is a fully independent clustering algorithm that only requires a distance threshold (Th(d)) to perform clustering.
## MATLAB Implementation
### Requirements
A MATLAB version of the algortihm is provided.
We have tested this on MATLAB versions 2020 and above and no additional libraries are required

### Installation
Once the codes have been downloaded onto your PC simply add the file paths to your MATLAB environment (home tab, environment section, click on the button "Set Path") and then you should be good to go!

### Use
There are two versions of dCrawler provided here.
* `mnl_Basic_dCrawler` - This is the simplest purest form of the crawler. Simply insert your values as the "InputMatrix", specify your threshold "EuThresh", and whether you want to have output figures ("FigYN").
* `mnl_Weighted_dCrawler` - An amended version where each point is given a particular weight when calculating the centroid positions.

## Python
### Requirements
We currently have a preliminary version in the form of a jupyter notebook. This was created in collaboration with [Biswanath Saha](https://github.com/Elsword016). There is an `environment.yaml` file at the moment but for more details see the `Readme.md` within the python implementation folder

### Deployment
To come...
* pip install
* Docker implementation
* Web-based deployment

## Want to know more?
Please see the the pdf [Schema_For_dCrawler](https://github.com/mleiwe/dCrawler/blob/main/Schema_For_dCrawler.pdf). Or watch the [attached movie](https://github.com/mleiwe/dCrawler/blob/main/SupplementaryVideo1_dCrawlerDemo.avi).

### When you should you use dCrawler?
There's no one perfect clustering algorithm (or at least to date!). Usually, it depends on the distribution of your points. In our examples below we have a comparison between DBSCAN and dCrawler. 

Top row (A) shows a situation where dCrawler is advantageous. I.e. if two clusters are nearby or potentially overlapping DBSCAN will allocate them as a single cluster, whereas because dCrawler's distance is derived from the putative cluster thresholds and applies the additional steps of adjustment and merge it can derive the position of the two clusters.

The middle row (B) shows where DBSCAN performs better; where clusters are in non-gaussian distributions (e.g. the smile). dCrawler fails in this situation as it splits the smile into three. If the threshold (d) was increased then the eyes would be included too.

Typically the situation is usually seen where the data is often continuous, such as sorting colours such as in the `peppers.png` image provided by MATLAB (C). DBSCAN will assign all the pixels to the same cluster, while dCrawler can identify 9 different colours.

![FigSX_ClusterAlgorithmComparison](https://github.com/mleiwe/dCrawler/assets/29621219/9ba7aa6e-cfd1-425c-b0e7-0fd3c800e96f)

## Any more questions?
Please contact me (mleiwe), you can find my contact information here --> https://github.com/mleiwe.
