from dCrawler_api import *

import pandas as pd
data_path = "../DemoData/demodata.csv"
df = pd.read_csv(data_path)
data = df.to_numpy()[:,:2]

crawler = dCrawler(threshold=1.0)
crawler.fit(data)
centroids,clusters = crawler.centroids,crawler.clusters
print(clusters)


