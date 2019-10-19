from sklearn.cluster import AgglomerativeClustering
import numpy as np
from numpy.random import *
import time
from datetime import datetime
from numpy import *

#X = rand(1797, 64)
def read_file(file_name):
    with open(file_name) as f:
        content = f.readline()
        content = content.split(' ')
        n = int(content[0])
        d = int(content[1])
        #k = int(content[2])
        X = zeros((n, d))
        for i in range(n):
            content = f.readline()
            content = content.split(' ')
            for j in range(d):
                X[i, j] = float(content[j])
        return X
#X = rand(15000, 2)
X = read_file("./news_11314_10.in")
print (X.shape)
n_clusters = 1
ward = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward', connectivity=None)
start = time.time()
ward.fit(X)
end = time.time()
print(end - start)
