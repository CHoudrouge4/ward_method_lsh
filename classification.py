### for example on digits
import read_tree
from sklearn.metrics.classification import accuracy_score, log_loss
from sklearn.datasets import load_digits
from sklearn.preprocessing import scale
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn import datasets

iris = datasets.load_iris()
data = scale(iris.data)
labels = iris.target

tree = read_file('output.in');
result_of_alg = clusters(tree, k)
def convert(clusters, n):
    clustering_vect = [0]*n
    for i in range(len(clusters)):
        for p in clusters[i]:
            clustering_vect[p] = i
    return clustering_vect



#test = [[3,4,5], [1,8,0], [2,6,7]]
#print(convert(test, 9))

normalized_mutual_info_score(convert(result_of_alg, len(labels)), labels)
