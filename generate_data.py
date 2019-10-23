from sklearn.datasets.samples_generator import make_blobs
from random import *
from numpy import *

range_from = -100
range_to = 100

# c is the number of cluster
# d is the dimension
def generate_random_clusters(c, d):
    centers = []
    for i in range(c):
        p = []
        for j in range(d):
            p.append(randint(range_from, range_to))
        centers.append(p)
    return centers

#Drawing a chart for our generated dataset
# import matplotlib.pyplot as plt
#set colors for the clusters
# colors = ['r','g','b','c','k','y','m']
# c = []
# for i in y:
#     c.append(colors[i])
# Plot the training points
# plt.scatter(X[:, 0], X[:, 1], c= c)
# plt.gray()
# plt.xlabel('X axis')
# plt.ylabel('Y axis')
# plt.show()

"""
    we can loop over the dimension, then over the number of points
"""
number_of_data = 1
dimensions = {200}
max_n = 20000
nb_center = 20
for d in dimensions:
    n = max_n
    #for n in range(20000, max_n, 1000):
    for i in range(number_of_data):
        centers = generate_random_clusters(nb_center, d)
        X, y = make_blobs(n_samples=n, n_features=d, centers=centers, cluster_std=0.8, center_box=(range_from, range_to), shuffle=True, random_state=0)
        random.shuffle(X)
        file_name = './data' + str(n) + '_' + str(i) + '_' + str(d) + '_' + str(nb_center) + '.in'
        head = str(n) + ' ' + str(d)
        savetxt(file_name, X, delimiter=' ', newline='\n', comments='', header=head)
        head = str(n)
        savetxt(file_name + '_labels', y, delimiter=' ', newline='\n', comments='', header=head)

# nb_center = 50
# centers = generate_random_clusters(nb_center, max_d)
# X, y = make_blobs(n_samples=max_n, n_features=max_d, centers=centers, cluster_std=0.8, center_box=(range_from, range_to), shuffle=True, random_state=0)
# file_name = './data/data' + str(max_n) + '_' + str(max_d) + '_' + str(nb_center) + '.in'
# n , d = X.shape
# random.shuffle(X)
# head = str(n) + ' ' + str(d) # ' ' + str(nb_center)
# savetxt(file_name, X, delimiter=' ', newline='\n', comments='', header=head)
