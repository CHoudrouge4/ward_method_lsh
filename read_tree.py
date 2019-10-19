from sklearn import datasets
from sklearn.cluster import AgglomerativeClustering
from numpy import *
import time
from numpy import *

import sys
sys.setrecursionlimit(12000)

def get_news_group_pca(dimension):
    from sklearn.datasets import fetch_20newsgroups
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.decomposition import PCA

    newsgroups_train = fetch_20newsgroups(subset='train')
    vectorizer = TfidfVectorizer(min_df=0.01, max_df=0.95)
    train_data = vectorizer.fit_transform(newsgroups_train.data)
    train_data = train_data.todense()
    train_labels = newsgroups_train.target;

    # initialise PCA with n_components = dimension
    pca = PCA(n_components=dimension)
    # apply pca
    pca.fit(train_data)
    Y = pca.transform(train_data)
    n_elements = len(unique(train_labels))
    return Y, n_elements, train_labels, len(set(train_labels))

def get_news_group(size):
    from sklearn.datasets import fetch_20newsgroups
    from sklearn.feature_extraction.text import TfidfVectorizer
    newsgroups_train = fetch_20newsgroups(subset='train')
    vectorizer = TfidfVectorizer(min_df=0.01, max_df=0.95)
    train_data = vectorizer.fit_transform(newsgroups_train.data)
    train_data = train_data.todense()
    train_labels = newsgroups_train.target;
    # initialise PCA with n_components = dimension
    n_elements = len(unique(train_labels[0:size]))
    return train_data[0:size, :], n_elements, train_labels[0:size], len(set(train_labels[0:size]))

def read_file(filename):
    f = open(filename, "r")
    lines = f.readlines()
    n = len(lines)+1
    nb_clust = n
    #print(n)
    clusters = {3*i*n+1: i for i in range(n)}
    T = [[i,-1] for i in range(n)]
    for l in lines:
        words = l.split(";")
        res, c1, c2 = words[:3]
        idres_str = res.split(",")
        idres = int(idres_str[0])*3*n + int(idres_str[1])

        idc1_str = c1.split(",")
        idc1 = int(idc1_str[0])*3*n + int(idc1_str[1])

        idc2_str = c2.split(",")
        idc2 = int(idc2_str[0])*3*n + int(idc2_str[1])

        #print(idres_str, idc2_str, idc1_str)
        #print(idres, idc1, idc2)

        clusters[idres] = nb_clust
        T.append([clusters[idc1], clusters[idc2]])
        nb_clust+=1

    return T


def get_cluster(tree, index):
    if tree[index][0] == index:
        return [index]
    return get_cluster(tree, tree[index][0]) + get_cluster(tree, tree[index][1])

def clusters(tree, k):
    nums = []
    i = 1
    while len(nums) != k:
        current = len(tree)-i
        if current in nums:
            nums.remove(current)
        c1 = tree[-i][0]
        c2 = tree[-i][1]
        nums.append(c1)
        if c1 != current:
            nums.append(c2)
        i+=1

    return [get_cluster(tree, nums[i]) for i in range(k)]

def get_dataset(name):
    from sklearn.preprocessing import scale
    data = []
    if name == "cancer":
        from sklearn.datasets import load_breast_cancer
        dataset = load_breast_cancer()
    elif name == "digits":
        from sklearn.datasets import load_digits
        dataset = load_digits()
    elif name == "iris":
        from sklearn.datasets import load_iris
        dataset = load_iris()
    elif name == "boston":
        from sklearn.datasets import load_boston
        dataset = load_boston()
    elif name == "KDD":
        from sklearn.datasets import fetch_kddcup99
        dataset = fetch_kddcup99(subset='SF')
        data = dataset.data[:2000, [0,2,3]]
    else:
        print("Unknown name of dataset")
        exit(-1)


    labels = dataset.target
    if data == []:
        data = scale(dataset.data)
        n_samples, n_features = data.shape
        n_elements = len(unique(labels))
    return data, n_elements, labels, len(set(labels))

# example
# t = [[0, -2], [1, -1], [2, -1], [3, -1], [4, -1], [0, 1], [2, 3], [5, 6], [7, 4]]
# print(clusters(t, 3))
#from sklearn.datasets import load_iris
from sklearn.metrics.cluster import normalized_mutual_info_score

def convert(clusters, n):
    clustering_vect = [0]*n
    for i in range(len(clusters)):
        for p in clusters[i]:
            clustering_vect[p] = i
    return clustering_vect


def approx_vs_ward(e, number_of_visited_leafs, numebr_of_trees, name, dimension, size):
#print("epsilon ", "0." + str(epsilon), numebr_of_trees, number_of_visited_leafs)
# #eps = [25, 50, 75, 85, 95, 200, 400];

    output_file = 'result_epsilons.txt'
    with open(output_file, 'a') as file:
        file.write('epsilon 0.' + str(e) + ' ' + str(numebr_of_trees) + ' ' + str(number_of_visited_leafs) + ' ' + str(size) + ' ' + str(dimension) +'\n')
        data, n, labels, k = get_news_group(size)
        file_name = './' + name + '_' + str(dimension) + '_' + str(e) + '_' + str(numebr_of_trees) + '_' + str(number_of_visited_leafs) + ".out"
        T = read_file(file_name)
        print (k)
        clust = clusters(T, k)
        print (len(clust))
        file.write('Algo ' + str(normalized_mutual_info_score(convert(clust, len(labels)), labels)) + '\n')
        ward = AgglomerativeClustering(n_clusters=k, linkage='ward', connectivity=None)
        data = loadtxt('news_' + str(size) + '.in', skiprows = 1)
        print (data.shape)
        clustering = ward.fit(data)
        clust = clustering.labels_
        file.write('std_ward ' + str(normalized_mutual_info_score(clust, labels)) + '\n')



def app_lsh_ward(data_file, input_file, size):
    data, n, labels, k = get_news_group(size)
    T = read_file(input_file)
    print (k)
    clust = clusters(T, k)
    print (len(clust))
    print('Algo ' + str(normalized_mutual_info_score(convert(clust, len(labels)), labels)) + '\n')
    ward = AgglomerativeClustering(n_clusters=k, linkage='ward', connectivity=None)
    data = loadtxt(data_file, skiprows = 1)
    print (data.shape)
    clustering = ward.fit(data)
    clust = clustering.labels_
    print('std_ward ' + str(normalized_mutual_info_score(clust, labels)) + '\n')

#
#         for name in data_sets:
#             data, n, labels, k = get_dataset(name)
#             ward = AgglomerativeClustering(n_clusters=k, linkage='ward', connectivity=None)
#             clustering = ward.fit(data)
#             clust = clustering.labels_
#             file.write('Ward ' + str(normalized_mutual_info_score(clust, labels)) + '\n')
#
#         for name in data_sets:
#             data, n, labels, k = get_dataset(name)
#             ward = AgglomerativeClustering(n_clusters=k, linkage='average', connectivity=None)
#             clustering = ward.fit(data)
#             clust = clustering.labels_
#             file.write('Average ' + str(normalized_mutual_info_score(clust, labels)) + '\n')
#
#         for name in data_sets:
#             data, n, labels, k = get_dataset(name)
#             ward = AgglomerativeClustering(n_clusters=k, linkage='single', connectivity=None)
#             clustering = ward.fit(data)
#             clust = clustering.labels_
#             file.write('Single ' + str(normalized_mutual_info_score(clust, labels)) + '\n')



def readFILE(file_name):
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



  #
  # std::vector<int> trees = {4 , 16};
  # std::vector<int> leaves = {5, 128};
  # std::vector<float> epsilons = {0.5, 1, 7};

#trees = [2]
#leaves = [10]
#epsilons = [800]
#d = 20
#k = 10
#ward = AgglomerativeClustering(n_clusters=k, linkage='ward', connectivity=None)
# with open('ward_accuracy7.txt', 'w') as f:
#     for e in epsilons:
#         for t in trees:
#             for l in leaves:
#                 for i in range(10000, 20000, 1000):
#                     for j in range(1):
#                         if i == 18000:
#                             j = 1
#                         data_file = './data/data' + str(i) + '_' + str(j) + '_' + str(d) + '_' + str(k) + '.in'
#                         data = readFILE(data_file)
#                         start = time.time()
#                         clustering = ward.fit(data)
#                         end = time.time()
#                         labels = clustering.labels_
#                         res_file = 'data' + str(i) + '_' + str(j) + '_' + str(d) + '_' + str(k) + '_' + str(e) + '_' + str(t) + '_' + str(l) + '.out'
#                         T = read_file(res_file)
#                         clust = clusters(T, k);
#                         acc = normalized_mutual_info_score(convert(clust, len(labels)), labels)
#                         print(str(end - start))
#                         f.write(str(end - start) + ' ' + str(acc) + ' ')
#i                     f.write('\n')

#
#data = readFILE("./data/iris.in")
#clustering = ward.fit(data)
#random.shuffle(data)
#labels = clustering.labels_
# data, n, labels, k = get_dataset('boston')
# file_name = './data/boston800_2_10.out'
# T = read_file(file_name)
# clust = clusters(T, k)
# print(normalized_mutual_info_score(convert(clust, len(labels)), labels))
#data, n, labels, k = get_dataset('boston')
#print(k)
## e psilon, number_of_visited_leafs, number_of_trees, dimension

def exp_sizes_acc():
    sizes = 100
    while sizes < 11314:
        approx_vs_ward(100, 128, 25, 'news' + str(sizes), 2164, sizes)
        sizes = sizes * 2

def tree_sizes_perf():
    sizes = 11314
    tree = [1, 4, 8, 12, 16, 20, 24, 28, 32]
    for t in tree:
        approx_vs_ward(200, 128, t, 'news' + str(sizes), 2164, sizes)

def leaves_perf():
    sizes = 11314
    leaves = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    for l in leaves:
        approx_vs_ward(200, l, 16, 'news' + str(sizes), 2164, sizes)

def epsilons_perf():
    sizes = 11314
    epsilons = [50, 100, 200, 400, 800, 1000]
    for e in epsilons:
        approx_vs_ward(e, 32, 16, 'news' + str(sizes), 2164, sizes)

#epsilons_perf()
#leaves_perf()

#tree_sizes_perf()
#_, _, labels, _ = get_news_group(2)
#for u in labels:
#    print (u)
def simple_exp():
    data_name = ['iris', 'cancer', 'boston', 'digits']
    for name in data_name:
        print(name)
        data, n, labels, k = get_dataset(name)
        ward = AgglomerativeClustering(n_clusters=k, linkage='ward', connectivity=None)
        clustering = ward.fit(data)
        clust = clustering.labels_
        print('Ward ' + str(normalized_mutual_info_score(clust, labels)) + '\n')
        T = read_file(name + '.out')
        clust = clusters(T, k)
        print('approx' + str(normalized_mutual_info_score(convert(clust, len(labels)), labels)))

#app_lsh_ward('news_5000_100.in', 'news_5000_100.out', 5000)
simple_exp()
