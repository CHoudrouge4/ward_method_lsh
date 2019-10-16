from sklearn import datasets
from numpy import *

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
    elif name == "newsgroup":
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn.datasets import fetch_20newsgroups
        dataset = fetch_20newsgroups(subset='train')
        vectorizer = TfidfVectorizer()
        data = vectorizer.fit_transform(dataset.data)
    #    data = vectors.nnz / float(vectors.shape[0])
        labels = dataset.target
        return data.toarray(), 1, labels, len(set(labels))
    else:
        print("Unknown name of dataset")
        exit(-1)

    labels = dataset.target
    if data == []:
        data = scale(dataset.data)
    n_samples, n_features = data.shape
    #n_elements = len(unique(labels))
    return data, 1, labels, len(set(labels))

#data_name = ["iris", "cancer", "digits", "boston", "KDD"]

data_name = ["iris"]
for name in data_name:
    print(name)
    x, _, _, _ = get_dataset(name)
    n, d = x.shape
    print (name, " shape: ",  n, ' ', d)
    file_name = name + ".in"
    savetxt(file_name, x, delimiter=' ')
