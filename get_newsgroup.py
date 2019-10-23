from numpy import *
from sklearn import datasets
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import PCA
from sklearn import *
from sklearn import preprocessing

def get_news_group_pca(size, dimension):
    pca = PCA(n_components=dimension)
    newsgroups_train = fetch_20newsgroups(subset='train')
    vectorizer = TfidfVectorizer(min_df=0.01, max_df=0.95)
    train_data = vectorizer.fit_transform(newsgroups_train.data)
    train_data = train_data.todense()
    pca.fit(train_data[:size, :])
    Y = pca.transform(train_data[:size, :])
    S = 10000.0 * Y
    print(train_data.shape)
    size, d = S.shape
    savetxt('news_' + str(size) + '_' + str(dimension) +'.in', S[:size, :], delimiter=' ', comments='',  header=str(size) + ' ' + str(d))

def get_news_group(size):
    newsgroups_train = fetch_20newsgroups(subset='train')
    vectorizer = TfidfVectorizer(min_df=0.01, max_df=0.95)
    train_data = vectorizer.fit_transform(newsgroups_train.data)
    train_data = train_data.todense()
    print(train_data.shape)
    Y = train_data[0:size, :]
    _, d = Y.shape
    savetxt("news_" + str(size) + ".in", Y, delimiter=' ', comments='',  header=str(size) + ' ' + str(d))


def generate_sizes():
    sizes = 100
    while sizes < 11314:
        get_news_group(sizes)
        sizes = sizes * 2

#get_news_group(5000)

#generate_sizes()
get_news_group_pca(11314, 10)
