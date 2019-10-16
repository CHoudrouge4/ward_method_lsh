from numpy import *
from sklearn import datasets
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import PCA
from sklearn import *

def get_news_group_pca(dimension):
    pca = PCA(n_components=dimension)
    newsgroups_train = fetch_20newsgroups(subset='train')
    vectorizer = TfidfVectorizer(min_df=0.01, max_df=0.95)
    train_data = vectorizer.fit_transform(newsgroups_train.data)
    train_data = train_data.todense()
    pca.fit(train_data)
    Y = pca.transform(train_data)
    print(train_data.shape)
    savetxt("news_1000.in", Y, delimiter=' ')

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

get_news_group(11314)

#generate_sizes()
#get_news_group_pca(2)
