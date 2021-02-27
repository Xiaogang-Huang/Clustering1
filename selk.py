import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
import pandas as pd


def centroid_init(data, k):
    n = data.shape[0]
    index = np.random.randint(0, n, k)
    centroids = data[index, :]
    return centroids


def centroid_pairwise_dist(X, centroids):
    return pairwise_distances(X, centroids, metric="euclidean")


def distEclud(vecA, vecB):
    # vecA,vecB是数组形式，列表形式不行
    return sum((vecA-vecB)**2)**0.5


def revise_centroids(data, k, cluster_assignment):
    new_centroids = []
    for i in range(k):
        # Select all data points that belong to cluster i. Fill in the blank (RHS only)
        member_data_points = data[cluster_assignment == i]
        # Compute the mean of the data points. Fill in the blank (RHS only)
        centroid = member_data_points.mean(axis=0)
        new_centroids.append(centroid)
    new_centroids = np.array(new_centroids)

    return new_centroids


def selk(data, k):
    n = data.shape[0]
    centroids = centroid_init(data, k)
    bound_low = centroid_pairwise_dist(data, centroids)
    cluster_assignment = np.argmin(bound_low, axis=1)
    bound_up = np.min(bound_low, axis=1)
    clusterChanged = True
    r = [True]*n
    t = 0
    while clusterChanged:
        d = centroid_pairwise_dist(centroids, centroids)
        row, col = np.diag_indices_from(d)
        d[row, col] = float("inf")
        s = 0.5*np.min(d, axis=1)
        for i in range(n):
            c = cluster_assignment[i]
            if bound_up[i] > s[c]:
                for j in range(k):
                    if j != c and bound_up[i] > bound_low[i, j] and\
                      bound_up[i] > 0.5*d[c, j]:
                        if r[i]:
                            temp = distEclud(data[i, :], centroids[c, :])
                            bound_up[i] = temp
                            r[i] = False
                        else:
                            temp = bound_up[i]
                        if temp > bound_low[i, j]:
                            bound_low[i, j] = distEclud(data[i, :],
                                                        centroids[j, :])
                            if bound_low[i, j] < temp:
                                cluster_assignment[i] = j
            new_centroids = revise_centroids(data, k, cluster_assignment)
            delta = []
        for i in range(k):
            delta.append(distEclud(centroids[i], new_centroids[i]))
        for i in range(n):
            for j in range(k):
                bound_low[i, j] = max(bound_low[i, j] - delta[j], 0)
            bound_up[i] = bound_up[i] + delta[cluster_assignment[i]]
            r[i] = True
        centroids = new_centroids
        t += 1
        if sum(delta) < 10**-3 or t > 100:
            clusterChanged = False

    return centroids, cluster_assignment


df_news = pd.read_table('G:/聚类分析/数据集/gaussiandata.txt', header=None)
data = np.array(df_news)
data = data[:, :2]
centroids, clusters = selk(data, 5)
mark = ['or', 'ob', 'og', 'ok', '^r', '+r', 'sr', 'dr', '<r', 'pr']
for i in range(5):
    d = data[np.array(clusters) == i, :]
    plt.plot([d[:, 0]], [d[:, 1]], mark[i], markersize=5)
plt.show()
