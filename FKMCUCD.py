# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 13:29:23 2020

@author: Huang
"""
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt


def distEclud(vecA, vecB):
    # vecA,vecB是数组形式，列表形式不行
    return sum((vecA-vecB)**2)**0.5


def FKMCUCD(X, k, n):
    N, d = np.shape(X)

    # 初始类中心
    index = random.sample(range(0, N), k)
    centroids = X[index, :]

    # 计算每个点到类中心的距离，并记录距离该点最近的n个中心
    distMat = np.zeros((N, k))
    W = []
    for i in range(N):
        for j in range(k):
            distMat[i, j] = distEclud(X[i, :], centroids[j, :])
        W.append(list(distMat[i, :].argsort()[:n]))
    r_max = N*[0]
    for i in range(N):
        indx = [x[n-1] for x in W]
        r_max[i] = distMat[i, W[i][n-1]]

    # 类中心的更新，并记录那些类中心未收敛
    centroids_prime = np.zeros((k, d))
    D = k*[0]
    sc_active = []
    sc_static = []
    active = k*[1]
    for i in range(k):
        indx = [x[0] for x in W]
        temp = X[np.array(indx) == i]
        centroids_prime[i, :] = np.mean(temp, axis=0)
        D[i] = distEclud(centroids[i, :], centroids_prime[i, :])
        if D[i] < 10 ** -3:
            sc_static.append(i)
            active[i] = 0
        else:
            sc_active.append(i)
    clusterChanged = True

    # 迭代直至收敛
    while clusterChanged:
        for i in range(N):
            W_active = []
            W_static = []
            W_t = []
            for j in W[i]:
                if active[j] == 1:
                    W_active.append(j)
                else:
                    W_static.append(j)
            for j in W_active:
                dt = distEclud(X[i, :], centroids_prime[j, :])
                if dt < r_max[i]:
                    W_t.append(j)
                if not W_static:
                    for t in W_static:
                        W_t.append(t)
            if len(W_t):
                rm = 0
                temp = []
                for j in W_t:
                    dt = distEclud(X[i, :], centroids_prime[j, :])
                    temp.append(dt)
                    if dt > rm:
                        rm = dt
                for j in sc_active:
                    if j not in W_active and r_max[i]-D[j] < rm:
                        dt = distEclud(X[i, :], centroids_prime[j, :])
                        if dt <= rm:
                            W_t.append(j)
                            temp.append(dt)
                W[i] = []
                for m in list(np.array(temp).argsort()):
                    W[i].append(W_t[m])
                r_max[i] = rm
            else:
                for j in range(k):
                    distMat[i, j] = distEclud(X[i, :], centroids_prime[j, :])
                W[i] = list(distMat[i, :].argsort()[:n])
                r_max[i] = distMat[i, W[i][n-1]]
            centroids = centroids_prime
            for i in range(k):
                indx = [x[0] for x in W]
                temp = X[np.array(indx) == i]
                centroids_prime[i, :] = np.mean(temp, axis=0)
                D[i] = distEclud(centroids[i, :], centroids_prime[i, :])
                if D[i] < 10 ** -8:
                    sc_static.append(i)
                    active[i] = 0
                else:
                    sc_active.append(i)
        if sum(D) < 10 ** -8:
            clusterChanged = False
    return centroids, W


df_news = pd.read_table('G:/聚类分析/数据集/gaussiandata.txt', header=None)
data = np.array(df_news)
data = data[:, :2]
centroids, W = FKMCUCD(data, 5, 4)
clusters = [x[0] for x in W]
mark = ['or', 'ob', 'og', 'ok', '^r', '+r', 'sr', 'dr', '<r', 'pr']
for i in range(5):
    d = data[np.array(clusters) == i, :]
    plt.plot([d[:, 0]], [d[:, 1]], mark[i], markersize=5)
plt.show()
