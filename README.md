常用算法R代码实现
===
聚类算法：
------
`kernel k means`算法         kernel_k_means.R <br>
`投影聚类算法`               PROCLUS.R<br>
`密度寻峰`聚类               clustering_by_fast_search.R <br>
`基于路径的密度寻峰`聚类(Floyd算法)     path_based_floyd.R <br>
`基于路径的密度寻峰`聚类(Dijkstra算法)      path_based_dijkstra.R <br>
`meanshift`聚类算法          meanshift.R <br>
负荷数据meanshift聚类效果     08-06-01data.R <br>

曲线距离：
------
曲线的`Merge_Split_Move距离`     MSM_Distance.R <br>
曲线的`shapedtw距离`和`dtw距离`   shapeDTW.R <br>
`b样条`拟合曲线                   bspline.R <br>
距离矩阵的计算             dist_mat.R

优化：
------
二次优化        01youhua.R

生成数据集：
------
'月牙型'数据集    moon_like_data.R

其他：
------
判断一个向量是否在矩阵的行(列)中    vec_in_matrix.R	
