import numpy as np
from scipy.special import comb
import f_index_dict_two


def f_gen_high_mat(madj, n):
    """
    生成high_mat矩阵。

    参数:
        madj (numpy array): 邻接矩阵。
        n (int): 节点数。

    返回:
        numpy array: high_mat矩阵。
    """
    # 计算组合数 nC2
    size = int(comb(n, 2))
    
    # 获取两两组合的索引
    index_dict = f_index_dict_two(n)
    
    # 初始化 high_mat 矩阵为零矩阵
    high_mat = np.zeros((n, size), dtype=int)
    
    # 计算 high_mat
    for i in range(n):
        for j in range(size):
            idx = index_dict[j, 0]
            idy = index_dict[j, 1]
            high_mat[i, j] = madj[i, idx] * madj[idx, idy] * madj[idy, i]
    
    return high_mat

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function high_mat = f_gen_high_mat(madj, n)
size = nchoosek(n, 2);
index_dict = f_index_dict_two(n);
high_mat = zeros(n, size);
for i = 1: n
    for j = 1: size
        idx = index_dict(j, 1); idy = index_dict(j, 2);
        high_mat(i, j) = madj(i, idx) * madj(idx, idy) * madj(idy, i);
    end
end
"""