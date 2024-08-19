import numpy as np
from scipy.special import comb

def f_index_dict_three(n):
    """
    生成包含所有三三组合的索引矩阵。

    参数:
        n (int): 节点数。

    返回:
        numpy array: 包含所有三三组合的索引矩阵。
    """
    order = 3
    num_of_index = int(comb(n, order))
    
    # 初始化 index_dict_three 矩阵
    index_dict_three = np.zeros((num_of_index, order), dtype=int)
    
    tik = 0
    for i in range(0, n):  # 注意Python索引从1开始
        for j in range(i, n):
            for k in range(j, n):
                index_dict_three[tik, 0] = i
                index_dict_three[tik, 1] = j
                index_dict_three[tik, 2] = k
                tik += 1
    
    return index_dict_three
"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function index_dict_three = f_index_dict_three(n)
order = 3;
num_of_index = nchoosek(n, order);
index_dict_three = zeros(num_of_index, order);
tik = 1;
for i = 1: n
    for j = i+1: n
        for k = j+1: n
            index_dict_three(tik, 1) = i;
            index_dict_three(tik, 2) = j;
            index_dict_three(tik, 3) = k;
            tik = tik + 1;
        end
    end
end
"""
