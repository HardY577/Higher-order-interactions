import numpy as np
from scipy.special import comb

def f_index_dict_two(n):
    """
    生成包含所有两两组合的索引矩阵。

    参数:
        n (int): 节点数。

    返回:
        numpy array: 包含所有两两组合的索引矩阵。
    """
    order = 2
    num_of_index = int(comb(n, order))
    
    # 初始化 index_dict_two 矩阵
    index_dict_two = np.zeros((num_of_index, order), dtype=int)
    
    tik = 0
    for i in range(1, n + 1):  # 注意Python索引从1开始
        for j in range(i + 1, n + 1):
            index_dict_two[tik, 0] = i
            index_dict_two[tik, 1] = j
            tik += 1
    
    return index_dict_two

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function index_dict_two = f_index_dict_two(n)
order = 2;
num_of_index = nchoosek(n, order);
index_dict_two = zeros(num_of_index, order);
tik = 1;
for i = 1: n
    for j = i+1: n
        index_dict_two(tik, 1) = i;
        index_dict_two(tik, 2) = j;
        tik = tik + 1;
    end
end

"""
