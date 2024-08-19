import numpy as np
from scipy.special import comb

def f_index_dict_four(n):
    """
    生成包含所有四四组合的索引矩阵。

    参数:
        n (int): 节点数。

    返回:
        numpy array: 包含所有四四组合的索引矩阵。
    """
    order = 4
    num_of_index = int(comb(n, order))  # 使用scipy.special.comb计算组合数
    
    # 初始化 index_dict_four 矩阵，类型为整数
    index_dict_four = np.zeros((num_of_index, order), dtype=int)
    
    tik = 0
    # 遍历所有可能的组合，填充矩阵
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            for k in range(j + 1, n + 1):
                for l in range(k + 1, n + 1):
                    # 将组合的索引放入矩阵的相应位置
                    index_dict_four[tik, 0] = i
                    index_dict_four[tik, 1] = j
                    index_dict_four[tik, 2] = k
                    index_dict_four[tik, 3] = l
                    tik += 1  # 更新tik索引

    return index_dict_four

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function index_dict_four = f_index_dict_four(n)
order = 4;
num_of_index = nchoosek(n, order);
index_dict_four = zeros(num_of_index, order);
tik = 1;
for i = 1: n
    for j = i+1: n
        for k = j+1: n
            for l = k+1: n
                index_dict_four(tik, 1) = i;
                index_dict_four(tik, 2) = j;
                index_dict_four(tik, 3) = k;
                index_dict_four(tik, 4) = l;
                tik = tik + 1;
            end
        end
    end
end
"""
