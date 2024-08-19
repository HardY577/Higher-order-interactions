import numpy as np

def f_gen_trans_mat(madj, n):
    """
    生成转换矩阵，将邻接矩阵的每一行标准化为概率分布。

    参数:
        madj (numpy array): 邻接矩阵。
        n (int): 节点数。

    返回:
        numpy array: 转换后的概率矩阵。
    """
    # 将邻接矩阵转换为浮点类型
    trans_mat = madj.astype(float)
    
    # 对每一行进行标准化，使其元素之和为1（即概率分布）
    for i in range(n):
        row_sum = np.sum(trans_mat[i, :])
        if row_sum > 0:  # 避免除以零
            trans_mat[i, :] /= row_sum
        else:
            trans_mat[i, :] = 0  # 如果该行和为零，则设置为全零行
    
    return trans_mat

"""
function trans_mat = f_gen_trans_mat(madj, n)
trans_mat = double(madj);
for i = 1: n
    trans_mat(i, :) = trans_mat(i, :) / sum(trans_mat(i, :));
end
"""