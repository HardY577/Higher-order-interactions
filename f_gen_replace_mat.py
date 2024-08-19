import numpy as np
from scipy.special import comb
from f_index_dict_two import f_index_dict_two  # 假设该方法已经定义并导入

def f_gen_replace_mat(madj2, madj3, n):
    """
    生成替换矩阵，基于madj3中的高阶链接信息，更新madj2。

    参数:
        madj2 (numpy array): 原始邻接矩阵。
        madj3 (numpy array): 高阶邻接矩阵。
        n (int): 节点数。

    返回:
        numpy array: 更新后的替换矩阵。
    """
    # 初始化替换矩阵
    replace_mat = np.zeros((n, n))
    
    # 计算高阶连接数
    high_link_num = int(comb(n, 2))
    
    # 获取索引字典，假设f_index_dict_two已经导入并可用
    index_dict = f_index_dict_two(n)
    
    # 遍历所有节点和高阶链接，更新替换矩阵
    for i in range(n):
        for j in range(high_link_num):
            if madj3[i, j] > 0:
                idx = index_dict[j, 0]  # Python的索引从0开始
                idy = index_dict[j, 1]
                replace_mat[i, idx] += 1
                replace_mat[i, idy] += 1
    
    # 将madj2与生成的replace_mat相加
    replace_mat += madj2
    
    return replace_mat

"""
function replace_mat = f_gen_raplace_mat(madj2, madj3, n)
replace_mat = zeros(n, n);
high_link_num = nchoosek(n, 2);
index_dict = f_index_dict_two(n);
for i = 1: n
    for j = 1: high_link_num
        if madj3(i, j) > 0
            idx = index_dict(j, 1); idy = index_dict(j, 2);
            replace_mat(i, idx) = replace_mat(i, idx) + 1;
            replace_mat(i, idy) = replace_mat(i, idy) + 1;
        end
    end
end
replace_mat = replace_mat + madj2;
"""