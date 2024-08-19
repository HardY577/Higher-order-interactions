import numpy as np

def f_node_neibour(madj, n):
    """
    计算每个节点的邻居记录及邻居数量。

    参数:
        madj (numpy array): 邻接矩阵。
        n (int): 节点数。

    返回:
        tuple: 包含邻居记录和邻居数量的元组 (neibour_record, number_record)。
    """
    # 初始化记录数组
    number_record = np.zeros(n, dtype=int)
    neibour_record = np.full((n, n), -1, dtype=int)
    
    # 遍历邻接矩阵，记录邻居节点和数量
    for i in range(n):
        for j in range(i, n):
            if madj[i, j] > 0:
                number_record[i] += 1
                number_record[j] += 1
                neibour_record[i, number_record[i] - 1] = j  # 索引从0开始
                neibour_record[j, number_record[j] - 1] = i
    
    return neibour_record, number_record

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function [neibour_record, number_record] = f_node_neibour(madj, n)
number_record = zeros(1, n);
neibour_record = zeros(n, n) - 1;
for i = 1 : n
    for j = i : n
        if madj(i, j) > 0
            number_record(i) = number_record(i) + 1; 
            number_record(j) = number_record(j) + 1;
            neibour_record(i, number_record(i)) = j;
            neibour_record(j, number_record(j)) = i;
        end
    end
end
"""