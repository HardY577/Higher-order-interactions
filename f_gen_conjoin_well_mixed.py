import numpy as np

def f_gen_conjoin_well_mixed(n1, n2):
    """
    生成一个混合连接矩阵，表示两个社团之间的连接。

    参数:
        n1 (int): 第一个社团的节点数。
        n2 (int): 第二个社团的节点数。

    返回:
        numpy array: 生成的混合连接矩阵。
    """
    size = n1 + n2
    
    # 初始化矩阵为零矩阵
    mat = np.zeros((size, size), dtype=int)
    
    # 构造第一个社团的完全连边
    for i in range(n1):
        for j in range(i + 1, n1):
            mat[i, j] = 1
            mat[j, i] = 1
    
    # 构造第二个社团的完全连边
    for i in range(n1, size):
        for j in range(i + 1, size):
            mat[i, j] = 1
            mat[j, i] = 1
    
    # 连接两个社团的一个节点
    if n1 > 0 and n2 > 0:
        mat[0, n1] = 1
        mat[n1, 0] = 1
    
    return mat

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function mat = f_gen_conjoin_well_mixed(n1, n2)
size = n1 + n2;
mat = zeros(size, size);
for i = 1: n1
    for j = i+1: n1
        mat(i, j) = 1;
        mat(j, i) = 1;
    end
end

for i = n1+1: size
    for j = i+1: size
        mat(i, j) = 1;
        mat(j, i) = 1;
    end
end

if n1 > 0 && n2 > 0
    mat(1, n1+1) = 1;
    mat(n1+1, 1) = 1;
end
"""