import numpy as np

def f_gen_conjoin_rich_club(m1, n1, m2, n2):
    """
    生成一个连接两个rich-club社团的矩阵。

    参数:
        m1 (int): 第一个社团的节点数。
        n1 (int): 第一个社团的外围节点数。
        m2 (int): 第二个社团的节点数。
        n2 (int): 第二个社团的外围节点数。

    返回:
        numpy array: 生成的矩阵，表示两个rich-club社团的连接。
    """
    size1 = m1 + n1
    size2 = m2 + n2
    size = size1 + size2
    
    # 初始化矩阵为零矩阵
    mat = np.zeros((size, size), dtype=int)
    
    # 构造第一个社团的连边
    for i in range(m1):
        for j in range(i + 1, size1):
            mat[i, j] = 1
            mat[j, i] = 1
    
    # 构造第二个社团的连边
    for i in range(size1, size1 + m2):
        for j in range(i + 1, size):
            mat[i, j] = 1
            mat[j, i] = 1
    
    # 连接两个社团的外围节点
    if n1 > 0 and n2 > 0:
        for i in range(m1):
            mat[i, size1 + i] = 1
            mat[size1 + i, i] = 1
    
    return mat

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
    function mat = f_gen_conjoin_rich_club(m1, n1, m2, n2)
    size1 = m1 + n1;
    size2 = m2 + n2;
    size = size1 + size2;
    mat = zeros(size, size);
    
    for i = 1: m1
        for j = i+1: size1
            mat(i, j) = 1;
            mat(j, i) = 1;
        end
    end

    for i = size1+1: (size1+m2)
        for j = i+1: size
            mat(i, j) = 1;
            mat(j, i) = 1;
        end
    end

    if n1 > 0 && n2 > 0
        for i = 1: m1
            mat(i, size1+i) = 1;
            mat(size1+i, i) = 1;
        end
    end
end
"""