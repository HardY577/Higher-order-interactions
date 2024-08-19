import numpy as np

def f_gen_star_net(n):
    """
    生成一个星形网络矩阵。节点1与所有其他节点相连，偶数节点与相邻奇数节点相连。

    参数:
        n (int): 星形网络的中心节点数量。

    返回:
        numpy array: 生成的网络邻接矩阵。
    """
    # 矩阵的大小为 2*n + 1
    size = 2 * n + 1
    
    # 初始化邻接矩阵
    mat = np.zeros((size, size), dtype=int)
    
    # 遍历节点，构建连接关系
    for i in range(size):
        if i == 0:  # 第一个节点与所有其他节点相连
            for j in range(i + 1, size):
                mat[i, j] = 1
                mat[j, i] = 1
        else:
            if i % 2 == 1:  # 偶数节点与其后的奇数节点相连
                mat[i, i + 1] = 1
                mat[i + 1, i] = 1
    
    return mat

"""
奇偶节点的连接：在遍历过程中，确保偶数节点（从0开始的奇数索引节点）与其后的节点相连，
这一逻辑与MATLAB代码一致，但需要确保这种连接方式适合你的应用场景。
function mat = f_gen_star_net(n)
size = 2 * n + 1;
mat = zeros(size, size);
for i = 1: size
    if i == 1
        for j = i+1: size
            mat(i, j) = 1;
            mat(j, i) = 1;
        end
    else
        if mod(i, 2) == 0
            mat(i, i+1) = 1;
            mat(i+1, i) = 1;
        end
    end
end
"""