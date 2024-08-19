import numpy as np
from scipy.special import comb

from f_index_dict_two import f_index_dict_two
"""
重要的更改说明：
1.组合数计算：使用 scipy.special.comb 来计算组合数 nchoosek(n, 2)，即从 n 个元素中选出 2 个的组合数目，并将结果存储在 size 变量中。
2.注释细化：每个逻辑块都有详细的注释，解释了其计算过程及判断条件。
3.函数调用：假设 f_index_dict_two 函数已经存在，并按照 MATLAB 版本的功能进行了实现。

潜在问题：
1.如果 madj2 和 madj3 具有稀疏矩阵的性质，可以考虑使用 scipy.sparse 来优化存储和计算，以提高效率。
"""
def f_gen_stru_info_average_type2(madj2, madj3, n):
    # 初始化t12, t13为大小为(n, 1)的零数组
    t12 = np.zeros(n)
    t13 = np.zeros(n)
    
    # 初始化t22为大小为(n, n)的零数组
    t22 = np.zeros((n, n))
    t23 = np.zeros((n, n))
    
    # 计算组合数nchoosek(n, 2)，即从n个元素中选2个的组合数目
    size = comb(n, 2, exact=True)
    
    # 初始化t33为大小为(n, size)的零数组
    t33 = np.zeros((n, size))
    
    # 调用f_index_dict_two函数生成索引字典
    index_dict = f_index_dict_two(n)
    
    # 计算t12和t13
    for i in range(n):
        if np.sum(madj2[i, :]) > 0:
            # 如果madj2第i行的和大于0，则计算t12[i]
            t12[i] = np.sum(madj2[i, :]) / np.sum(madj2[i, :])
        
        if np.sum(madj3[i, :]) > 0:
            # 如果madj3第i行的和大于0，则计算t13[i]
            t13[i] = np.sum(madj3[i, :]) / np.sum(madj3[i, :])
    
    # 计算t22
    for i in range(n):
        for j in range(n):
            if np.sum(madj2[i, :]) > 0:
                # 如果madj2第i行的和大于0，则计算t22[i, j]
                t22[i, j] = madj2[i, j] / np.sum(madj2[i, :])
    
    # 计算t23
    for i in range(n):
        for j in range(size):
            if madj3[i, j] > 0 and np.sum(madj3[i, :]) > 0:
                # 如果madj3第i行第j列大于0且madj3第i行的和大于0
                idx = index_dict[j, 0]
                idy = index_dict[j, 1]
                t23[i, idx] += 1 / np.sum(madj3[i, :])
                t23[i, idy] += 1 / np.sum(madj3[i, :])
    
    # 计算t33
    for i in range(n):
        for j in range(size):
            if madj3[i, j] > 0 and np.sum(madj3[i, :]) > 0:
                # 如果madj3第i行第j列大于0且madj3第i行的和大于0
                t33[i, j] = madj3[i, j] / np.sum(madj3[i, :])
    
    # 返回计算的t12, t13, t22, t23, t33
    return t12, t13, t22, t23, t33


"""
function [t12, t13, t22, t23, t33] = f_gen_stru_info_average_type2(madj2, madj3, n)
t12 = zeros(n, 1); t13 = zeros(n, 1);
t22 = zeros(n, n); t23 = zeros(n, n);
size = nchoosek(n, 2);
t33 = zeros(n, size);
index_dict = f_index_dict_two(n);

for i = 1: n
    if sum(madj2(i, :)) > 0
        t12(i) = sum(madj2(i, :)) / sum(madj2(i, :));
    end
    if sum(madj3(i, :)) > 0
        t13(i) = sum(madj3(i, :)) / sum(madj3(i, :));
    end
end

for i = 1: n
    for j = 1: n
        if sum(madj2(i, :)) > 0
            t22(i, j) = madj2(i, j) / sum(madj2(i, :));
        end
    end
end

for i = 1: n
    for j = 1: size
        if madj3(i, j) > 0 && sum(madj3(i, :)) > 0
            idx = index_dict(j, 1); idy = index_dict(j, 2);
            t23(i, idx) = t23(i, idx) + 1 / sum(madj3(i, :));
            t23(i, idy) = t23(i, idy) + 1 / sum(madj3(i, :));
        end
    end
end

for i = 1: n
    for j = 1: size
        if madj3(i, j) > 0 && sum(madj3(i, :)) > 0
            t33(i, j) = madj3(i, j) / sum(madj3(i, :));
        end
    end
end
"""