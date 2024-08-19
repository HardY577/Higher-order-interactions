import numpy as np

from f_index_dict_two import f_index_dict_two

"""
重要更改说明：
1.索引转换：MATLAB的索引从1开始，Python从0开始，因此在index_dict的访问时做了相应调整。
2.T_tol初始化：MATLAB中的zeros(n)+1在Python中被替换为np.ones(n)。
3.np.sum()替代MATLAB的sum()：Python中使用np.sum()来进行数组元素的求和操作。

潜在问题：
1.数值稳定性：与之前的代码类似，数值计算的稳定性可能需要在实际数据上进行验证。
2.性能优化：根据数据规模，可能需要对大规模矩阵计算进行进一步优化。
"""
def f_gen_stru_info_accumulate(madj2, madj3, n):
    t12 = np.zeros(n)
    t13 = np.zeros(n)
    t22 = np.zeros((n, n))
    t23 = np.zeros((n, n))
    size = n * (n - 1) // 2  # 对应 nchoosek(n, 2)
    t33 = np.zeros((n, size))
    index_dict = f_index_dict_two(n)
    T_tol = np.ones(n)  # 等效于 MATLAB 中的 zeros(n) + 1
    
    # 计算 t12 和 t13
    for i in range(n):
        if np.sum(madj2[i, :]) > 0:
            t12[i] = np.sum(madj2[i, :]) / T_tol[i]
        if np.sum(madj3[i, :]) > 0:
            t13[i] = np.sum(madj3[i, :]) / T_tol[i]
    
    # 计算 t22
    for i in range(n):
        for j in range(n):
            if np.sum(madj2[i, :]) > 0:
                t22[i, j] = madj2[i, j] / T_tol[i]
    
    # 计算 t23
    for i in range(n):
        for j in range(size):
            if madj3[i, j] > 0:
                idx = index_dict[j, 0]  # MATLAB 索引从 1 开始，Python 从 0 开始
                idy = index_dict[j, 1]
                t23[i, idx] += 1 / T_tol[i]
                t23[i, idy] += 1 / T_tol[i]
    
    # 计算 t33
    for i in range(n):
        for j in range(size):
            if madj3[i, j] > 0:
                t33[i, j] = madj3[i, j] / T_tol[i]
    
    return t12, t13, t22, t23, t33

"""
function [t12, t13, t22, t23, t33] = f_gen_stru_info_accumulate(madj2, madj3, n)
t12 = zeros(n, 1); t13 = zeros(n, 1);
t22 = zeros(n, n); t23 = zeros(n, n);
size = nchoosek(n, 2);
t33 = zeros(n, size);
index_dict = f_index_dict_two(n);
T_tol = zeros(n)+1;


for i = 1: n
    if sum(madj2(i, :)) > 0
        t12(i) = sum(madj2(i, :)) / T_tol(i);
    end
    if sum(madj3(i, :)) > 0
        t13(i) = sum(madj3(i, :)) / T_tol(i);
    end
end

for i = 1: n
    for j = 1: n
        if sum(madj2(i, :)) > 0
            t22(i, j) = madj2(i, j) / T_tol(i);
        end
    end
end

for i = 1: n
    for j = 1: size
        if madj3(i, j) > 0
            idx = index_dict(j, 1); idy = index_dict(j, 2);
            t23(i, idx) = t23(i, idx) + 1 / T_tol(i);
            t23(i, idy) = t23(i, idy) + 1 / T_tol(i);
        end
    end
end

for i = 1: n
    for j = 1: size
        if madj3(i, j) > 0
            t33(i, j) = madj3(i, j) / T_tol(i);
        end
    end
end
"""