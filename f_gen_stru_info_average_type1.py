import numpy as np
from f_index_dict_two import f_index_dict_two

"""
重要更改说明：
1.索引转换：MATLAB中的index_dict(j, 1)和index_dict(j, 2)在Python中对应index_dict[j, 0]和index_dict[j, 1]，索引从1调整为0。
2.np.sum()替代MATLAB的sum()：Python使用np.sum()来处理数组求和。
3.T_tol的计算：为保证除法操作不出错，添加了对T_tol[i]是否为0的检查。
潜在问题：
1.除以零的处理：在计算t12, t13, t22, t23, t33时，我们假设T_tol[i]不会为零。如果有可能为零，可能需要处理分母为零的情况，避免出现异常。
"""
def f_gen_stru_info_average_type1(madj2, madj3, n):
    t12 = np.zeros(n)
    t13 = np.zeros(n)
    t22 = np.zeros((n, n))
    t23 = np.zeros((n, n))
    size = n * (n - 1) // 2  # 对应 nchoosek(n, 2)
    t33 = np.zeros((n, size))
    index_dict = f_index_dict_two(n)
    T_tol = np.zeros(n)
    
    # 计算 T_tol
    for i in range(n):
        T_tol[i] = np.sum(madj2[i, :]) + np.sum(madj3[i, :])
    
    # 计算 t12 和 t13
    for i in range(n):
        if T_tol[i] != 0:
            t12[i] = np.sum(madj2[i, :]) / T_tol[i]
            t13[i] = np.sum(madj3[i, :]) / T_tol[i]
    
    # 计算 t22
    for i in range(n):
        for j in range(n):
            if T_tol[i] != 0:
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
function [t12, t13, t22, t23, t33] = f_gen_stru_info_average_type1(madj2, madj3, n)
t12 = zeros(n, 1); t13 = zeros(n, 1);
t22 = zeros(n, n); t23 = zeros(n, n);
size = nchoosek(n, 2);
t33 = zeros(n, size);
index_dict = f_index_dict_two(n);
T_tol = zeros(n);
for i = 1: n
    T_tol(i) = sum(madj2(i, :)) + sum(madj3(i, :));
end

for i = 1: n
    t12(i) = sum(madj2(i, :)) / T_tol(i);
    t13(i) = sum(madj3(i, :)) / T_tol(i);
end

for i = 1: n
    for j = 1: n
        t22(i, j) = madj2(i, j) / T_tol(i);
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