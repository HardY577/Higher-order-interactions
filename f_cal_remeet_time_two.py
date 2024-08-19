import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import bicgstab
from f_find_index import f_find_index
from f_index_dict_two import f_index_dict_two

"""
重要更改说明：
1.索引调整：MATLAB中的索引是从1开始的，转换为Python代码时，所有涉及到数组索引的部分需要减1或加1来符合Python的索引规范。
2.稀疏矩阵构建：MATLAB的sparse函数在Python中使用scipy.sparse.coo_matrix来构建稀疏矩阵，并在后续转换为csc_matrix格式，以便与bicgstab函数兼容。
3.求解器：MATLAB的bicgstab在Python中使用scipy.sparse.linalg.bicgstab来替代，接口相似，但输出可能需要进一步确认。
潜在问题：
1.性能：Python中的稀疏矩阵处理速度可能较MATLAB有所不同，尤其在处理大型稀疏矩阵时需要留意性能。
2.数值稳定性：bicgstab求解器的数值稳定性可能因实现差异而有所变化，需要在实际运行时验证求解的精度和速度。
"""

def f_cal_remeet_time_two(trans_mat, n):
    index_dict = f_index_dict_two(n)
    dict_size = n * (n - 1) // 2  # nchoosek(n, 2) in Python
    storage_size = n
    
    # 初始化稀疏矩阵相关的数组
    id_x_arr = np.zeros(storage_size, dtype=int)
    id_y_arr = np.zeros(storage_size, dtype=int)
    val_arr = np.zeros(storage_size)
    
    b_arr = np.full(dict_size, -0.5)
    
    tik = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            seq = [i + 1, j + 1]  # Python中的索引需要加1
            index_idx = f_find_index(seq, index_dict)
            id_x_arr[tik] = index_idx
            id_y_arr[tik] = index_idx
            val_arr[tik] = -1
            tik += 1
            
            # 遍历转移矩阵
            for k in range(n):
                if k == j or trans_mat[i, k] == 0:
                    continue
                seq = [k + 1, j + 1]
                index_idy = f_find_index(seq, index_dict)
                id_x_arr[tik] = index_idx
                id_y_arr[tik] = index_idy
                val_arr[tik] = trans_mat[i, k] / 2
                tik += 1
                
            for k in range(n):
                if k == i or trans_mat[j, k] == 0:
                    continue
                seq = [i + 1, k + 1]
                index_idy = f_find_index(seq, index_dict)
                id_x_arr[tik] = index_idx
                id_y_arr[tik] = index_idy
                val_arr[tik] = trans_mat[j, k] / 2
                tik += 1
    
    # 构建稀疏矩阵
    adj_mat = coo_matrix((val_arr, (id_x_arr, id_y_arr)), shape=(dict_size, dict_size)).tocsc()
    
    # 使用bicgstab求解
    retime, flag = bicgstab(adj_mat, b_arr)
    
    return retime, flag

"""
function [retime, flag]= f_cal_remeet_time_two(trans_mat, n)
index_dict = f_index_dict_two(n);
dict_size = nchoosek(n, 2);
storage_size = n;
id_x_arr = zeros(storage_size, 1); id_y_arr = zeros(storage_size, 1);
val_arr = zeros(storage_size, 1);
b_arr = zeros(dict_size, 1) - 1/2;
tik = 1;
for i = 1: n
    for j = i+1: n
        seq = [i, j];
        index_idx = f_find_index(seq, index_dict);
        id_x_arr(tik) = index_idx;
        id_y_arr(tik) = index_idx;
        val_arr(tik) = -1;
        tik = tik + 1;
        for k = 1: n
            if k == j || trans_mat(i, k) == 0
                continue
            end
            seq = [k, j];
            index_idy = f_find_index(seq, index_dict);
            id_x_arr(tik) = index_idx;
            id_y_arr(tik) = index_idy;
            val_arr(tik) = trans_mat(i, k) / 2;
            tik = tik + 1;
        end
        for k = 1: n
            if k == i || trans_mat(j, k) == 0
                continue
            end
            seq = [i, k];
            index_idy = f_find_index(seq, index_dict);
            id_x_arr(tik) = index_idx;
            id_y_arr(tik) = index_idy;
            val_arr(tik) = trans_mat(j, k) / 2;
            tik = tik + 1;
        end
    end
end

adj_mat = -sparse(id_x_arr, id_y_arr, val_arr, dict_size, dict_size);
clear id_x_arr id_y_arr val_arr index_dict
b_arr = -sparse(b_arr);
[retime, flag] = bicgstab(adj_mat, b_arr);
"""

