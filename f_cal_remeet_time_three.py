import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import bicgstab

from f_find_index import f_find_index
from f_index_dict_three import f_index_dict_three
from f_index_dict_two import f_index_dict_two

"""
重要更改说明：
1.sorted(set(seq)) 替代 unique(seq)：在Python中使用set去重，并用sorted来保证顺序一致性。
2.索引调整：同样地，所有数组索引从1开始改为从0开始。
3.b_arr更新逻辑：直接用b_arr[tik1] -= val来减少赋值次数。
4.稀疏矩阵构建：使用scipy.sparse.coo_matrix进行稀疏矩阵构建，并转换为csc格式以便与bicgstab函数兼容。
5.求解参数调整：增加tol=1e-5和maxiter=100000，与MATLAB中bicgstab默认的10^5次迭代匹配。

潜在问题：
1.性能与稀疏矩阵规模：根据矩阵规模大小，稀疏矩阵的存储方式和求解过程可能在性能上有所不同，实际运行时需要根据数据规模调优。
2.数值稳定性：仍需要在实际数据中验证求解器的数值稳定性。
"""
def f_cal_remeet_time_three(trans_mat, retime2, n):
    index_dict_two = f_index_dict_two(n)
    index_dict_three = f_index_dict_three(n)
    dict_size = n * (n - 1) * (n - 2) // 6  # 对应 nchoosek(n, 3)
    
    # 初始化存储空间
    storage_size = 1  # 根据代码，初始大小为1
    id_x_arr = np.zeros(storage_size, dtype=int)
    id_y_arr = np.zeros(storage_size, dtype=int)
    val_arr = np.zeros(storage_size)
    b_arr = np.full(dict_size, -1/3)
    
    tik1 = 0
    tik2 = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                index_idx = tik1
                id_x_arr[tik2] = index_idx
                id_y_arr[tik2] = index_idx
                val_arr[tik2] = -1
                tik2 += 1
                
                # 遍历i行的邻接矩阵
                for l in range(n):
                    if trans_mat[i, l] == 0:
                        continue
                    seq = sorted(set([l, j, k]))  # Python集合自动去重并排序
                    len_seq = len(seq)
                    
                    if len_seq == 1:
                        continue
                    elif len_seq == 2:
                        index = f_find_index(seq, index_dict_two)
                        val = retime2[index] * trans_mat[i, l] / 3
                        b_arr[tik1] -= val
                    else:
                        index_idy = f_find_index(seq, index_dict_three)
                        id_x_arr[tik2] = index_idx
                        id_y_arr[tik2] = index_idy
                        val_arr[tik2] = trans_mat[i, l] / 3
                        tik2 += 1
                
                # 遍历j行的邻接矩阵
                for l in range(n):
                    if trans_mat[j, l] == 0:
                        continue
                    seq = sorted(set([i, l, k]))
                    len_seq = len(seq)
                    
                    if len_seq == 1:
                        continue
                    elif len_seq == 2:
                        index = f_find_index(seq, index_dict_two)
                        val = retime2[index] * trans_mat[j, l] / 3
                        b_arr[tik1] -= val
                    else:
                        index_idy = f_find_index(seq, index_dict_three)
                        id_x_arr[tik2] = index_idx
                        id_y_arr[tik2] = index_idy
                        val_arr[tik2] = trans_mat[j, l] / 3
                        tik2 += 1
                
                # 遍历k行的邻接矩阵
                for l in range(n):
                    if trans_mat[k, l] == 0:
                        continue
                    seq = sorted(set([i, j, l]))
                    len_seq = len(seq)
                    
                    if len_seq == 1:
                        continue
                    elif len_seq == 2:
                        index = f_find_index(seq, index_dict_two)
                        val = retime2[index] * trans_mat[k, l] / 3
                        b_arr[tik1] -= val
                    else:
                        index_idy = f_find_index(seq, index_dict_three)
                        id_x_arr[tik2] = index_idx
                        id_y_arr[tik2] = index_idy
                        val_arr[tik2] = trans_mat[k, l] / 3
                        tik2 += 1
                
                tik1 += 1
    
    # 构建稀疏矩阵
    adj_mat = coo_matrix((val_arr, (id_x_arr, id_y_arr)), shape=(dict_size, dict_size)).tocsc()
    
    # 使用bicgstab求解
    retime, flag = bicgstab(adj_mat, b_arr, tol=1e-5, maxiter=100000)
    
    return retime, flag


"""
function [retime, flag] = f_cal_remeet_time_three(trans_mat, retime2, n)
index_dict_two = f_index_dict_two(n);
index_dict_three = f_index_dict_three(n);
dict_size = nchoosek(n, 3);
storage_size = 1;
id_x_arr = zeros(storage_size, 1); id_y_arr = zeros(storage_size, 1);
val_arr = zeros(storage_size, 1);
b_arr = zeros(dict_size, 1) - 1/3;
tik1 = 1;
tik2 = 1;
for i = 1: n
    for j = i+1: n
        for k = j+1: n
            index_idx = tik1;
            id_x_arr(tik2) = index_idx;
            id_y_arr(tik2) = index_idx;
            val_arr(tik2) = -1;
            tik2 = tik2 + 1;
            for l = 1: n
                if trans_mat(i, l) == 0
                    continue
                end
                seq = [l, j, k];
                seq = unique(seq);
                len = length(seq);
                if len == 1
                    continue
                elseif len == 2
                    index = f_find_index(seq, index_dict_two);
                    val = retime2(index) * trans_mat(i, l) / 3;
                    b_arr(tik1) = b_arr(tik1) - val;
                else
                    index_idy = f_find_index(seq, index_dict_three);
                    id_x_arr(tik2) = index_idx;
                    id_y_arr(tik2) = index_idy;
                    val_arr(tik2) = trans_mat(i, l) / 3;
                    tik2 = tik2 + 1;
                end
            end
            for l = 1: n
                if trans_mat(j, l) == 0
                    continue
                end
                seq = [i, l, k];
                seq = unique(seq);
                len = length(seq);
                if len == 1
                    continue
                elseif len == 2
                    index = f_find_index(seq, index_dict_two);
                    val = retime2(index) * trans_mat(j, l) / 3;
                    b_arr(tik1) = b_arr(tik1) - val;
                else
                    index_idy = f_find_index(seq, index_dict_three);
                    id_x_arr(tik2) = index_idx;
                    id_y_arr(tik2) = index_idy;
                    val_arr(tik2) = trans_mat(j, l) / 3;
                    tik2 = tik2 + 1;
                end
            end
            for l = 1: n
                if trans_mat(k, l) == 0
                    continue
                end
                seq = [i, j, l];
                seq = unique(seq);
                len = length(seq);
                if len == 1
                    continue
                elseif len == 2
                    index = f_find_index(seq, index_dict_two);
                    val = retime2(index) * trans_mat(k, l) / 3;
                    b_arr(tik1) = b_arr(tik1) - val;
                else
                    index_idy = f_find_index(seq, index_dict_three);
                    id_x_arr(tik2) = index_idx;
                    id_y_arr(tik2) = index_idy;
                    val_arr(tik2) = trans_mat(k, l) / 3;
                    tik2 = tik2 + 1;
                end
            end
            tik1 = tik1 + 1;
        end
    end
end

adj_mat = -sparse(id_x_arr, id_y_arr, val_arr, dict_size, dict_size);
clear id_x_arr id_y_arr val_arr index_dict
b_arr = -sparse(b_arr);
[retime, flag] = bicgstab(adj_mat, b_arr, [], 10^5);
"""