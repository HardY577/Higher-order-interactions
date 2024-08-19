import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import bicgstab

from f_find_index import f_find_index
from f_index_dict_four import f_index_dict_four
from f_index_dict_three import f_index_dict_three
from f_index_dict_two import f_index_dict_two

"""
1.动态数组扩展：由于id_x_arr、id_y_arr和val_arr数组的大小不确定，我在每次添加数据时检查数组是否需要扩展，并使用np.append动态增加数组大小。这在Python中比预先定义一个固定大小的数组更为灵活。
2.组合计算：np.math.comb(n, 4)用于计算组合数量，替代了MATLAB中的nchoosek函数。
3.稀疏矩阵：使用scipy.sparse.coo_matrix来构建稀疏矩阵，并使用scipy.sparse.linalg.bicgstab求解线性方程组。
4.代码逻辑：四重循环的逻辑保持与原始MATLAB代码一致，主要是计算不同的组合并更新b_arr和稀疏矩阵adj_mat。
"""
def f_cal_remeet_time_four(trans_mat, retime2, retime3, n):
    # 获取不同组合的索引字典
    index_dict_two = f_index_dict_two(n)
    index_dict_three = f_index_dict_three(n)
    index_dict_four = f_index_dict_four(n)
    
    # 计算4元素组合的数量
    dict_size = int(np.math.comb(n, 4))  # 使用Python中的组合公式
    
    # 初始化存储数组
    storage_size = 1
    id_x_arr = np.zeros(storage_size, dtype=int)
    id_y_arr = np.zeros(storage_size, dtype=int)
    val_arr = np.zeros(storage_size)
    
    # 初始化b数组，并将其设置为-1/4
    b_arr = np.full(dict_size, -1/4)
    
    # tik1用于追踪当前索引，tik2用于追踪存储数组的索引
    tik1 = 0
    tik2 = 0
    
    # 四重循环用于遍历i, j, k, g组合
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for g in range(k + 1, n):
                    index_idx = tik1  # 当前索引
                    
                    # 更新id_x_arr和id_y_arr数组，表示对角元素，值为-1
                    if tik2 >= storage_size:
                        # 动态增加存储空间，如果tik2超出当前数组大小
                        id_x_arr = np.append(id_x_arr, np.zeros(storage_size, dtype=int))
                        id_y_arr = np.append(id_y_arr, np.zeros(storage_size, dtype=int))
                        val_arr = np.append(val_arr, np.zeros(storage_size))
                    
                    id_x_arr[tik2] = index_idx
                    id_y_arr[tik2] = index_idx
                    val_arr[tik2] = -1
                    tik2 += 1
                    
                    # 对每个l值进行遍历，更新b_arr和存储数组
                    for l in range(n):
                        if trans_mat[i, l] == 0:
                            continue  # 如果转移矩阵的值为0，跳过
                        
                        # 生成唯一组合
                        seq = np.unique([l, j, k, g])
                        len_seq = len(seq)
                        
                        if len_seq == 1:
                            continue  # 如果组合长度为1，跳过
                        elif len_seq == 2:
                            index = f_find_index(seq, index_dict_two)
                            val = retime2[index] * trans_mat[i, l] / 4
                            b_arr[tik1] -= val  # 更新b数组
                        elif len_seq == 3:
                            index = f_find_index(seq, index_dict_three)
                            val = retime3[index] * trans_mat[i, l] / 4
                            b_arr[tik1] -= val
                        else:
                            index_idy = f_find_index(seq, index_dict_four)
                            if tik2 >= len(id_x_arr):
                                # 动态增加存储空间
                                id_x_arr = np.append(id_x_arr, np.zeros(storage_size, dtype=int))
                                id_y_arr = np.append(id_y_arr, np.zeros(storage_size, dtype=int))
                                val_arr = np.append(val_arr, np.zeros(storage_size))
                            
                            id_x_arr[tik2] = index_idx
                            id_y_arr[tik2] = index_idy
                            val_arr[tik2] = trans_mat[i, l] / 4
                            tik2 += 1
                    
                    # 同样的操作应用于j, k, g的位置，分别遍历l值
                    # 循环逻辑与i位置相同，这里省略重复注释
                    for l in range(n):
                        if trans_mat[j, l] == 0:
                            continue
                        seq = np.unique([i, l, k, g])
                        len_seq = len(seq)
                        
                        if len_seq == 1:
                            continue
                        elif len_seq == 2:
                            index = f_find_index(seq, index_dict_two)
                            val = retime2[index] * trans_mat[j, l] / 4
                            b_arr[tik1] -= val
                        elif len_seq == 3:
                            index = f_find_index(seq, index_dict_three)
                            val = retime3[index] * trans_mat[j, l] / 4
                            b_arr[tik1] -= val
                        else:
                            index_idy = f_find_index(seq, index_dict_four)
                            if tik2 >= len(id_x_arr):
                                id_x_arr = np.append(id_x_arr, np.zeros(storage_size, dtype=int))
                                id_y_arr = np.append(id_y_arr, np.zeros(storage_size, dtype=int))
                                val_arr = np.append(val_arr, np.zeros(storage_size))
                            id_x_arr[tik2] = index_idx
                            id_y_arr[tik2] = index_idy
                            val_arr[tik2] = trans_mat[j, l] / 4
                            tik2 += 1
                    
                    for l in range(n):
                        if trans_mat[k, l] == 0:
                            continue
                        seq = np.unique([i, j, l, g])
                        len_seq = len(seq)
                        
                        if len_seq == 1:
                            continue
                        elif len_seq == 2:
                            index = f_find_index(seq, index_dict_two)
                            val = retime2[index] * trans_mat[k, l] / 4
                            b_arr[tik1] -= val
                        elif len_seq == 3:
                            index = f_find_index(seq, index_dict_three)
                            val = retime3[index] * trans_mat[k, l] / 4
                            b_arr[tik1] -= val
                        else:
                            index_idy = f_find_index(seq, index_dict_four)
                            if tik2 >= len(id_x_arr):
                                id_x_arr = np.append(id_x_arr, np.zeros(storage_size, dtype=int))
                                id_y_arr = np.append(id_y_arr, np.zeros(storage_size, dtype=int))
                                val_arr = np.append(val_arr, np.zeros(storage_size))
                            id_x_arr[tik2] = index_idx
                            id_y_arr[tik2] = index_idy
                            val_arr[tik2] = trans_mat[k, l] / 4
                            tik2 += 1
                    
                    for l in range(n):
                        if trans_mat[g, l] == 0:
                            continue
                        seq = np.unique([i, j, k, l])
                        len_seq = len(seq)
                        
                        if len_seq == 1:
                            continue
                        elif len_seq == 2:
                            index = f_find_index(seq, index_dict_two)
                            val = retime2[index] * trans_mat[g, l] / 4
                            b_arr[tik1] -= val
                        elif len_seq == 3:
                            index = f_find_index(seq, index_dict_three)
                            val = retime3[index] * trans_mat[g, l] / 4
                            b_arr[tik1] -= val
                        else:
                            index_idy = f_find_index(seq, index_dict_four)
                            if tik2 >= len(id_x_arr):
                                id_x_arr = np.append(id_x_arr, np.zeros(storage_size, dtype=int))
                                id_y_arr = np.append(id_y_arr, np.zeros(storage_size, dtype=int))
                                val_arr = np.append(val_arr, np.zeros(storage_size))
                            id_x_arr[tik2] = index_idx
                            id_y_arr[tik2] = index_idy
                            val_arr[tik2] = trans_mat[g, l] / 4
                            tik2 += 1
                    
                    tik1 += 1  # 更新tik1计数器
    
    # 构建稀疏矩阵adj_mat
    adj_mat = coo_matrix((val_arr[:tik2], (id_x_arr[:tik2], id_y_arr[:tik2])), shape=(dict_size, dict_size))
    
    # 将b_arr转换为稀疏矩阵格式
    b_arr_sparse = coo_matrix(-b_arr).T
    
    # 使用双共轭梯度稳定法(bicgstab)求解线性方程组
    retime, _ = bicgstab(adj_mat, b_arr_sparse.toarray().flatten(), tol=1e-5, maxiter=10**5)
    
    return retime


"""
function [retime] = f_cal_remeet_time_four(trans_mat, retime2, retime3, n)
index_dict_two = f_index_dict_two(n);
index_dict_three = f_index_dict_three(n);
index_dict_four = f_index_dict_four(n);
dict_size = nchoosek(n, 4);
storage_size = 1;
id_x_arr = zeros(storage_size, 1); id_y_arr = zeros(storage_size, 1);
val_arr = zeros(storage_size, 1);
b_arr = zeros(dict_size, 1) - 1/4;
tik1 = 1;
tik2 = 1;
for i = 1: n
    for j = i+1: n
        for k = j+1: n
            for g = k+1: n
                index_idx = tik1;
                id_x_arr(tik2) = index_idx;
                id_y_arr(tik2) = index_idx;
                val_arr(tik2) = -1;
                tik2 = tik2 + 1;
                for l = 1: n
                    if trans_mat(i, l) == 0
                        continue
                    end
                    seq = [l, j, k, g];
                    seq = unique(seq);
                    len = length(seq);
                    if len == 1
                        continue
                    elseif len == 2
                        index = f_find_index(seq, index_dict_two);
                        val = retime2(index) * trans_mat(i, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    elseif len == 3
                        index = f_find_index(seq, index_dict_three);
                        val = retime3(index) * trans_mat(i, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    else
                        index_idy = f_find_index(seq, index_dict_four);
                        id_x_arr(tik2) = index_idx;
                        id_y_arr(tik2) = index_idy;
                        val_arr(tik2) = trans_mat(i, l) / 4;
                        tik2 = tik2 + 1;
                    end
                end
                for l = 1: n
                    if trans_mat(j, l) == 0
                        continue
                    end
                    seq = [i, l, k, g];
                    seq = unique(seq);
                    len = length(seq);
                    if len == 1
                        continue
                    elseif len == 2
                        index = f_find_index(seq, index_dict_two);
                        val = retime2(index) * trans_mat(j, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    elseif len == 3
                        index = f_find_index(seq, index_dict_three);
                        val = retime3(index) * trans_mat(j, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    else
                        index_idy = f_find_index(seq, index_dict_four);
                        id_x_arr(tik2) = index_idx;
                        id_y_arr(tik2) = index_idy;
                        val_arr(tik2) = trans_mat(j, l) / 4;
                        tik2 = tik2 + 1;
                    end
                end
                for l = 1: n
                    if trans_mat(k, l) == 0
                        continue
                    end
                    seq = [i, j, l, g];
                    seq = unique(seq);
                    len = length(seq);
                    if len == 1
                        continue
                    elseif len == 2
                        index = f_find_index(seq, index_dict_two);
                        val = retime2(index) * trans_mat(k, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    elseif len == 3
                        index = f_find_index(seq, index_dict_three);
                        val = retime3(index) * trans_mat(k, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    else
                        index_idy = f_find_index(seq, index_dict_four);
                        id_x_arr(tik2) = index_idx;
                        id_y_arr(tik2) = index_idy;
                        val_arr(tik2) = trans_mat(k, l) / 4;
                        tik2 = tik2 + 1;
                    end
                end
                for l = 1: n
                    if trans_mat(g, l) == 0
                        continue
                    end
                    seq = [i, j, k, l];
                    seq = unique(seq);
                    len = length(seq);
                    if len == 1
                        continue
                    elseif len == 2
                        index = f_find_index(seq, index_dict_two);
                        val = retime2(index) * trans_mat(g, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    elseif len == 3
                        index = f_find_index(seq, index_dict_three);
                        val = retime3(index) * trans_mat(g, l) / 4;
                        b_arr(tik1) = b_arr(tik1) - val;
                    else
                        index_idy = f_find_index(seq, index_dict_four);
                        id_x_arr(tik2) = index_idx;
                        id_y_arr(tik2) = index_idy;
                        val_arr(tik2) = trans_mat(g, l) / 4;
                        tik2 = tik2 + 1;
                    end
                end
            tik1 = tik1 + 1;
            end
        end
    end
end

adj_mat = -sparse(id_x_arr, id_y_arr, val_arr, dict_size, dict_size);
clear id_x_arr id_y_arr val_arr index_dict
b_arr = -sparse(b_arr);
retime = bicgstab(adj_mat, b_arr, [], 10^5);
"""
