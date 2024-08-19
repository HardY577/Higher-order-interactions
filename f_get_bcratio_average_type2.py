import numpy as np

from f_find_index import f_find_index
from f_gen_stru_info_average_type2 import f_gen_stru_info_average_type2
from f_index_dict_four import f_index_dict_four
from f_index_dict_three import f_index_dict_three
from f_index_dict_two import f_index_dict_two

"""
重要的更改说明：
1.矩阵幂运算：使用 np.linalg.matrix_power 计算 trans_mat 的幂次。
2.单位矩阵：使用 np.eye(n) 生成单位矩阵 trans_mat0。
3.条件语句：在计算 bcratio 时，增加了对 tik2 为 0 的检查，防止出现除零错误。
4.注释细化：在每个步骤都添加了详细的中文注释，以便于理解代码逻辑。
潜在问题：
1.如果 trans_mat、madj2 和 madj3 是稀疏矩阵，考虑使用 scipy.sparse 进行优化。
2.需要确保 f_index_dict_two、f_index_dict_three、f_index_dict_four、f_find_index 
    和 f_gen_stru_info_average_type2 函数在 Python 代码中已经正确实现。
"""
def f_get_bcratio_average_type2(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2):
    # 计算转移矩阵的平方和单位矩阵
    trans_mat2 = np.linalg.matrix_power(trans_mat, 2)
    trans_mat0 = np.eye(n)  # 单位矩阵

    # 获取索引字典
    index_dict2 = f_index_dict_two(n)
    index_dict3 = f_index_dict_three(n)
    index_dict4 = f_index_dict_four(n)

    # 调用 f_gen_stru_info_average_type2 获取结构信息
    t12, t13, t22, t23, t33 = f_gen_stru_info_average_type2(madj2, madj3, n)

    # 初始化 tik1 和 tik2
    tik1 = 0
    tik2 = 0

    # 计算 tik1
    for j in range(n):
        for m in range(n):
            seq = np.unique([j, m])  # 生成唯一序列
            len_seq = len(seq)
            if len_seq == 1:
                val = 0
            else:
                index = f_find_index(seq, index_dict2)
                val = retime2[index]
            tik1 += repro_val[j] * (t12[m] + t13[m]) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 计算初始的 tik2
    for j in range(n):
        for m in range(n):
            seq = np.unique([j, m])  # 生成唯一序列
            len_seq = len(seq)
            if len_seq == 1:
                val = 0
            else:
                index = f_find_index(seq, index_dict2)
                val = retime2[index]
            tik2 += repro_val[j] * (t12[m] / 2 + t13[m] / 3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 计算含有 t22 和 t23 的 tik2
    for j in range(n):
        for m in range(n):
            for j1 in range(n):
                seq = np.unique([j, j1])  # 生成唯一序列
                len_seq = len(seq)
                if len_seq == 1:
                    val = 0
                else:
                    index = f_find_index(seq, index_dict2)
                    val = retime2[index]
                tik2 += repro_val[j] * (t22[m, j1] / 2 + t23[m, j1] / 3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 计算含有 disc1 和 disc2 的 tik2
    for j in range(n):
        for m in range(n):
            for j1 in range(n):
                seq = np.unique([j, m, j1])  # 生成唯一序列
                len_seq = len(seq)
                if len_seq == 1:
                    val = 0
                elif len_seq == 2:
                    index = f_find_index(seq, index_dict2)
                    val = retime2[index]
                else:
                    index = f_find_index(seq, index_dict3)
                    val = retime3[index]
                tik2 += repro_val[j] * ((disc1 - 1) * t22[m, j1] / 2 + (disc2 - 1) * t23[m, j1] / 3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 检查 madj3 的元素和是否大于 0
    if np.sum(madj3) > 0:
        # 计算含有 t33 的 tik2
        for j in range(n):
            for m in range(n):
                for j1 in range(n):
                    for j2 in range(n):
                        seq = np.unique([j, j1, j2])  # 生成唯一序列
                        len_seq = len(seq)
                        if len_seq == 1:
                            val = 0
                        elif len_seq == 2:
                            index = f_find_index(seq, index_dict2)
                            val = retime2[index]
                        else:
                            index = f_find_index(seq, index_dict3)
                            val = retime3[index]
                        seq = np.unique([m, j1, j2])  # 生成唯一序列
                        if len(seq) == 3:
                            index = f_find_index([j1, j2], index_dict2)
                            tik2 += repro_val[j] * (disc2 - 1) * t33[m, index] / 6 * (trans_mat2[j, m] - trans_mat0[j, m]) * val

        # 计算含有 t33 和 disc2^2 的 tik2
        for j in range(n):
            for m in range(n):
                for j1 in range(n):
                    for j2 in range(n):
                        seq = np.unique([j, m, j1, j2])  # 生成唯一序列
                        len_seq = len(seq)
                        if len_seq == 1:
                            val = 0
                        elif len_seq == 2:
                            index = f_find_index(seq, index_dict2)
                            val = retime2[index]
                        elif len_seq == 3:
                            index = f_find_index(seq, index_dict3)
                            val = retime3[index]
                        else:
                            index = f_find_index(seq, index_dict4)
                            val = retime4[index]
                        seq = np.unique([m, j1, j2])  # 生成唯一序列
                        if len(seq) == 3:
                            index = f_find_index([j1, j2], index_dict2)
                            tik2 += repro_val[j] * (disc2 - 1)**2 * t33[m, index] / 6 * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 计算 bcratio 比率
    bcratio = tik1 / tik2 if tik2 != 0 else 0  # 防止tik2为0导致除零错误

    return bcratio

"""
function bcratio = f_get_bcratio_average_type2(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
trans_mat2 = trans_mat ^ 2;
trans_mat0 = trans_mat ^ 0;
index_dict2 = f_index_dict_two(n);
index_dict3 = f_index_dict_three(n);
index_dict4 = f_index_dict_four(n);
[t12, t13, t22, t23, t33] = f_gen_stru_info_average_type2(madj2, madj3, n);

tik1 = 0;
for j = 1: n
    for m = 1: n
        seq = [j, m];
        seq = unique(seq);
        len = length(seq);
        if len == 1
            val = 0;
        else
            index = f_find_index(seq, index_dict2);
            val = retime2(index);
        end
        tik1 = tik1 + repro_val(j) * (t12(m) + t13(m)) * (trans_mat2(j, m) - trans_mat0(j, m)) * val;
    end
end

tik2 = 0;
for j = 1: n
    for m = 1: n
        seq = [j, m];
        seq = unique(seq);
        len = length(seq);
        if len == 1
            val = 0;
        else
            index = f_find_index(seq, index_dict2);
            val = retime2(index);
        end
        tik2 = tik2 + repro_val(j) * (t12(m)/2 + t13(m)/3) * (trans_mat2(j, m) - trans_mat0(j, m)) * val;
    end
end

for j = 1: n
    for m = 1: n
        for j1 = 1: n
            seq = [j, j1];
            seq = unique(seq);
            len = length(seq);
            if len == 1
                val = 0;
            else
                index = f_find_index(seq, index_dict2);
                val = retime2(index);
            end
            tik2 = tik2 + repro_val(j) * (t22(m, j1)/2 + t23(m, j1)/3) * (trans_mat2(j, m) - trans_mat0(j, m)) * val;
        end
    end
end

for j = 1: n
    for m = 1: n
        for j1 = 1: n
            seq = [j, m, j1];
            seq = unique(seq);
            len = length(seq);
            if len == 1
                val = 0;
            elseif len == 2
                index = f_find_index(seq, index_dict2);
                val = retime2(index);
            else
                index = f_find_index(seq, index_dict3);
                val = retime3(index);
            end
            tik2 = tik2 + repro_val(j) * ((disc1-1)*t22(m, j1)/2 + (disc2-1)*t23(m, j1)/3) * (trans_mat2(j, m) - trans_mat0(j, m)) * val;
        end
    end
end

if sum(sum(madj3)) > 0

for j = 1: n
    for m = 1: n
        for j1 = 1: n
            for j2 = 1: n
                seq = [j, j1, j2];
                seq = unique(seq);
                len = length(seq);
                if len == 1
                    val = 0;
                elseif len == 2
                    index = f_find_index(seq, index_dict2);
                    val = retime2(index);
                else
                    index = f_find_index(seq, index_dict3);
                    val = retime3(index);
                end
                seq = [m, j1, j2];
                seq = unique(seq);
                len = length(seq);
                if len == 3
                    index = f_find_index([j1, j2], index_dict2);
                    tik2 = tik2 + repro_val(j) * (disc2-1)*t33(m, index)/6 * (trans_mat2(j, m) - trans_mat0(j, m)) * val;
                end
            end
        end
    end
end

for j = 1: n
    for m = 1: n
        for j1 = 1: n
            for j2 = 1: n
                seq = [j, m, j1, j2];
                seq = unique(seq);
                len = length(seq);
                if len == 1
                    val = 0;
                elseif len == 2
                    index = f_find_index(seq, index_dict2);
                    val = retime2(index);
                elseif len == 3
                    index = f_find_index(seq, index_dict3);
                    val = retime3(index);
                else
                    index = f_find_index(seq, index_dict4);
                    val = retime4(index);
                end
                seq = [m, j1, j2];
                seq = unique(seq);
                len = length(seq);
                if len == 3
                    index = f_find_index([j1, j2], index_dict2);
                    tik2 = tik2 + repro_val(j) * (disc2-1)^2*t33(m, index)/6 * (trans_mat2(j, m) - trans_mat0(j, m)) * val;
                end
            end
        end
    end
end

end

bcratio = tik1 / tik2;
"""