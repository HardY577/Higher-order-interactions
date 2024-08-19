import numpy as np

from f_find_index import f_find_index
from f_gen_stru_info_accumulate import f_gen_stru_info_accumulate
from f_index_dict_four import f_index_dict_four
from f_index_dict_three import f_index_dict_three
from f_index_dict_two import f_index_dict_two


"""
重要更改说明：
1.矩阵幂运算：trans_mat ^ 2 在Python中对应np.linalg.matrix_power(trans_mat, 2)，而trans_mat ^ 0对应单位矩阵np.identity(n)。
2.索引转换：MATLAB中的unique()函数生成唯一的元素集合，Python中使用list(set([...]))。
3.矩阵元素访问：MATLAB使用索引访问矩阵元素，Python对应使用二维数组的下标访问方式。
4.逻辑结构保持一致：所有的for循环和条件判断语句结构与MATLAB一致，确保了代码逻辑不变。

潜在问题：
1.除以零问题：在计算bcratio时，如果tik2为零，则会出现除以零的错误，需要对这种情况进行异常处理。
"""
def f_get_bcratio_accumulate(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2):
    trans_mat2 = np.linalg.matrix_power(trans_mat, 2)
    trans_mat0 = np.identity(n)  # 对应于 trans_mat ^ 0
    index_dict2 = f_index_dict_two(n)
    index_dict3 = f_index_dict_three(n)
    index_dict4 = f_index_dict_four(n)
    
    # 生成结构信息
    t12, t13, t22, t23, t33 = f_gen_stru_info_accumulate(madj2, madj3, n)
    
    tik1 = 0
    for j in range(n):
        for m in range(n):
            seq = list(set([j, m]))
            if len(seq) == 1:
                val = 0
            else:
                index = f_find_index(seq, index_dict2)
                val = retime2[index]
            tik1 += repro_val[j] * (t12[m] + t13[m]) * (trans_mat2[j, m] - trans_mat0[j, m]) * val
    
    tik2 = 0
    for j in range(n):
        for m in range(n):
            seq = list(set([j, m]))
            if len(seq) == 1:
                val = 0
            else:
                index = f_find_index(seq, index_dict2)
                val = retime2[index]
            tik2 += repro_val[j] * (t12[m]/2 + t13[m]/3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val
    
    for j in range(n):
        for m in range(n):
            for j1 in range(n):
                seq = list(set([j, j1]))
                if len(seq) == 1:
                    val = 0
                else:
                    index = f_find_index(seq, index_dict2)
                    val = retime2[index]
                tik2 += repro_val[j] * (t22[m, j1]/2 + t23[m, j1]/3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val
    
    for j in range(n):
        for m in range(n):
            for j1 in range(n):
                seq = list(set([j, m, j1]))
                if len(seq) == 1:
                    val = 0
                elif len(seq) == 2:
                    index = f_find_index(seq, index_dict2)
                    val = retime2[index]
                else:
                    index = f_find_index(seq, index_dict3)
                    val = retime3[index]
                tik2 += repro_val[j] * ((disc1 - 1) * t22[m, j1]/2 + (disc2 - 1) * t23[m, j1]/3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val
    
    if np.sum(madj3) > 0:
        for j in range(n):
            for m in range(n):
                for j1 in range(n):
                    for j2 in range(n):
                        seq = list(set([j, j1, j2]))
                        if len(seq) == 1:
                            val = 0
                        elif len(seq) == 2:
                            index = f_find_index(seq, index_dict2)
                            val = retime2[index]
                        else:
                            index = f_find_index(seq, index_dict3)
                            val = retime3[index]
                        seq_m = list(set([m, j1, j2]))
                        if len(seq_m) == 3:
                            index_m = f_find_index([j1, j2], index_dict2)
                            tik2 += repro_val[j] * (disc2 - 1) * t33[m, index_m]/6 * (trans_mat2[j, m] - trans_mat0[j, m]) * val

        for j in range(n):
            for m in range(n):
                for j1 in range(n):
                    for j2 in range(n):
                        seq = list(set([j, m, j1, j2]))
                        if len(seq) == 1:
                            val = 0
                        elif len(seq) == 2:
                            index = f_find_index(seq, index_dict2)
                            val = retime2[index]
                        elif len(seq) == 3:
                            index = f_find_index(seq, index_dict3)
                            val = retime3[index]
                        else:
                            index = f_find_index(seq, index_dict4)
                            val = retime4[index]
                        seq_m = list(set([m, j1, j2]))
                        if len(seq_m) == 3:
                            index_m = f_find_index([j1, j2], index_dict2)
                            tik2 += repro_val[j] * (disc2 - 1)**2 * t33[m, index_m]/6 * (trans_mat2[j, m] - trans_mat0[j, m]) * val
    
    bcratio = tik1 / tik2
    return bcratio

"""
function bcratio = f_get_bcratio_accumulate(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
trans_mat2 = trans_mat ^ 2;
trans_mat0 = trans_mat ^ 0;
index_dict2 = f_index_dict_two(n);
index_dict3 = f_index_dict_three(n);
index_dict4 = f_index_dict_four(n);
[t12, t13, t22, t23, t33] = f_gen_stru_info_accumulate(madj2, madj3, n);

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