import numpy as np

from f_find_index import f_find_index
from f_gen_stru_info_average_type1 import f_gen_stru_info_average_type1
from f_index_dict_four import f_index_dict_four
from f_index_dict_three import f_index_dict_three
from f_index_dict_two import f_index_dict_two

"""
重要更改说明：
1.矩阵幂计算：在Python中，矩阵幂使用np.linalg.matrix_power方法，而在MATLAB中可以直接使用^符号。
2.生成单位矩阵：trans_mat0是单位矩阵，可以用np.eye(n)生成。
3.np.unique：用于生成唯一的组合，这是MATLAB中unique函数的Python替代。
4.注释：详细描述了每个步骤的目的，方便后续优化与理解。

潜在问题：
1.代码中的函数f_find_index、f_index_dict_two等函数没有提供，这些函数的实现可能会影响代码的执行，需要确保它们在Python环境中正确实现。
"""
def f_get_bcratio_average_type1(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2):
    # 计算trans_mat的平方和0次幂
    trans_mat2 = np.linalg.matrix_power(trans_mat, 2)
    trans_mat0 = np.eye(n)  # 相当于trans_mat的0次幂，生成单位矩阵

    # 获取不同组合的索引字典
    index_dict2 = f_index_dict_two(n)
    index_dict3 = f_index_dict_three(n)
    index_dict4 = f_index_dict_four(n)

    # 生成结构信息，使用average_type1的方式
    t12, t13, t22, t23, t33 = f_gen_stru_info_average_type1(madj2, madj3, n)

    tik1 = 0  # 初始化tik1，用于存储加权和

    # 第一个双重循环，计算tik1的值
    for j in range(n):
        for m in range(n):
            seq = np.unique([j, m])  # 生成唯一的组合
            if len(seq) == 1:
                val = 0  # 如果组合长度为1，值设为0
            else:
                index = f_find_index(seq, index_dict2)  # 查找组合在字典中的索引
                val = retime2[index]  # 获取对应的时间值
            # 计算tik1的加权和
            tik1 += repro_val[j] * (t12[m] + t13[m]) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    tik2 = 0  # 初始化tik2

    # 第二个双重循环，计算tik2的初始值
    for j in range(n):
        for m in range(n):
            seq = np.unique([j, m])
            if len(seq) == 1:
                val = 0
            else:
                index = f_find_index(seq, index_dict2)
                val = retime2[index]
            # 计算tik2的加权和，t12和t13有不同的权重
            tik2 += repro_val[j] * (t12[m]/2 + t13[m]/3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 第三个三重循环，进一步累加tik2的值
    for j in range(n):
        for m in range(n):
            for j1 in range(n):
                seq = np.unique([j, j1])
                if len(seq) == 1:
                    val = 0
                else:
                    index = f_find_index(seq, index_dict2)
                    val = retime2[index]
                # 根据t22和t23的值进一步累加tik2
                tik2 += repro_val[j] * (t22[m, j1]/2 + t23[m, j1]/3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 第四个三重循环，包含组合判断
    for j in range(n):
        for m in range(n):
            for j1 in range(n):
                seq = np.unique([j, m, j1])
                if len(seq) == 1:
                    val = 0
                elif len(seq) == 2:
                    index = f_find_index(seq, index_dict2)
                    val = retime2[index]
                else:
                    index = f_find_index(seq, index_dict3)
                    val = retime3[index]
                # 使用disc1和disc2作为不同权重的系数
                tik2 += repro_val[j] * ((disc1-1)*t22[m, j1]/2 + (disc2-1)*t23[m, j1]/3) * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 判断madj3矩阵是否非零，如果非零，进入更复杂的四重循环
    if np.sum(madj3) > 0:
        for j in range(n):
            for m in range(n):
                for j1 in range(n):
                    for j2 in range(n):
                        seq = np.unique([j, j1, j2])
                        if len(seq) == 1:
                            val = 0
                        elif len(seq) == 2:
                            index = f_find_index(seq, index_dict2)
                            val = retime2[index]
                        else:
                            index = f_find_index(seq, index_dict3)
                            val = retime3[index]
                        seq = np.unique([m, j1, j2])
                        if len(seq) == 3:
                            index = f_find_index([j1, j2], index_dict2)
                            # 进一步累加tik2，使用t33矩阵和disc2权重系数
                            tik2 += repro_val[j] * (disc2-1) * t33[m, index]/6 * (trans_mat2[j, m] - trans_mat0[j, m]) * val

        # 再次使用四重循环，计算tik2的最终累加值
        for j in range(n):
            for m in range(n):
                for j1 in range(n):
                    for j2 in range(n):
                        seq = np.unique([j, m, j1, j2])
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
                        seq = np.unique([m, j1, j2])
                        if len(seq) == 3:
                            index = f_find_index([j1, j2], index_dict2)
                            # 最终tik2累加，使用平方的disc2系数
                            tik2 += repro_val[j] * (disc2-1)**2 * t33[m, index]/6 * (trans_mat2[j, m] - trans_mat0[j, m]) * val

    # 计算bcratio的比值
    bcratio = tik1 / tik2
    return bcratio

"""
function bcratio = f_get_bcratio_average_type1(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
trans_mat2 = trans_mat ^ 2;
trans_mat0 = trans_mat ^ 0;
index_dict2 = f_index_dict_two(n);
index_dict3 = f_index_dict_three(n);
index_dict4 = f_index_dict_four(n);
[t12, t13, t22, t23, t33] = f_gen_stru_info_average_type1(madj2, madj3, n);

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