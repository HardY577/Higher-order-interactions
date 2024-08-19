import numpy as np

from f_baseline_disc1_disc2 import f_baseline_disc1_disc2
from f_cal_remeet_time_four import f_cal_remeet_time_four
from f_cal_remeet_time_three import f_cal_remeet_time_three
from f_cal_remeet_time_two import f_cal_remeet_time_two
from f_gen_conjoin_rich_club import f_gen_conjoin_rich_club
from f_gen_high_mat import f_gen_high_mat
from f_gen_replace_mat import f_gen_replace_mat
from f_gen_trans_mat import f_gen_trans_mat
from f_get_bcratio_accumulate import f_get_bcratio_accumulate

"""
重要更改说明：
1.矩阵初始化：使用 np.zeros_like() 替代 MATLAB 中的 0*madj3 来生成全零矩阵。
2.繁殖值计算：使用 np.sum(replace_mat, axis=0) / np.sum(replace_mat) 计算繁殖值。
3.避免传递空值：对于 f_get_bcratio_accumulate() 中的 retime4 参数，在纯二阶交互中传递 None，以表明不计算四维重新会合时间。

潜在问题：
1.如果 madj2 和 madj3 是稀疏矩阵，并且网络规模较大，考虑使用 scipy.sparse 来优化内存和计算效率。
2.确保 f_baseline_disc1_disc2、f_gen_conjoin_rich_club、f_gen_high_mat、f_gen_raplace_mat、f_gen_trans_mat、f_cal_remeet_time_two、
    f_cal_remeet_time_three、f_cal_remeet_time_four 和 f_get_bcratio_accumulate 函数在 Python 中已经正确实现并经过测试。
"""
# 初始化参数
m1 = 2
n1 = 10
m2 = 2
n2 = 10
n = n1 + n2 + m1 + m2  # 计算总的节点数
disc1 = 1
disc2 = f_baseline_disc1_disc2(disc1)  # 计算 disc2 的值

# 生成连体网络
madj2 = f_gen_conjoin_rich_club(m1, n1, m2, n2)
madj3 = f_gen_high_mat(madj2, n)

# ===== 纯二阶交互 =====
# 生成替代矩阵，madj3 被乘以 0，相当于不考虑第三阶的交互
replace_mat = f_gen_replace_mat(madj2, np.zeros_like(madj3), n)

# 计算繁殖值
repro_val = np.sum(replace_mat, axis=0) / np.sum(replace_mat)

# 生成转移矩阵
trans_mat = f_gen_trans_mat(replace_mat, n)

# 计算两维的重新会合时间
retime2 = f_cal_remeet_time_two(trans_mat, n)

# 计算三维的重新会合时间
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n)

# 计算二阶交互下的 bcratio
bcratio = f_get_bcratio_accumulate(madj2, np.zeros_like(madj3), trans_mat, repro_val, n, retime2, retime3, None, disc1, disc2)
print("Purely pairwise interactions bcratio:", bcratio)

# ===== 高阶交互 (二阶和三阶交互的混合) =====
# 生成新的替代矩阵，考虑到二阶和三阶的交互
replace_mat = f_gen_replace_mat(madj2, madj3, n)

# 重新计算繁殖值
repro_val = np.sum(replace_mat, axis=0) / np.sum(replace_mat)

# 生成新的转移矩阵
trans_mat = f_gen_trans_mat(replace_mat, n)

# 重新计算两维的重新会合时间
retime2 = f_cal_remeet_time_two(trans_mat, n)

# 重新计算三维的重新会合时间
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n)

# 计算四维的重新会合时间
retime4 = f_cal_remeet_time_four(trans_mat, retime2, retime3, n)

# 计算高阶交互下的 bcratio
bcratio = f_get_bcratio_accumulate(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
print("Higher-order interactions bcratio:", bcratio)

"""
m1 = 2; n1 = 10;
m2 = 2; n2 = 10;
n = n1 + n2 + m1 + m2;
disc1 = 1; disc2 = f_baseline_disc1_disc2(disc1);
madj2 = f_gen_conjoin_rich_club(m1, n1, m2, n2);  %% generate conjoined network
madj3 = f_gen_high_mat(madj2, n);

%%%% Purely pairwise interactions
replace_mat = f_gen_raplace_mat(madj2, 0*madj3, n);
repro_val = sum(replace_mat) / sum(sum(replace_mat));
trans_mat = f_gen_trans_mat(replace_mat, n);  
retime2 = f_cal_remeet_time_two(trans_mat, n);    %% calculate two-dimensional coalescence time
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n);    %% calculate three-dimensional coalescence time
bcratio = f_get_bcratio_accumulate(madj2, 0*madj3, trans_mat, repro_val, n, retime2, retime3, 0, disc1, disc2)

%%%% Higher-order interactions (mixture of second-order and third-order interactions)
replace_mat = f_gen_raplace_mat(madj2, madj3, n);
repro_val = sum(replace_mat) / sum(sum(replace_mat));
trans_mat = f_gen_trans_mat(replace_mat, n);
retime2 = f_cal_remeet_time_two(trans_mat, n);
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n);
retime4 = f_cal_remeet_time_four(trans_mat, retime2, retime3, n);    %% calculate four-dimensional coalescence time
bcratio = f_get_bcratio_accumulate(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
"""