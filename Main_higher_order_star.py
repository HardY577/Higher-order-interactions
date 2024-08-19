import numpy as np

from f_baseline_disc1_disc2 import f_baseline_disc1_disc2
from f_cal_remeet_time_four import f_cal_remeet_time_four
from f_cal_remeet_time_three import f_cal_remeet_time_three
from f_cal_remeet_time_two import f_cal_remeet_time_two
from f_gen_high_mat import f_gen_high_mat
from f_gen_replace_mat import f_gen_replace_mat
from f_gen_star_net import f_gen_star_net
from f_gen_trans_mat import f_gen_trans_mat
from f_get_bcratio_accumulate import f_get_bcratio_accumulate


"""
重要更改说明：
1.矩阵初始化：使用 np.zeros_like() 替代 MATLAB 中的 0*madj3 和 0*madj2 来生成全零矩阵。
2.繁殖值计算：使用 np.sum(replace_mat, axis=0) / np.sum(replace_mat) 计算繁殖值。
3.避免传递空值：在 f_get_bcratio_accumulate() 调用时，对于不涉及四维重新会合时间的情况，传递 None 表示不计算该值。

潜在问题：
1.确保函数 f_baseline_disc1_disc2、f_gen_star_net、f_gen_high_mat、f_gen_raplace_mat、f_gen_trans_mat、f_cal_remeet_time_two、
    f_cal_remeet_time_three、f_cal_remeet_time_four 和 f_get_bcratio_accumulate 均已正确实现并经过测试。
"""
# 初始化参数
leaf_num = 10
n = leaf_num * 2 + 1  # 计算总的节点数
disc1 = 1
disc2 = f_baseline_disc1_disc2(disc1)  # 计算 disc2 的值

# ===== 纯二阶交互 =====
# 生成二阶交互的星型网络结构
madj2 = f_gen_star_net(leaf_num)

# 生成三阶交互结构
madj3 = f_gen_high_mat(madj2, n)

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

# 计算累积收益下的临界比率 bcratio
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

# 计算高阶交互下的临界比率 bcratio
bcratio = f_get_bcratio_accumulate(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
print("Higher-order interactions bcratio:", bcratio)

# ===== 纯三阶交互 =====
# 重新生成二阶交互的星型网络结构
madj2 = f_gen_star_net(leaf_num)

# 重新生成三阶交互结构
madj3 = f_gen_high_mat(madj2, n)

# 生成替代矩阵，只考虑三阶交互，madj2 被乘以 0
replace_mat = f_gen_replace_mat(np.zeros_like(madj2), madj3, n)

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

# 计算累积收益下的临界比率 bcratio
bcratio = f_get_bcratio_accumulate(np.zeros_like(madj2), madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
print("Purely third-order interactions bcratio:", bcratio)

"""
leaf_num = 10;
n = leaf_num * 2 + 1;
disc1 = 1; disc2 = f_baseline_disc1_disc2(disc1);

%%%% Purely pairwise interactions
madj2 = f_gen_star_net(leaf_num);     %% generate pairwise structure
madj3 = f_gen_high_mat(madj2, n);     %% generate third-order structure
replace_mat = f_gen_raplace_mat(madj2, 0*madj3, n);      %% generate replacement graph
repro_val = sum(replace_mat) / sum(sum(replace_mat));    %% calculate reproductive values
trans_mat = f_gen_trans_mat(replace_mat, n);             
retime2 = f_cal_remeet_time_two(trans_mat, n);    %% calculate two-dimensional coalescence time
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n);         %% calculate three-dimensional coalescence time
bcratio = f_get_bcratio_accumulate(madj2, 0*madj3, trans_mat, repro_val, n, retime2, retime3, 0, disc1, disc2)   %% calculate bcratio under accumulated payoff

%%%% Higher-order interactions (mixture of second-order and third-order interactions)
replace_mat = f_gen_raplace_mat(madj2, madj3, n);
repro_val = sum(replace_mat) / sum(sum(replace_mat));
trans_mat = f_gen_trans_mat(replace_mat, n);
retime2 = f_cal_remeet_time_two(trans_mat, n);
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n);
retime4 = f_cal_remeet_time_four(trans_mat, retime2, retime3, n);    %% calculate four-dimensional coalescence time
bcratio = f_get_bcratio_accumulate(madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)

%%%% Purely third-order interactions
madj2 = f_gen_star_net(leaf_num);
madj3 = f_gen_high_mat(madj2, n);
replace_mat = f_gen_raplace_mat(0*madj2, madj3, n);
repro_val = sum(replace_mat) / sum(sum(replace_mat));
trans_mat = f_gen_trans_mat(replace_mat, n);
retime2 = f_cal_remeet_time_two(trans_mat, n);
retime3 = f_cal_remeet_time_three(trans_mat, retime2, n);
retime4 = f_cal_remeet_time_four(trans_mat, retime2, retime3, n);
bcratio = f_get_bcratio_accumulate(0*madj2, madj3, trans_mat, repro_val, n, retime2, retime3, retime4, disc1, disc2)
"""