def f_ana_two(leaf_num, disc1):
    """
    计算f_ana_two的值。
    
    参数:
        leaf_num (float): 叶子数量。
        disc1 (float): 用于计算的输入值。
    
    返回:
        float: 计算的结果值。
    """
    n = leaf_num
    # 计算f1
    f1 = 36 * n * (-7 + 12 * n) * (2 + 23 * n + 14 * n**2)
    # 计算f2
    f2 = 2 * (6 - 29 * n - 1267 * n**2 + 818 * n**3 + 1057 * n**4) + \
         disc1 * (-12 - 338 * n - 1534 * n**2 + 1181 * n**3 + 1288 * n**4)
    # 返回计算结果
    return f1 / f2

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function val = f_ana_two(leaf_num, disc1)
n = leaf_num;
f1 = (36 * n * (-7+12*n) * (2+23*n+14*n^2));
f2 = 2 * (6-29*n-1267*n^2+818*n^3+1057*n^4)+disc1*(-12-338*n-1534*n^2+1181*n^3+1288*n^4);
val = f1 / f2;
"""

