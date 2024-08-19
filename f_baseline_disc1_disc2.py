import math

def f_baseline_disc1_disc2(disc1):
    """
    基于输入的disc1值计算disc2。
    
    参数:
        disc1 (float): 用于计算的输入值。
    
    返回:
        float: 计算得到的disc2值。
    """
    # 进行计算，相当于MATLAB的sqrt函数
    disc2 = math.sqrt(4 + 3 * disc1) - 1
    return disc2

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function disc2 = f_baseline_disc1_disc2(disc1)
disc2 = (4 + 3 * disc1)^(1/2) - 1;
"""