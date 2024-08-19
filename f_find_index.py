import numpy as np

def f_find_index(seq, index_dict):
    """
    找到seq在index_dict中的序号。

    参数:
        seq (list or numpy array): 需要排序的序列。
        index_dict (numpy array): 包含多个序列的二维数组。

    返回:
        int: 序号，如果找不到则返回-1。
    """
    # 将序列排序
    seq = sorted(seq)
    
    # 获取index_dict的大小
    len_index_dict = index_dict.shape[0]
    
    # 初始化返回值
    serial_num = -1
    
    # 遍历index_dict，查找与排序后的seq匹配的序列
    for i in range(len_index_dict):
        norm_seq = index_dict[i, :]
        # 判断两个序列是否相同
        if np.array_equal(seq, norm_seq):
            serial_num = i  # MATLAB索引从1开始，Python从0开始，所以需要加1
            break
    
    return serial_num

"""
原Matlab代码中涉及到索引的需注意，Python索引从0开始，matlab索引从1开始
function serial_num = f_find_index(seq, index_dict)
seq = sort(seq);
size_record = size(index_dict);
len = size_record(1);
serial_num = -1;
for i = 1: len
    norm_seq = index_dict(i, :);
    if isequal(seq, norm_seq)
        serial_num = i;
        break
    end
end
"""