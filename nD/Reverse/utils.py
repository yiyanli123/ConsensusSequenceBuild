import edlib
from func_timeout import func_set_timeout
from dna_features_viewer import GraphicFeature, GraphicRecord
from Levenshtein import distance
import os
import matplotlib.pyplot as plt
import numpy as np 

#每次读取n行文本
def read(filepath, nRows):
    i = 0
    lines = []  # a buffer to cache lines
    with open(filepath,'r',encoding="UTF-8" ) as f:
        for line in f:
            i += 1
            lines.append(line.strip())  # append a line
            if i >= nRows:
                yield lines
                # reset buffer
                i = 0
                lines.clear()
    # remaining lines
    if i > 0:
        yield lines

def check_dir_exist1(res_dir,hairpin_num):
    if not os.path.exists(res_dir + "/hairpin"+str(hairpin_num)+"_detail"):
        os.makedirs(res_dir + "/hairpin"+str(hairpin_num)+"_detail")
        print("%s 目录创建成功"%"/hairpin"+str(hairpin_num)+"_detail")
        
def check_dir_exist2(res_dir,hairpin_num,useful_num):
    if not os.path.exists(res_dir + "/hairpin"+str(hairpin_num)+"_detail" + "/subread_use" + str(useful_num)):
        os.makedirs(res_dir + "/hairpin"+str(hairpin_num)+"_detail"+ "/subread_use" + str(useful_num))
        print("%s 目录创建成功"%"/hairpin"+str(hairpin_num)+"_detail"+ "/subread_use" + str(useful_num))
    #else:
        #print("%s 目录已存在"%"/hairpin"+str(hairpin_num)+"_detail")
        
def DrawGeneFig1(start_list,readID,lens,hairpin,num_id,num_id2,res_dir):  #loc用这个adjusted_locs 
    features = []
    for start in start_list:
        features.append(GraphicFeature(start=start, end=start+len(hairpin), strand=+1, color="#ffd700",label=str(start)))
    record = GraphicRecord(sequence_length=lens, features=features)
    ax, _ = record.plot(figure_width=10)
    dic = res_dir+"/hairpin" + str(num_id)+ "_detail/"+"subread_use" + str(num_id2) + "/"
    PngName = readID +".png"
    ax.figure.savefig(dic + PngName, bbox_inches='tight')
    plt.close(ax.figure)
    
def DrawGeneFig2(start_list,readID,lens,hairpin,num_id,res_dir):  #loc用这个adjusted_locs 
    features = []
    for start in start_list:
        features.append(GraphicFeature(start=start, end=start+len(hairpin), strand=+1, color="#ffd700",label=str(start)))
    record = GraphicRecord(sequence_length=lens, features=features)
    ax, _ = record.plot(figure_width=10)
    dic = res_dir+"/hairpin" + str(num_id)+ "_detail/"
    PngName = readID +"lt3.png"
    ax.figure.savefig(dic + PngName, bbox_inches='tight')
    plt.close(ax.figure)
    
def DNA_complement2(sequence):
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')     # trantab = str.maketrans(intab, outtab)   # 制作翻译表
    string = sequence.translate(trantab)     # str.translate(trantab)  # 转换字符
    return string[::-1]


def Compute_N50(read_len_list):
    read_len_list.sort()
    read_len_list.reverse()
    BaseSum = np.sum(read_len_list)
    N50_pos = BaseSum / 2.0 
    ValueSum = 0
    for value in read_len_list:
        ValueSum += value
        if N50_pos <= ValueSum:
            N50 = value
            break 
    return N50 

