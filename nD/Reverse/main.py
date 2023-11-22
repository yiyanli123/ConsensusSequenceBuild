import pandas as pd 
import numpy as np
from Levenshtein import distance
from math import ceil
from random import choices, randint
import random
from dna_features_viewer import GraphicFeature, GraphicRecord
import os
import pyalign
from collections import Counter
from args_parser import set_parser
from func_timeout import func_set_timeout
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import edlib
from utils import read, check_dir_exist1, check_dir_exist2, DrawGeneFig1,+ DNA_complement2, Compute_N50

# 导入time模块
import time
# 程序开始时间
begin_time = time.time()


args = set_parser()
res_dir = args.res_dir
res =  args.fq_file
subread_len_filter_thre = args.subread_len_filter_thre

hairpin = args.hairpin_seq
windows_len = len(hairpin)
step_len = args.step_len
LD_thre = np.ceil(windows_len*args.err_threshold)
adjust_dis = args.adjusted_len
hairpin_len = len(args.hairpin_seq)


@func_set_timeout(5)
def test_edlib(seq1,seq2):
    test_res = edlib.align(seq1,seq2, mode="global", task="path")
    res_ld = 1- (test_res["editDistance"]/len(seq1))
    return np.round(res_ld,3)

if not os.path.exists(res_dir):
    os.makedirs(res_dir)
    print("%s 目录创建成功"%res_dir)
else:
    print("%s 目录已经存在"%res_dir)

raw_hairpin_distribution_dic = {}  #放到read函数循环体外面

raw_read_len_dict = {}             #记录原始check到hairpin的read长度  表格第一部分
raw_subreads_1_10_len = {}          #记录subread1-8各自的subread长度   表格第四部分

subread_filter_2_read_dict = {}    #记录第二次筛选过后的保留下来的read长度，表格第二部分

first_filt_dic = {}
second_filt_dic = {}


check_read_count = 0   #记录check了多少条序列
lines_gen = read(res, 4)
read_num = 1
draw_count = 0
draw_count2 = 0
for read_ele in lines_gen:
    check_read_count += 1
    locs = []
    for i in range((np.ceil(len(read_ele[1]) - windows_len)/step_len).astype("int")):
        ld = distance(read_ele[1][i*step_len:i*step_len+len(hairpin)], hairpin)
        if ld <=LD_thre :
            locs.append(i*step_len)

    adjusted_locs = []
    if locs == [] :   #没有hairpin被识别到
        adjusted_hairpin_num = 0
        if adjusted_hairpin_num not in raw_hairpin_distribution_dic:
            raw_hairpin_distribution_dic[adjusted_hairpin_num] = 1
        else:
            raw_hairpin_distribution_dic[adjusted_hairpin_num] += 1
    else:
        adjusted_locs.append(locs[0])
        for j in locs:
            if j <= adjusted_locs[-1] + adjust_dis:
                continue
            else:
                adjusted_locs.append(j)

    if adjusted_locs != [] :
        adjusted_hairpin_num = len(adjusted_locs)
        if adjusted_hairpin_num not in raw_hairpin_distribution_dic:
            raw_hairpin_distribution_dic[adjusted_hairpin_num] = 1
        else:
            raw_hairpin_distribution_dic[adjusted_hairpin_num] += 1

        if adjusted_hairpin_num not in raw_read_len_dict:   #记录原始read的长度
            raw_read_len_dict[adjusted_hairpin_num] = []
            raw_read_len_dict[adjusted_hairpin_num].append(len(read_ele[1]))
        else:
            raw_read_len_dict[adjusted_hairpin_num].append(len(read_ele[1]))
    #print(adjusted_locs)

    ##输出原始识别到的subraed到subread1~8的文件中，对应输出表格的第四部分。
    #判断第一个hairpin的位置
    
    if adjusted_locs != [] :
        if adjusted_locs[0] >0:                           #原始识别到的subread1>0
            subreads_len_list_raw = []                    #记录每个subread的长度
            subreads_list_raw = []                        #subread具体序列  

            for ll in range(len(adjusted_locs)):
                if ll == 0 :
                    subreads_len_list_raw.append(adjusted_locs[0] - 0)
                    subreads_list_raw.append(read_ele[1][:adjusted_locs[0]])  #添加具体序列信息
                elif 0< ll < len(adjusted_locs)-1 :
                    subreads_len_list_raw.append(adjusted_locs[ll] - (adjusted_locs[ll-1]+hairpin_len))
                    subreads_list_raw.append(read_ele[1][adjusted_locs[ll-1]+hairpin_len:adjusted_locs[ll]])
                else:
                    subreads_len_list_raw.append(adjusted_locs[ll] - (adjusted_locs[ll-1]+hairpin_len))
                    subreads_list_raw.append(read_ele[1][adjusted_locs[ll-1]+hairpin_len:adjusted_locs[ll]])
                    subreads_len_list_raw.append(len(read_ele[1])-(adjusted_locs[ll]+hairpin_len))
                    subreads_list_raw.append(read_ele[1][adjusted_locs[ll]+hairpin_len:])

        else:
            adjusted_locs_new = adjusted_locs[1:]
            subreads_len_list_raw = []                    #记录每个subread的长度
            subreads_list_raw = []                        #subread具体序列

            for ll in range(len(adjusted_locs_new)):
                if ll == 0 :
                    subreads_len_list_raw.append(adjusted_locs_new[0] - 0)
                    subreads_list_raw.append(read_ele[1][:adjusted_locs_new[0]])  #添加具体序列信息
                elif 0< ll < len(adjusted_locs_new)-1 :
                    subreads_len_list_raw.append(adjusted_locs_new[ll] - (adjusted_locs_new[ll-1]+hairpin_len))
                    subreads_list_raw.append(read_ele[1][adjusted_locs_new[ll-1]+hairpin_len:adjusted_locs_new[ll]])
                else:
                    subreads_len_list_raw.append(adjusted_locs_new[ll] - (adjusted_locs_new[ll-1]+hairpin_len))
                    subreads_list_raw.append(read_ele[1][adjusted_locs_new[ll-1]+hairpin_len:adjusted_locs_new[ll]])
                    subreads_len_list_raw.append(len(read_ele[1])-(adjusted_locs_new[ll]+hairpin_len))
                    subreads_list_raw.append(read_ele[1][adjusted_locs_new[ll]+hairpin_len:])

        #把的subreadx输出到.fa文件中,记录每个subread的长度
        for subread_ele in range(len(subreads_list_raw)):
            if subread_ele+1 not in raw_subreads_1_10_len:
                raw_subreads_1_10_len[subread_ele+1] = []
                raw_subreads_1_10_len[subread_ele+1].append(len(subreads_list_raw[subread_ele]))
            else:
                raw_subreads_1_10_len[subread_ele+1].append(len(subreads_list_raw[subread_ele]))
            
            if subread_ele <= 9:
                with open(res_dir + "/" + "subread" + str(subread_ele+1) + ".fa" , "a",encoding="utf-8") as f1:
                    f1.write(">" + read_ele[0][1:] + "_subread" + str(subread_ele + 1))
                    f1.write("\n")
                    f1.write(subreads_list_raw[subread_ele])
                    f1.write("\n")
                f1.close()


    #print(len(adjusted_locs))
    if 9>=adjusted_hairpin_num>=3:
        if adjusted_locs[0] >0 :#原始识别到的hairpin数目
            subreads_len_list_tmp = []
            subreads_list_tmp = []
            for m in range(len(adjusted_locs)):
                if m == 0:
                    subreads_len_list_tmp.append(adjusted_locs[0] - 0)
                    subreads_list_tmp.append(read_ele[1][:adjusted_locs[0]])  #添加具体序列信息  #之前有bug
                elif 0<m< len(adjusted_locs)-1 :
                    subreads_len_list_tmp.append(adjusted_locs[m] - (adjusted_locs[m-1]+hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs[m-1]+hairpin_len:adjusted_locs[m]])
                else :
                    subreads_len_list_tmp.append(adjusted_locs[m] - (adjusted_locs[m-1]+hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs[m-1]+hairpin_len:adjusted_locs[m]])
                    subreads_len_list_tmp.append(len(read_ele[1])-(adjusted_locs[m]+hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs[m]+hairpin_len:])
            if len(subreads_list_tmp[0]) < 0.5*(len(subreads_list_tmp[1])):  #第一个的hairpin长度<第二个hairpin的50%
                subreads_list_tmp_1 = subreads_list_tmp[1:]
                subread_idx = [xx for xx in range(len(subreads_list_tmp_1))]
                final_subread_idx = []
                for j in subread_idx:
                    if 0.8*(len(subreads_list_tmp_1[1])) <= len(subreads_list_tmp_1[j]) <= 1.2*(len(subreads_list_tmp_1[1])):
                        final_subread_idx.append(j)
                subreads_list_tmp_1_compli = []
                for i in range(len(subreads_list_tmp_1)):
                    if i % 2 == 0:
                        subreads_list_tmp_1_compli.append(subreads_list_tmp_1[i])
                    else:
                        subreads_list_tmp_1_compli.append(DNA_complement2(subreads_list_tmp_1[i]))
                        
                subreads_list_tmp_1 = subreads_list_tmp_1_compli
                
                filter_first_hairpin_num = len(final_subread_idx)-1
                check_dir_exist1(res_dir,filter_first_hairpin_num) 
                
                if filter_first_hairpin_num >=2:
                    if filter_first_hairpin_num not in first_filt_dic:
                        first_filt_dic[filter_first_hairpin_num] = 1
                    else:
                        first_filt_dic[filter_first_hairpin_num] += 1
                
                #获取经过第一次过滤后的subread
                subread_filtered1_list = [subreads_list_tmp_1[jj] for jj in final_subread_idx]
                if len(subread_filtered1_list) >=2:
                    subread_filtered2_idx = []
                    for ff in final_subread_idx:
                        try:
                            res_Identity = test_edlib(subreads_list_tmp_1[1], subreads_list_tmp_1[ff])
                        except:
                            res_Identity = 0
                            print("Time error"+read_ele[0])
                        if res_Identity > 0.8:
                            subread_filtered2_idx.append(ff) #记录经过第二次过滤的subread的绝对idx，而不是相对idx
                            #subread_filtered2.append(subread_filtered1_list[ff])#记录经过第二次过滤的subread的序列信息，长度和subread_filtered2_idx一样
                        del res_Identity
                    filter_second_subread_num = len(subread_filtered2_idx)
                    
                    if filter_second_subread_num >=3:
                        if filter_second_subread_num not in second_filt_dic:
                            second_filt_dic[filter_second_subread_num] = 1
                        else:
                            second_filt_dic[filter_second_subread_num] += 1

                        if filter_second_subread_num not in subread_filter_2_read_dict:
                            subread_filter_2_read_dict[filter_second_subread_num] = []
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))
                        else:
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))

                        if filter_second_subread_num >=3:

                            check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                            with open(res_dir+"/hairpin"+ str(filter_first_hairpin_num) +"_detail/"+ "subread_use" + str(filter_second_subread_num)+ "/" +  read_ele[0][1:], "w", encoding="utf-8") as f:
                                for kk in subread_filtered2_idx:
                                    f.write(">" + str(kk+1))
                                    f.write("\n")
                                    f.write(subreads_list_tmp_1[kk])
                                    f.write("\n")
                                    #f.flush()
                            f.close()
                            if draw_count <= 200:
                                DrawGeneFig1(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,filter_second_subread_num,res_dir)
                                draw_count += 1
                    else: #如果最终的subread数目≤3
                        #check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                        if draw_count2 <= 200:
                            DrawGeneFig2(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,res_dir)
                            draw_count2 += 1
                            
            else: #第一个的hairpin长度>第二个hairpin的50%
                subreads_list_tmp_1 = subreads_list_tmp
                subread_idx = [xx for xx in range(len(subreads_list_tmp_1))]
                final_subread_idx = []
                for j in subread_idx:
                    if 0.8*(len(subreads_list_tmp_1[0])) <= len(subreads_list_tmp_1[j]) <= 1.2*(len(subreads_list_tmp_1[0])):
                        final_subread_idx.append(j)
                subreads_list_tmp_1_compli = []
                for i in range(len(subreads_list_tmp_1)):
                    if i % 2 == 0:
                        subreads_list_tmp_1_compli.append(subreads_list_tmp_1[i])
                    else:
                        subreads_list_tmp_1_compli.append(DNA_complement2(subreads_list_tmp_1[i]))
                        
                subreads_list_tmp_1 = subreads_list_tmp_1_compli
                
                filter_first_hairpin_num = len(final_subread_idx)-1
                check_dir_exist1(res_dir,filter_first_hairpin_num)
                #获取经过第一次过滤后的subread
                subread_filtered1_list = [subreads_list_tmp_1[jj] for jj in final_subread_idx]
                
                if filter_first_hairpin_num >=2:
                        if filter_first_hairpin_num not in first_filt_dic:
                            first_filt_dic[filter_first_hairpin_num] = 1
                        else:
                            first_filt_dic[filter_first_hairpin_num] += 1
                
                if len(subread_filtered1_list) >=2:
                    subread_filtered2_idx = []
                    for ff in final_subread_idx:
                        try:
                            res_Identity = test_edlib(subreads_list_tmp_1[1], subreads_list_tmp_1[ff])
                        except:
                            res_Identity = 0
                            print("Time error"+read_ele[0])
                        if res_Identity > 0.8:
                            subread_filtered2_idx.append(ff) #记录经过第二次过滤的subread的绝对idx，而不是相对idx
                            #subread_filtered2.append(subread_filtered1_list[ff])#记录经过第二次过滤的subread的序列信息，长度和subread_filtered2_idx一样
                        del res_Identity
                    filter_second_subread_num = len(subread_filtered2_idx)
                    
                    if filter_second_subread_num >=3:
                        if filter_second_subread_num not in second_filt_dic:
                            second_filt_dic[filter_second_subread_num] = 1
                        else:
                            second_filt_dic[filter_second_subread_num] += 1

                        if filter_second_subread_num not in subread_filter_2_read_dict:
                            subread_filter_2_read_dict[filter_second_subread_num] = []
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))
                        else:
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))

                        if filter_second_subread_num >=3:
                            check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                            with open(res_dir+"/hairpin"+ str(filter_first_hairpin_num) +"_detail/"+ "subread_use" + str(filter_second_subread_num)+ "/" +  read_ele[0][1:], "w", encoding="utf-8") as f:
                                for kk in subread_filtered2_idx:
                                    f.write(">" + str(kk+1))
                                    f.write("\n")
                                    f.write(subreads_list_tmp_1[kk])
                                    f.write("\n")
                                    #f.flush()
                            f.close()
                            if draw_count <= 200:
                                DrawGeneFig1(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,filter_second_subread_num,res_dir)
                                draw_count += 1
                    else: #如果最终的subread数目≤3
                        #check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                        if draw_count2 <= 200:
                            DrawGeneFig2(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,res_dir)
                            draw_count2 += 1
                
        else:  #第一个hairpin位置<0
            adjusted_locs_new = adjusted_locs[1:]
            subreads_len_list_tmp = []
            subreads_list_tmp = []
            for m in range(len(adjusted_locs_new)):
                if m == 0:
                    subreads_len_list_tmp.append(adjusted_locs_new[0] - (adjusted_locs[0] + hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs[0] + hairpin_len:adjusted_locs_new[0]])  #添加具体序列信息
                elif 0<m< len(adjusted_locs_new)-1 :
                    subreads_len_list_tmp.append(adjusted_locs_new[m] - (adjusted_locs_new[m-1]+hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs_new[m-1]+hairpin_len:adjusted_locs_new[m]])
                else:
                    subreads_len_list_tmp.append(adjusted_locs_new[m] - (adjusted_locs_new[m-1]+hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs_new[m-1]+hairpin_len:adjusted_locs_new[m]])
                    subreads_len_list_tmp.append(len(read_ele[1])-(adjusted_locs_new[m]+hairpin_len))
                    subreads_list_tmp.append(read_ele[1][adjusted_locs_new[m]+hairpin_len:])
            #print(subreads_list_tmp)
            if len(subreads_list_tmp[0]) < 0.5*(len(subreads_list_tmp[1])):  #第一个的hairpin长度<第二个hairpin的50%
                subreads_list_tmp_1 = subreads_list_tmp[1:]
                subread_idx = [xx for xx in range(len(subreads_list_tmp_1))]
                final_subread_idx = []
                for j in subread_idx:
                    if 0.8*(len(subreads_list_tmp_1[1])) <= len(subreads_list_tmp_1[j]) <= 1.2*(len(subreads_list_tmp_1[1])):
                        final_subread_idx.append(j)
                subreads_list_tmp_1_compli = []
                for i in range(len(subreads_list_tmp_1)):
                    if i % 2 == 0:
                        subreads_list_tmp_1_compli.append(subreads_list_tmp_1[i])
                    else:
                        subreads_list_tmp_1_compli.append(DNA_complement2(subreads_list_tmp_1[i]))

                subreads_list_tmp_1 = subreads_list_tmp_1_compli

                filter_first_hairpin_num = len(final_subread_idx)-1
                check_dir_exist1(res_dir,filter_first_hairpin_num)
                
                if filter_first_hairpin_num >=2:
                        if filter_first_hairpin_num not in first_filt_dic:
                            first_filt_dic[filter_first_hairpin_num] = 1
                        else:
                            first_filt_dic[filter_first_hairpin_num] += 1

                #获取经过第一次过滤后的subread
                subread_filtered1_list = [subreads_list_tmp_1[jj] for jj in final_subread_idx]
                if len(subread_filtered1_list) >=2:
                    subread_filtered2_idx = []
                    for ff in final_subread_idx:
                        try:
                            res_Identity = test_edlib(subreads_list_tmp_1[1], subreads_list_tmp_1[ff])
                        except:
                            res_Identity = 0
                            print("Time error"+read_ele[0])
                        if res_Identity > 0.8:
                            subread_filtered2_idx.append(ff) #记录经过第二次过滤的subread的绝对idx，而不是相对idx
                            #subread_filtered2.append(subread_filtered1_list[ff])#记录经过第二次过滤的subread的序列信息，长度和subread_filtered2_idx一样
                        del res_Identity
                    filter_second_subread_num = len(subread_filtered2_idx)

                    if filter_second_subread_num >=3:
                        if filter_second_subread_num not in second_filt_dic:
                            second_filt_dic[filter_second_subread_num] = 1
                        else:
                            second_filt_dic[filter_second_subread_num] += 1

                        if filter_second_subread_num not in subread_filter_2_read_dict:
                            subread_filter_2_read_dict[filter_second_subread_num] = []
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))
                        else:
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))

                        if filter_second_subread_num >=3:
                            check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                            with open(res_dir+"/hairpin"+ str(filter_first_hairpin_num) +"_detail/"+ "subread_use" + str(filter_second_subread_num)+ "/" +  read_ele[0][1:], "w", encoding="utf-8") as f:
                                for kk in subread_filtered2_idx:
                                    f.write(">" + str(kk+1))
                                    f.write("\n")
                                    f.write(subreads_list_tmp_1[kk])
                                    f.write("\n")
                                    #f.flush()
                            f.close()
                            if draw_count <= 200:
                                DrawGeneFig1(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,filter_second_subread_num,res_dir)
                                draw_count += 1
                    else: #如果最终的subread数目≤3
                        #check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                        if draw_count2 <= 200:
                            DrawGeneFig2(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,res_dir)
                            draw_count2 += 1


            else:#第一个的hairpin长度>第二个hairpin的50%
                subreads_list_tmp_1 = subreads_list_tmp
                subread_idx = [xx for xx in range(len(subreads_list_tmp_1))]
                final_subread_idx = []
                for j in subread_idx:
                    if 0.8*(len(subreads_list_tmp_1[1])) <= len(subreads_list_tmp_1[j]) <= 1.2*(len(subreads_list_tmp_1[1])):
                        final_subread_idx.append(j)
                subreads_list_tmp_1_compli = []
                for i in range(len(subreads_list_tmp_1)):
                    if i % 2 == 0:
                        subreads_list_tmp_1_compli.append(subreads_list_tmp_1[i])
                    else:
                        subreads_list_tmp_1_compli.append(DNA_complement2(subreads_list_tmp_1[i]))

                subreads_list_tmp_1 = subreads_list_tmp_1_compli

                filter_first_hairpin_num = len(final_subread_idx)-1
                check_dir_exist1(res_dir,filter_first_hairpin_num)
                #获取经过第一次过滤后的subread
                subread_filtered1_list = [subreads_list_tmp_1[jj] for jj in final_subread_idx]
                
                if filter_first_hairpin_num >=2:
                        if filter_first_hairpin_num not in first_filt_dic:
                            first_filt_dic[filter_first_hairpin_num] = 1
                        else:
                            first_filt_dic[filter_first_hairpin_num] += 1

                if len(subread_filtered1_list) >=2:
                    subread_filtered2_idx = []
                    for ff in final_subread_idx:
                        try:
                            res_Identity = test_edlib(subreads_list_tmp_1[1], subreads_list_tmp_1[ff])
                        except:
                            res_Identity = 0
                            print("Time error"+read_ele[0])
                        if res_Identity > 0.8:
                            subread_filtered2_idx.append(ff) #记录经过第二次过滤的subread的绝对idx，而不是相对idx
                            #subread_filtered2.append(subread_filtered1_list[ff])#记录经过第二次过滤的subread的序列信息，长度和subread_filtered2_idx一样
                        del res_Identity
                    filter_second_subread_num = len(subread_filtered2_idx)

                    if filter_second_subread_num >=3:
                        if filter_second_subread_num not in second_filt_dic:
                            second_filt_dic[filter_second_subread_num] = 1
                        else:
                            second_filt_dic[filter_second_subread_num] += 1

                        #记录subread筛选后的read长度
                        if filter_second_subread_num not in subread_filter_2_read_dict:
                            subread_filter_2_read_dict[filter_second_subread_num] = []
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))
                        else:
                            subread_filter_2_read_dict[filter_second_subread_num].append(len(read_ele[1]))

                        if filter_second_subread_num >=3:
                            check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                            with open(res_dir+"/hairpin"+ str(filter_first_hairpin_num) +"_detail/"+ "subread_use" + str(filter_second_subread_num)+ "/" +  read_ele[0][1:], "w", encoding="utf-8") as f:
                                for kk in subread_filtered2_idx:
                                    f.write(">" + str(kk+1))
                                    f.write("\n")
                                    f.write(subreads_list_tmp_1[kk])
                                    f.write("\n")
                                    #f.flush()
                            f.close()
                            if draw_count <= 200:
                                DrawGeneFig1(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,filter_second_subread_num,res_dir)
                                draw_count += 1
                    else: #如果最终的subread数目≤3
                        #check_dir_exist2(res_dir,filter_first_hairpin_num,filter_second_subread_num)
                        if draw_count2 <= 200:
                            DrawGeneFig2(adjusted_locs,read_ele[0][1:],len(read_ele[1]),hairpin,filter_first_hairpin_num,res_dir)
                            draw_count2 += 1
                    

#统计subread识别后的read数据量、长度等信息
#base数据量
read_gt3 = 0
read_gt4 = 0
read_gt5 = 0
read_gt6 = 0
read_gt7 = 0 
read_gt8 = 0 
read_gt9 = 0 
read_gt10 = 0 

read_gt3_read_num = 0
read_gt4_read_num = 0
read_gt5_read_num = 0
read_gt6_read_num = 0
read_gt7_read_num = 0 
read_gt8_read_num = 0
read_gt9_read_num = 0
read_gt10_read_num = 0

read_gt3_len_list = []
read_gt4_len_list = []
read_gt5_len_list = []
read_gt6_len_list = []
read_gt7_len_list = []
read_gt8_len_list = []
read_gt9_len_list = []
read_gt10_len_list = []

for i,j in raw_read_len_dict.items():
    if i >= 2:
        read_gt3 += np.sum(j)
        read_gt3_read_num += len(j)
        #read_gt4_len.append(np.round(np.mean(j),3))
        read_gt3_len_list.extend(j)                    #计算N50

for i,j in raw_read_len_dict.items():
    if i >= 3:
        read_gt4 += np.sum(j)
        read_gt4_read_num += len(j)
        #read_gt4_len.append(np.round(np.mean(j),3))
        read_gt4_len_list.extend(j)                    #计算N50

for i,j in raw_read_len_dict.items():
    if i >= 4:
        read_gt5 += np.sum(j)
        read_gt5_read_num += len(j)
        #read_gt5_len.append(np.round(np.mean(j)))
        read_gt5_len_list.extend(j)                    #计算N50

for i,j in raw_read_len_dict.items():
    if i >= 5:
        read_gt6 += np.sum(j)
        read_gt6_read_num += len(j)
        #read_gt6_len.append(np.round(np.mean(j)))
        read_gt6_len_list.extend(j)                    #计算N50

for i,j in raw_read_len_dict.items():
    if i >= 6:
        read_gt7 += np.sum(j)
        read_gt7_read_num += len(j)
        #read_gt7_len.append(np.round(np.mean(j)))
        read_gt7_len_list.extend(j)                    #计算N50

for i,j in raw_read_len_dict.items():
    if i >= 7:
        read_gt8 += np.sum(j)
        read_gt8_read_num += len(j)
        #read_gt8_len.append(np.round(np.mean(j)))
        read_gt8_len_list.extend(j)                    #计算N50
        
for i,j in raw_read_len_dict.items():
    if i >= 8:
        read_gt9 += np.sum(j)
        read_gt9_read_num += len(j)
        #read_gt8_len.append(np.round(np.mean(j)))
        read_gt9_len_list.extend(j)
for i,j in raw_read_len_dict.items():
    if i >= 9:
        read_gt10 += np.sum(j)
        read_gt10_read_num += len(j)
        #read_gt8_len.append(np.round(np.mean(j)))
        read_gt10_len_list.extend(j)

        
        
#计算N50
read_gt3_len_N50 = Compute_N50(read_gt3_len_list)
read_gt4_len_N50 = Compute_N50(read_gt4_len_list)
read_gt5_len_N50 = Compute_N50(read_gt5_len_list)
read_gt6_len_N50 = Compute_N50(read_gt6_len_list)
read_gt7_len_N50 = Compute_N50(read_gt7_len_list)
read_gt8_len_N50 = Compute_N50(read_gt8_len_list)
read_gt9_len_N50 = Compute_N50(read_gt9_len_list)
read_gt10_len_N50 = Compute_N50(read_gt10_len_list)

print("---------------------------------------------------")

print("Subread识别后subread≥3的Read Num为: " +str(read_gt3_read_num))
print("Subread识别后subread≥4的Read Num为: " +str(read_gt4_read_num))
print("Subread识别后subread≥5的Read Num为: " +str(read_gt5_read_num))
print("Subread识别后subread≥6的Read Num为: " +str(read_gt6_read_num))
print("Subread识别后subread≥7的Read Num为: " +str(read_gt7_read_num))
print("Subread识别后subread≥8的Read Num为: " +str(read_gt8_read_num))
print("Subread识别后subread≥9的Read Num为: " +str(read_gt9_read_num))
print("Subread识别后subread≥10的Read Num为: " +str(read_gt10_read_num))

print("---------------------------------------------------")

print("Subread识别后subread≥3的Base Num为: " +str(read_gt3))
print("Subread识别后subread≥4的Base Num为: " +str(read_gt4))
print("Subread识别后subread≥5的Base Num为: " +str(read_gt5))
print("Subread识别后subread≥6的Base Num为: " +str(read_gt6))
print("Subread识别后subread≥7的Base Num为: " +str(read_gt7))
print("Subread识别后subread≥8的Base Num为: " +str(read_gt8))
print("Subread识别后subread≥7的Base Num为: " +str(read_gt9))
print("Subread识别后subread≥8的Base Num为: " +str(read_gt10))

print("---------------------------------------------------")

print("Subread识别后subread≥3的Read AVG Length为: " +str(np.round(np.mean(read_gt3_len_list))))
print("Subread识别后subread≥4的Read AVG Length为: " +str(np.round(np.mean(read_gt4_len_list))))
print("Subread识别后subread≥5的Read AVG Length为: " +str(np.round(np.mean(read_gt5_len_list))))
print("Subread识别后subread≥6的Read AVG Length为: " +str(np.round(np.mean(read_gt6_len_list))))
print("Subread识别后subread≥7的Read AVG Length为: " +str(np.round(np.mean(read_gt7_len_list))))
print("Subread识别后subread≥8的Read AVG Length为: " +str(np.round(np.mean(read_gt8_len_list))))
print("Subread识别后subread≥7的Read AVG Length为: " +str(np.round(np.mean(read_gt9_len_list))))
print("Subread识别后subread≥8的Read AVG Length为: " +str(np.round(np.mean(read_gt10_len_list))))

print("---------------------------------------------------")

print("Subread识别后subread≥3的Read N50为: " +str(read_gt3_len_N50))
print("Subread识别后subread≥4的Read N50为: " +str(read_gt4_len_N50))
print("Subread识别后subread≥5的Read N50为: " +str(read_gt5_len_N50))
print("Subread识别后subread≥6的Read N50为: " +str(read_gt6_len_N50))
print("Subread识别后subread≥7的Read N50为: " +str(read_gt7_len_N50))
print("Subread识别后subread≥8的Read N50为: " +str(read_gt8_len_N50))
print("Subread识别后subread≥7的Read N50为: " +str(read_gt9_len_N50))
print("Subread识别后subread≥8的Read N50为: " +str(read_gt10_len_N50))

print("---------------------------------------------------")
#打印raw data的subread长度

for subread_idx, subread_len_tmp in  raw_subreads_1_10_len.items():
    if subread_idx <= 10:
        print("第" + str(subread_idx) +"条subread的平均长度为: " + str(np.round(np.mean(subread_len_tmp),3)))


#把所有subread的长度放到一起
total_raw_subread_len = []
for i in raw_subreads_1_10_len.values():
    total_raw_subread_len.extend(i)

print("总体subread的平均长度为: " + str(np.round(np.mean(total_raw_subread_len),3)))

print("---------------------------------------------------")


#统计subread筛选后的read数据量、长度等信息
#base数据量
read_gt3_2 = 0
read_gt4_2 = 0
read_gt5_2 = 0
read_gt6_2 = 0
read_gt7_2 = 0 
read_gt8_2 = 0 
read_gt9_2 = 0 
read_gt10_2 = 0 

read_gt3_read_num_2 = 0
read_gt4_read_num_2 = 0
read_gt5_read_num_2 = 0
read_gt6_read_num_2 = 0
read_gt7_read_num_2 = 0
read_gt8_read_num_2 = 0
read_gt9_read_num_2 = 0
read_gt10_read_num_2 = 0

#read长度
read_gt3_len_2 = []
read_gt4_len_2 = []
read_gt5_len_2 = []
read_gt6_len_2 = []
read_gt7_len_2 = []
read_gt8_len_2 = []
read_gt9_len_2 = []
read_gt10_len_2 = []

read_gt3_len_list_2 = []
read_gt4_len_list_2 = []
read_gt5_len_list_2 = []
read_gt6_len_list_2 = []
read_gt7_len_list_2 = []
read_gt8_len_list_2 = []
read_gt9_len_list_2 = []
read_gt10_len_list_2 = []

if 4 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 3:
            read_gt3_2 += np.sum(j)                          #计算base
            read_gt3_read_num_2 += len(j)                    #计算read数目
            read_gt3_len_2.append(np.round(np.mean(j)))      #计算平均长度
            read_gt3_len_list_2.extend(j)                    #计算N50
    read_gt3_len_N50_2 = Compute_N50(read_gt3_len_list_2)
    
if 4 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 4:
            read_gt4_2 += np.sum(j)                          #计算base
            read_gt4_read_num_2 += len(j)                    #计算read数目
            read_gt4_len_2.append(np.round(np.mean(j)))      #计算平均长度
            read_gt4_len_list_2.extend(j)                    #计算N50
    read_gt4_len_N50_2 = Compute_N50(read_gt4_len_list_2)
    
if 5 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 5:
            read_gt5_2 += np.sum(j)
            read_gt5_read_num_2 += len(j)
            read_gt5_len_2.append(np.round(np.mean(j)))
            read_gt5_len_list_2.extend(j)                    #计算N50
    read_gt5_len_N50_2 = Compute_N50(read_gt5_len_list_2)
    
if 6 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 6:
            read_gt6_2 += np.sum(j)
            read_gt6_read_num_2 += len(j)
            read_gt6_len_2.append(np.round(np.mean(j)))
            read_gt6_len_list_2.extend(j)                    #计算N50
    read_gt6_len_N50_2 = Compute_N50(read_gt6_len_list_2)

if 7 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 7:
            read_gt7_2 += np.sum(j)
            read_gt7_read_num_2 += len(j)
            read_gt7_len_2.append(np.round(np.mean(j)))
            read_gt7_len_list_2.extend(j)                    #计算N50
    read_gt7_len_N50_2 = Compute_N50(read_gt7_len_list_2)
if 8 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 8:
            read_gt8_2 += np.sum(j)
            read_gt8_read_num_2 += len(j)
            read_gt8_len_2.append(np.round(np.mean(j)))
            read_gt8_len_list_2.extend(j)                    #计算N50
    read_gt8_len_N50_2 = Compute_N50(read_gt8_len_list_2)

if 9 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 9:
            read_gt9_2 += np.sum(j)
            read_gt9_read_num_2 += len(j)
            read_gt9_len_2.append(np.round(np.mean(j)))
            read_gt9_len_list_2.extend(j)                    #计算N50
    read_gt9_len_N50_2 = Compute_N50(read_gt9_len_list_2)
if 10 in subread_filter_2_read_dict:
    for i,j in subread_filter_2_read_dict.items():
        if i >= 10:
            read_gt10_2 += np.sum(j)
            read_gt10_read_num_2 += len(j)
            read_gt10_len_2.append(np.round(np.mean(j)))
            read_gt10_len_list_2.extend(j)                    #计算N50
    read_gt10_len_N50_2 = Compute_N50(read_gt10_len_list_2)




print("Subread筛选后subread≥3的Read Num为: " +str(read_gt3_read_num_2))
print("Subread筛选后subread≥4的Read Num为: " +str(read_gt4_read_num_2))
print("Subread筛选后subread≥5的Read Num为: " +str(read_gt5_read_num_2))
print("Subread筛选后subread≥6的Read Num为: " +str(read_gt6_read_num_2))
print("Subread筛选后subread≥7的Read Num为: " +str(read_gt7_read_num_2))
print("Subread筛选后subread≥8的Read Num为: " +str(read_gt8_read_num_2))
print("Subread筛选后subread≥9的Read Num为: " +str(read_gt9_read_num_2))
print("Subread筛选后subread≥10的Read Num为: " +str(read_gt10_read_num_2))

print("---------------------------------------------------")

print("Subread筛选后subread≥3的Base Num为: " +str(read_gt3_2))
print("Subread筛选后subread≥4的Base Num为: " +str(read_gt4_2))
print("Subread筛选后subread≥5的Base Num为: " +str(read_gt5_2))
print("Subread筛选后subread≥6的Base Num为: " +str(read_gt6_2))
print("Subread筛选后subread≥7的Base Num为: " +str(read_gt7_2))
print("Subread筛选后subread≥8的Base Num为: " +str(read_gt8_2))
print("Subread筛选后subread≥9的Base Num为: " +str(read_gt9_2))
print("Subread筛选后subread≥10的Base Num为: " +str(read_gt10_2))

print("---------------------------------------------------")

print("Subread筛选后subread≥3的Read AVG Length为: " +str(np.round(np.mean(read_gt3_len_2))))
print("Subread筛选后subread≥4的Read AVG Length为: " +str(np.round(np.mean(read_gt4_len_2))))
print("Subread筛选后subread≥5的Read AVG Length为: " +str(np.round(np.mean(read_gt5_len_2))))
print("Subread筛选后subread≥6的Read AVG Length为: " +str(np.round(np.mean(read_gt6_len_2))))
print("Subread筛选后subread≥7的Read AVG Length为: " +str(np.round(np.mean(read_gt7_len_2))))
print("Subread筛选后subread≥8的Read AVG Length为: " +str(np.round(np.mean(read_gt8_len_2))))
print("Subread筛选后subread≥9的Read AVG Length为: " +str(np.round(np.mean(read_gt9_len_2))))
print("Subread筛选后subread≥10的Read AVG Length为: " +str(np.round(np.mean(read_gt10_len_2))))

print("---------------------------------------------------")

print("Subread筛选后subread≥3的Read N50为: " +str(read_gt3_len_N50_2))
print("Subread筛选后subread≥4的Read N50为: " +str(read_gt4_len_N50_2))
if 5 in subread_filter_2_read_dict:
    print("Subread筛选后subread≥5的Read N50为: " +str(read_gt5_len_N50_2))
else:
    print("Subread筛选后subread≥5的Read N50为: " +str(0))

if 6 in subread_filter_2_read_dict:
    print("Subread筛选后subread≥6的Read N50为: " +str(read_gt6_len_N50_2))
else:
    print("Subread筛选后subread≥6的Read N50为: " +str(0))
    
if 7 in subread_filter_2_read_dict:
    print("Subread筛选后subread≥7的Read N50为: " +str(read_gt7_len_N50_2))
else:
    print("Subread筛选后subread≥7的Read N50为: " +str(0))

if 8 in subread_filter_2_read_dict:
    print("Subread筛选后subread≥8的Read N50为: " +str(read_gt8_len_N50_2))
else:
    print("Subread筛选后subread≥8的Read N50为: " +str(0))
    
if 9 in subread_filter_2_read_dict:
    print("Subread筛选后subread≥9的Read N50为: " +str(read_gt9_len_N50_2))
else:
    print("Subread筛选后subread≥9的Read N50为: " +str(0))
    
if 10 in subread_filter_2_read_dict:
    print("Subread筛选后subread≥10的Read N50为: " +str(read_gt10_len_N50_2))
else:
    print("Subread筛选后subread≥10的Read N50为: " +str(0))

print("---------------------------------------------------")



#print(raw_read_len_dict)






#raw_hairpin_distribution_dic为原始识别到的hairpin数目  出这个数目
total_read_num = np.sum(list(raw_hairpin_distribution_dic.values()))
print("Subread识别后subread≥3的Read Num占比为: " + str(np.round(read_gt3_read_num/total_read_num,3)))
print("Subread识别后subread≥4的Read Num占比为: " + str(np.round(read_gt4_read_num/total_read_num,3)))
print("Subread识别后subread≥5的Read Num占比为: " + str(np.round(read_gt5_read_num/total_read_num,3)))
print("Subread识别后subread≥6的Read Num占比为: " + str(np.round(read_gt6_read_num/total_read_num,3)))
print("Subread识别后subread≥7的Read Num占比为: " + str(np.round(read_gt7_read_num/total_read_num,3)))
print("Subread识别后subread≥8的Read Num占比为: " + str(np.round(read_gt8_read_num/total_read_num,3)))
print("Subread识别后subread≥9的Read Num占比为: " + str(np.round(read_gt9_read_num/total_read_num,3)))
print("Subread识别后subread≥10的Read Num占比为: " + str(np.round(read_gt10_read_num/total_read_num,3)))

print("---------------------------------------------------")
print("Subread筛选后subread≥3的Read Num占比为: " + str(np.round(read_gt3_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥4的Read Num占比为: " + str(np.round(read_gt4_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥5的Read Num占比为: " + str(np.round(read_gt5_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥6的Read Num占比为: " + str(np.round(read_gt6_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥7的Read Num占比为: " + str(np.round(read_gt7_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥8的Read Num占比为: " + str(np.round(read_gt8_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥9的Read Num占比为: " + str(np.round(read_gt9_read_num_2/total_read_num,3)))
print("Subread筛选后subread≥10的Read Num占比为: " + str(np.round(read_gt10_read_num_2/total_read_num,3)))

#first_filt_dic为第一次过滤得到的hairpin数目
#second_filt_dic为最终的数目  出这个数目



print(raw_hairpin_distribution_dic)
print(first_filt_dic)
print(second_filt_dic)

end_time = time.time()
# 运行时间run_time。round()函数取整
run_time = round(end_time-begin_time)
# 计算时分秒
hour = run_time//3600
minute = (run_time-3600*hour)//60
second = run_time-3600*hour-60*minute
# 输出
print (f'该程序运行时间：{hour}小时{minute}分钟{second}秒')