import numpy as np
import pandas as pd 
import os 
import argparse
import xlwt

parser = argparse.ArgumentParser(description='对高精度结果进行表格输出和可视化。')
parser.add_argument('--input_dir', type=str, help='输入目录')
args = parser.parse_args()


congnizaton1 = {"Read Number":[],"Base Number":[],"N50 Length":[], "AVG Length":[],"Read Ratio":[]}
filterd1 = {"Read Number":[],"Base Number":[],"N50 Length":[], "AVG Length":[],"Read Ratio":[]}
construction1 = {"Read Number":[],"Base Number":[],"AVG Length":[],"Identity Rate":[]}
raw_subread = {"AVG Length":[],"Identity Rate":[]}


##读取s.sh的log文件
sh_log = args.input_dir + ".out"
with open(sh_log, "r", encoding="utf-8") as f1:
    log_lines = f1.readlines()
    for i in log_lines:
        if "识别后" in i and "Read Num为" in i :
            #print(i)
            congnizaton1["Read Number"].append(float(i.split(":")[1].strip()))
        elif "识别后" in i and "Read AVG Length" in i :
            congnizaton1["AVG Length"].append(float(i.split(":")[1].strip()))
        elif "识别后" in i and "Base Num" in i :
            congnizaton1["Base Number"].append(float(i.split(":")[1].strip()))
        elif "识别后" in i and "Read N50" in i :
            congnizaton1["N50 Length"].append(float(i.split(":")[1].strip()))
        elif "识别后" in i and "Read Num占比" in i :
            congnizaton1["Read Ratio"].append(float(i.split(":")[1].strip()))
        elif "subread的平均长度" in i :
            raw_subread["AVG Length"].append(float(i.split(":")[1].strip()))
        elif "筛选后" in i and "Read Num为" in i :
            filterd1["Read Number"].append(float(i.split(":")[1].strip()))
        elif "筛选后" in i and "Base Num" in i :
            filterd1["Base Number"].append(float(i.split(":")[1].strip()))
        elif "筛选后" in i and "Read N50" in i :
            filterd1["N50 Length"].append(float(i.split(":")[1].strip()))
        elif "筛选后" in i and "Read AVG Length" in i :
            filterd1["AVG Length"].append(float(i.split(":")[1].strip()))
        elif "筛选后" in i and "Read Num占比" in i :
            filterd1["Read Ratio"].append(float(i.split(":")[1].strip()))

### 读取consensus的一些信息
with open(args.input_dir + "/ConsIden.txt","r", encoding="utf-8") as f2:
    consensus_lines = f2.readlines()
    for i in consensus_lines:
        if "Consensus" in i and "Identity" in i :
            construction1["Identity Rate"].append(float(i.split(":")[1].strip()))
        elif "Consensus" in i and "Read Number" in i :
            construction1["Read Number"].append(float(i.split(":")[1].strip()))
        elif "Consensus" in i and "Base Number" in i :
            construction1["Base Number"].append(float(i.split(":")[1].strip()))
        elif "Consensus" in i and "AVG Length" in i :
            construction1["AVG Length"].append(float(i.split(":")[1].strip()))
        elif "subread" in i and "Identity" in i:
            raw_subread["Identity Rate"].append(float(i.split(":")[1].strip()))

##标准输出表格信息，利用xlwt模块

#创建新的workbook（其实就是创建新的excel）
workbook = xlwt.Workbook(encoding= 'utf-8')
# 创建新的sheet表
worksheet = workbook.add_sheet("sheet1")
style = xlwt.XFStyle()
al = xlwt.Alignment()
al.horz = 0x02  # 设置水平居中
al.vert = 0x01  # 设置垂直居中
style.alignment = al
    
#worksheet.write(0,1, "subread识别后")
worksheet.write_merge(0, 0, 1, 8, 'subread识别后',style)
worksheet.write(1,1, "≥3D",style)
worksheet.write(1,2, "≥4D",style)
worksheet.write(1,3, "≥5D",style)
worksheet.write(1,4, "≥6D",style)
worksheet.write(1,5, "≥7D",style)
worksheet.write(1,6, "≥8D",style)
worksheet.write(1,7, "≥9D",style)
worksheet.write(1,8, "≥10D",style)

worksheet.write(2,0, "Read Number",style)
worksheet.write(3,0, "Base Number",style)
worksheet.write(4,0, "N50 Length",style)
worksheet.write(5,0, "AVG Length",style)
worksheet.write(6,0, "Read Ratio",style)

#part2
worksheet.write_merge(8, 8, 1, 8, 'subread筛选后',style)
worksheet.write(9,1, "≥3D",style)
worksheet.write(9,2, "≥4D",style)
worksheet.write(9,3, "≥5D",style)
worksheet.write(9,4, "≥6D",style)
worksheet.write(9,5, "≥7D",style)
worksheet.write(9,6, "≥8D",style)
worksheet.write(9,7, "≥9D",style)
worksheet.write(9,8, "≥10D",style)

worksheet.write(10,0, "Read Number",style)
worksheet.write(11,0, "Base Number",style)
worksheet.write(12,0, "N50 Length",style)
worksheet.write(13,0, "AVG Length",style)
worksheet.write(14,0, "Read Ratio",style)

#part3

worksheet.write_merge(16, 16, 1, 8, 'Consensus构建后',style)
worksheet.write(17,1, "≥3D",style)
worksheet.write(17,2, "≥4D",style)
worksheet.write(17,3, "≥5D",style)
worksheet.write(17,4, "≥6D",style)
worksheet.write(17,5, "≥7D",style)
worksheet.write(17,6, "≥8D",style)
worksheet.write(17,7, "≥9D",style)
worksheet.write(17,8, "≥10D",style)

worksheet.write(18,0, "Read Number",style)
worksheet.write(19,0, "Base Number",style)
worksheet.write(20,0, "AVG Length",style)
worksheet.write(21,0, "Identity",style)

#part4
worksheet.write_merge(24, 24, 1, 11, '原始subread统计',style)
worksheet.write(25,1, 1,style)
worksheet.write(25,2, 2,style)
worksheet.write(25,3, 3,style)
worksheet.write(25,4, 4,style)
worksheet.write(25,5, 5,style)
worksheet.write(25,6, 6,style)
worksheet.write(25,7, 7,style)
worksheet.write(25,8, 8,style)
worksheet.write(25,9, 9,style)
worksheet.write(25,10, 10,style)
worksheet.write(25,11, "Overall",style)


worksheet.write(26,0, "AVG Length",style)
worksheet.write(27,0, "Identity",style)

#workbook.save("test.xls")

#输出subread识别后内容
row_id = 2  #起始行从2开始
for i,j in congnizaton1.items():
    col_id = 1  #起始列从1开始
    for k in j :
        #print(row_id,col_id)
        worksheet.write(row_id, col_id, k,style)
        col_id += 1 
    row_id += 1

#输出subread筛选后内容
row_id = 10  #起始行从2开始
for i,j in filterd1.items():
    col_id = 1  #起始列从1开始
    for k in j :
        #print(row_id,col_id)
        worksheet.write(row_id, col_id, k,style)
        col_id += 1 
    row_id += 1
    
#输出consensus构建后内容
row_id = 18  #起始行从2开始
for i,j in construction1.items():
    col_id = 1  #起始列从1开始
    for k in j :
        #print(row_id,col_id)
        worksheet.write(row_id, col_id, k,style)
        col_id += 1 
    row_id += 1
    
#输出subread相关内容
row_id = 26  #起始行从2开始
for i,j in raw_subread.items():
    col_id = 1  #起始列从1开始
    for k in j :
        #print(row_id,col_id)
        worksheet.write(row_id, col_id, k,style)
        col_id += 1 
    row_id += 1
    
workbook.save(args.input_dir + "/Summary.xls")

### subrea识别后read数目占比
#congnizaton1_read_type = ["≥3","≥4","≥5","≥6","≥7","≥8"]
congnizaton1_read_type = ["≥3","≥4","≥5","≥6","≥7","≥8","≥9","≥10"]
congnizaton1_read_ratio = congnizaton1["Read Ratio"]
filterd_read_ratio = filterd1["Read Ratio"]

###重新画图 202310123
#两次过滤后的最终结果  base mean identity
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_context({"figure.figsize":(7,7)})

# 创建数据

#bar_data = [2813,4200,5265,4881,3320]

# 创建图形和轴对象

fig, ax1 = plt.subplots()
ax1.set_ylim(0,1)

# 绘制柱状图
g = sns.barplot(x=congnizaton1_read_type, y=congnizaton1_read_ratio,ax=ax1,width=0.8,palette=["gainsboro","lightgray","lightgray","silver","darkgray","gray","dimgray","black"])
ax1.tick_params(labelsize=20)
g.set_xlabel('Number of Useful Subreads',fontsize=20)
g.set_ylabel("Ratio",fontsize=20)
# 显示图形
#plt.show()
plt.savefig(args.input_dir + "/subread_check_read_ratio.png")
plt.close()



sns.set_context({"figure.figsize":(7,7)})
# 创建数据
#congnizaton1_read_type = ["≥3","≥4","≥5","≥6","≥7","≥8"]
#bar_data = [2813,4200,5265,4881,3320]

# 创建图形和轴对象

fig, ax1 = plt.subplots()
ax1.set_ylim(0,1)

# 绘制柱状图
g = sns.barplot(x=congnizaton1_read_type, y=filterd_read_ratio,ax=ax1,width=0.8,palette=["gainsboro","lightgray","lightgray","silver","darkgray","gray","dimgray","black"])
ax1.tick_params(labelsize=20)
g.set_xlabel('Number of Useful Subreads',fontsize=20)
g.set_ylabel("Ratio",fontsize=20)
# 显示图形
#plt.show()
plt.savefig(args.input_dir + "/filtered_read_ratio.png")
plt.close()

subread_id = ["1","2","3","4","5","6","7","8","9","10","Overall"]
subread_length = raw_subread["AVG Length"]
subread_length_df = pd.DataFrame({"subread_id":subread_id,"subread_length":subread_length})


sns.set_context({"figure.figsize":(10,10)})
g = sns.lineplot(x="subread_id",y="subread_length", markers="o",data=subread_length_df)

g.tick_params(labelsize=20)
g.set_xlabel('Subread ID',fontsize=20)
g.set_ylabel("Subread Length",fontsize=20)

plt.savefig(args.input_dir + "/Subread_Length.png")
plt.close()

subread_id = ["1","2","3","4","5","6","7","8","9","10","Overall"]
subread_identity = raw_subread["Identity Rate"]
subread_identity_df = pd.DataFrame({"subread_id":subread_id,"subread_identity":subread_identity})

sns.set_context({"figure.figsize":(10,10)})
g = sns.lineplot(x="subread_id",y="subread_identity", markers="o",data=subread_identity_df)

g.tick_params(labelsize=20)
g.set_xlabel('Subread ID',fontsize=20)
g.set_ylabel("Subread Identity",fontsize=20)

plt.savefig(args.input_dir + "/Subread_Identity.png")
plt.close()

consensus_id = ["≥3","≥4","≥5","≥6","≥7","≥8","≥9","≥10"]
consensus_identity = construction1["Identity Rate"]
consensus_identity_df = pd.DataFrame({"consensus_id":consensus_id,"consensus_identity":consensus_identity})

sns.set_context({"figure.figsize":(10,10)})
g = sns.lineplot(x="consensus_id",y="consensus_identity", markers="o",data=consensus_identity_df)

g.tick_params(labelsize=20)
g.set_xlabel('Consensus Type',fontsize=20)
g.set_ylabel("Consensus Identity",fontsize=20)
g.set_ylim(0.95,1)
plt.savefig(args.input_dir + "/Consensus_Identity.png")
plt.close()