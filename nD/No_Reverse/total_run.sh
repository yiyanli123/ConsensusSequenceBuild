#!/bin/bash

#echo "终端输入脚本,参数数量,参数: $0,$#,$@ 。"
#echo "终端输入参数,数量: $@,$#,"
help() {

    echo "可选参数可以下列参数组合:"
    # echo "注意 -r -u -f -c 参数请勿同时使用"
    echo "    --h          : 帮助信息"
    echo "    --fq_file fastq file pathway.      : 原始fastq文件地址"
    echo "    --ref_file fastq file pathway.      : Reference文件地址"
    echo "    --res_dir result directory.   : 输出目录"
    echo "    --hairpin_seq  hairpin1 sequence.  : 第一段hairpin序列"
    echo "    --err_threshold  erro sequence threshold. : 可接受错误率阈值"
    echo "    --step_len  step length.  : 扫描步长"
    echo "    --adjusted_len  adjusted length.  :hairpin矫正长度."
    echo "    --subread_len_filter_thre  subread length filter threshold.  :hairpin1和hairpin2之间的poly长度"

    exit 1

}
 
if [[ $# == 0 || "$1" == "--h" ]]; then
    help
    exit 1
fi

# Define default values
subread_len_filter_thre=0.8
adjusted_len=200
err_threshold=0.3
step_len=1

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fq_file) fq_file="$2"; shift ;;
    --ref_file) ref_file="$2"; shift ;;
    --res_dir) res_dir="$2"; shift;;
    --hairpin_seq) hairpin_seq="$2"; shift;;
    --err_threshold) err_threshold="$2"; shift;;
    --step_len) step_len="$2"; shift;;
    --adjusted_len) adjusted_len="$2"; shift;;
    --subread_len_filter_thre) subread_len_filter_thre="$2"; shift;;
    *)  echo "Unknown option: $1" ; exit 1 ;;

  esac
  shift
done

echo "fastq文件为:$fq_file"
echo "Reference文件为: $ref_file"
echo "res_dir为:$res_dir"
echo "hairpin_seq为:$hairpin_seq"
echo "err_threshold为:$err_threshold"
echo "step_len为:$step_len"
echo "adjusted_len为:$adjusted_len"
echo "subread_len_filter_thre为:$subread_len_filter_thre"

startTime=`date +"%Y-%m-%d %H:%M:%S"`

#核心拆分程序
python main.py --fq_file $fq_file --res_dir $res_dir --hairpin_seq $hairpin_seq 
#1>log 2>err.log 

#consensus操作
#进入目录,对各类subread进行consensus操作
cd $res_dir

#: << EOF
#EOF


for i in {2,3,4,5,6,7,8,9}; do 
  for j in {3,4,5,6,7,8,9,10}; do 
    if [ -d "hairpin${i}_detail/subread_use${j}" ]; then 
      
      if [ $j -eq 3 ]; then 
        cd hairpin${i}_detail/subread_use${j}
        #echo "hairpin${i}_detail/subread_use${j} 存在"
        #echo `pwd`
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons >> ../../Consensus_gt_${j}.fa && rm $f.cons ; done
        cd ../../
        #echo `pwd`
      elif [ $j -eq 4 ]; then 
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa ../../Consensus_gt_3.fa >/dev/null  && rm $f.cons ; done
        cd ../../
      elif [ $j -eq 5 ]; then 
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa  ../../Consensus_gt_3.fa ../../Consensus_gt_4.fa >/dev/null && rm $f.cons ; done
        cd ../../
      elif [ $j -eq 6 ]; then 
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa  ../../Consensus_gt_3.fa ../../Consensus_gt_4.fa ../../Consensus_gt_5.fa >/dev/null && rm $f.cons ; done
        cd ../../
      elif [ $j -eq 7 ]; then
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa  ../../Consensus_gt_3.fa ../../Consensus_gt_4.fa ../../Consensus_gt_5.fa ../../Consensus_gt_6.fa >/dev/null && rm $f.cons ; done
        cd ../../
      elif [ $j -eq 8 ]; then
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa  ../../Consensus_gt_3.fa ../../Consensus_gt_4.fa ../../Consensus_gt_5.fa ../../Consensus_gt_6.fa ../../Consensus_gt_7.fa >/dev/null && rm $f.cons ; done
        cd ../../
      elif [ $j -eq 9 ]; then
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa  ../../Consensus_gt_3.fa ../../Consensus_gt_4.fa ../../Consensus_gt_5.fa ../../Consensus_gt_6.fa ../../Consensus_gt_7.fa ../../Consensus_gt_8.fa >/dev/null && rm $f.cons ; done
        cd ../../
      elif [ $j -eq 10 ]; then
        cd hairpin${i}_detail/subread_use${j}
        for f in *[^.png]; do abpoa $f -o $f.cons 2>/dev/null  &&  sed "s/Consensus_sequence/$f/" $f.cons | tee -a ../../Consensus_gt_${j}.fa  ../../Consensus_gt_3.fa ../../Consensus_gt_4.fa ../../Consensus_gt_5.fa ../../Consensus_gt_6.fa ../../Consensus_gt_7.fa ../../Consensus_gt_8.fa ../../Consensus_gt_9.fa >/dev/null && rm $f.cons ; done
        cd ../../
      fi
    fi
  done 
done

#以上对consensus部分进行了构建
#对consensus序列进行评估
#构建sam
if [ -n "$ref_file" ]; then
    if [ -d ./evaluation_space ]; then 
      echo "evaluation_space目录已存在"
    else
      mkdir evaluation_space
    fi
fi

for i in {3,4,5,6,7,8,9,10}; do
  if [ -e Consensus_gt_${i}.fa ]; then
    if [ -n "$ref_file" ]; then
      minimap2 -ax map-ont --secondary=no --eqx -t 24 ../$ref_file Consensus_gt_${i}.fa  > evaluation_space/Consensus_total_${i}.cns.sam 2>/dev/null
      perl ../count_pass.pl evaluation_space/Consensus_total_${i}.cns.sam 
      #rm evaluation_space/Consensus_total_${i}.cns.sam 
      echo "≥${i}的Consensus序列Identity为:`less evaluation_space/Consensus_total_${i}.cns.total.stat  | grep total_iden_rate | awk -F  ':' '{print $2}'`" >>ConsIden.txt
    else 
      echo "≥${i}的Consensus序列Identity为:0" >>ConsIden.txt
    fi
  #读取*.total.stat里的关键词
  else 
    echo "≥${i}的Consensus序列Identity为:0" >>ConsIden.txt
  fi
done 


#echo "开始评估各个consensus序列特征"
#echo `pwd`

#评估gtx的一些指标
#cd evaluation_space
for i in {3,4,5,6,7,8,9,10}; do
  if [ -e Consensus_gt_${i}.fa ]; then
    echo "≥${i}的Consensus序列的Read Number为:`seqkit stat Consensus_gt_${i}.fa | awk '{printf $4 "\t"}' | awk '{printf $2 "\t"}' | sed 's/,//g'`" >> ConsIden.txt
  else 
    echo "≥${i}的Consensus序列的Read Number为:0" >> ConsIden.txt
  fi 
done 

for i in {3,4,5,6,7,8,9,10}; do
  if [ -e Consensus_gt_${i}.fa ]; then
    echo "≥${i}的Consensus序列的Base Number为:`seqkit stat Consensus_gt_${i}.fa | awk '{printf $5 "\t"}' | awk '{printf $2 "\t"}' | sed 's/,//g'`" >> ConsIden.txt
  else 
    echo "≥${i}的Consensus序列的Base Number为:0" >> ConsIden.txt
  fi
done 

for i in {3,4,5,6,7,8,9,10}; do
  if [ -e Consensus_gt_${i}.fa ]; then
    echo "≥${i}的Consensus序列AVG Length为:`seqkit stat Consensus_gt_${i}.fa | awk '{printf $7 "\t"}' | awk '{printf $2 "\t"}' | sed 's/,//g'`" >> ConsIden.txt
  else
    echo "≥${i}的Consensus序列AVG Length为:0" >> ConsIden.txt
  fi
done 
#cd ..


#文件在res_dir里的ConsIden.txt
#echo "Consensus序列评估完成"

#开始评估subread的Identity
#对subread进行下采样，减少评估时间，采样1/10
if [ -n "$ref_file" ]; then 
    for i in {1..10};do 
      if [ -e subread${i}.fa ]; then
        seqtk sample -s100 subread${i}.fa  10000 >subread${i}_tenth.fa
        #rm subread${i}.fa
        minimap2 -ax map-ont --secondary=no --eqx -t 24 ../$ref_file subread${i}_tenth.fa  > evaluation_space/subread${i}.sam 2>/dev/null
        perl ../count_pass.pl evaluation_space/subread${i}.sam
        rm evaluation_space/subread${i}.sam
        #读取*.total.stat里的关键词
        echo "subread${i}的Identity为:`less evaluation_space/subread${i}.total.stat  | grep total_iden_rate | awk -F  ':' '{print $2}'`" >>ConsIden.txt
      else
        echo "subread${i}的Identity为:0" >>ConsIden.txt
      fi
    done
else
  echo "不需要评估subread的Identity."
  for i in {1..10}; do
    echo "subread${i}的Identity为:0" >>ConsIden.txt
    done
fi 

if [ -n "$ref_file" ]; then
  #评估总体subread的Identity，用抽样的1w条来做
  cat subread?_tenth.fa > Overall_subread.fa
  minimap2 -ax map-ont --secondary=no --eqx -t 24 ../$ref_file Overall_subread.fa  > evaluation_space/Overall_subread.sam 2>/dev/null  
  perl ../count_pass.pl evaluation_space/Overall_subread.sam
  rm evaluation_space/subread${i}.sam
  echo "总体subread的Identity为:`less evaluation_space/Overall_subread.total.stat  | grep total_iden_rate | awk -F  ':' '{print $2}'`" >>ConsIden.txt
else
  echo "不需要评估总体subread的Identity."
  echo "总体subread的Identity为:0" >>ConsIden.txt
  
fi

cd ../
#输出excel表格以及可视化图
python settle_info_out_report.py  --input_dir $res_dir 





endTime=`date +"%Y-%m-%d %H:%M:%S"`
st=`date -d  "$startTime" +%s`
et=`date -d  "$endTime" +%s`
sumTime=$(($et-$st))
echo "Total time is : $sumTime second."