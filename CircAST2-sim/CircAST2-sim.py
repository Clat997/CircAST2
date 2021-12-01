#!usr/bin/python
# -*- coding: UTF-8 -*-
import optparse
import random
import numpy
import os
import shutil
import time
import re
import math


# inputGtfName="E:/乱七八糟/mouse/Temp/mm10.ncbiRefSeq.gtf"
# inputChromosomePath =  'E:/乱七八糟/mouse/Chromosomes/'
# inputCircRNAList="E:/乱七八糟/mouse/simulator_Data/ciri_chr1.txt"

inputGtfName='/home/xiehongwen/my_code_on_linux/temp/mm10.ncbiRefSeq.gtf'
inputChromosomePath = '/home/Pub/Index/Chromosomes_Mus'
inputCircRNAList='/home/xiehongwen/my_code_on_linux/restart_20200708/simulating/simulating_chr_all/overlap1/update_simulate_liner2_update_circ_generation/selectCiri.txt'
inputInsertLength=100
inputReadLength=100
inputTotalNumber=20000000
overlap_rate = 0.8##overlap的比例

reference = open(inputCircRNAList, 'r')
os.mkdir("simulate_temp")
tempfile_path = ("simulate_temp/")
circ_infor = open("circRNA_information_1.txt", "a")
gtf_circ = open(inputGtfName, 'r')  # split gtf files by chromosomes
chr_name_gtf = []
for line in gtf_circ:
    line = line.strip('\n')
    array = line.split()
    if len(array[0]) < 6:
        if array[2] == "exon":
            chr_name = array[0]
            out_gtf = open(tempfile_path + chr_name + ".gtf", "a")
            out_gtf.write(line + "\n")
            out_gtf.close()
            if chr_name not in chr_name_gtf:
                chr_name_gtf.append(chr_name)
gtf_circ.close()

circrna_name = []
gene_name_list = []
chr_name_record = []
chain_record = []
gene_name_set = []
ciri_input = open("ciri_input.txt", "a")
out_liner_gtf = open(tempfile_path + "simulating_liner.gtf", "a")
ciri_dict_name = []
i = 0
circ_bound = {}
min_circ_left = 0
max_circ_right = 0
for infor in reference:
    i = i + 1
    if i == 1:
        infor = infor.strip('\n')
        ciri_input.write(infor + "\n")
    if i > 1 :
        infor = infor.strip('\n')
        circ = infor.split('\t')
        if (circ[8] == 'exon') and (circ[1] in chr_name_gtf) and (float(circ[7]) >= 0.1):##修改只关注chr1染色体
            exon = {}
            exon_num = 1
            exon_num_choose = 0
            gene_gtf = open(tempfile_path + str(circ[1]) + ".gtf", "r")
            gennn = gene_gtf.readlines()
            for line in gennn:
                gene_information = line.strip('\n')
                gene_information = gene_information.split('\t')
                if circ[9][:-1] == gene_information[8].split('"')[1]:
                    if int(circ[2]) <= int(gene_information[3]) and int(circ[3]) >= int(gene_information[4]):  ##在GTF文件中找到找到在环内的exon
                        if exon_num == 1:
                            exon[exon_num] = str(gene_information[3]) + '-' + str(gene_information[4])
                            exon_num += 1
                        elif exon_num > 1:
                            if int(gene_information[3]) > int(exon[int(exon_num) - 1].split('-')[1]):  ##如果有多个exon，第二个exon的起点需大于第一个exon的终点
                                exon[exon_num] = str(gene_information[3]) + '-' + str(gene_information[4])
                                exon_num += 1
            gene_gtf.close()  ##得到了在环上的exon信息：exon
            exon_liner1_outside = []
            gene_name_list.append(circ[9][:-1])
            chr_name_record.append(circ[1])
            chain_record.append(circ[10])
            if len(gene_name_list) > 1:
                if gene_name_list[-1] != gene_name_list[-2]:##只要下一个基因跟上一个基因名不一样了，就结算上一个基因的线性转录本
                    # circ_bound.append(str(min_circ_left) + '-' + str(max_circ_right) + '-' + str(gene_name_list[-2]))
                    exon_record = {}
                    exon_num_record = 1
                    exon_liner1 = []
                    exon_liner2 = []
                    exon_liner_variable1 = []
                    exon_liner_variable2 = []
                    exon_num_liner1 = 0
                    exon_num_liner2 = 0
                    exon_num_choose = 0
                    for line in gennn:
                        gene_information = line.strip('\n')
                        gene_information = gene_information.split('\t')
                        if gene_name_list[-2] == gene_information[8].split('"')[1]:
                            if int(min_circ_left) <= int(gene_information[3]) and int(max_circ_right) >= int(gene_information[4]):  ##在GTF文件中找到找到在环内的exon
                                if exon_num_record == 1:
                                    exon_record[exon_num_record] = str(gene_information[3]) + '-' + str(gene_information[4])
                                    exon_num_record += 1
                                elif exon_num_record > 1:
                                    if int(gene_information[3]) > int(exon_record[int(exon_num_record) - 1].split('-')[1]):  ##如果有多个exon，第二个exon的起点需大于第一个exon的终点
                                        exon_record[exon_num_record] = str(gene_information[3]) + '-' + str(gene_information[4])
                                        exon_num_record += 1
                    gene_gtf.close()
                    
                    
                    if overlap_rate == 0:
                        min_exon_num = 0
                    else:
                        circ_exon_num = len(exon_record)
                        min_rate = 1
                        min_exon_num = 0
                        for i in range(circ_exon_num):
                            percentage = (i + 1)/circ_exon_num
                            if abs(percentage - overlap_rate) < min_rate:
                                min_rate = abs(percentage - overlap_rate)
                                min_exon_num = (i + 1)

                    # exon_count = 0
                    # for line in gennn:
                        # gene_information = line.strip('\n')
                        # gene_information = gene_information.split('\t')
                        # if circ[9][:-1] == gene_information[8].split('"')[1]:
                            # if gene_information[8].split('"')[3] not in transcript_id:
                                # transcript_id.append(gene_information[8].split('"')[3])
                            # if len(transcript_id) <= 1:
                                # if int(gene_information[3]) < int(circ[2]) or int(gene_information[4]) > int(circ[3]):
                                    # exon_liner[exon_num_liner] = str(gene_information[3]) + '-' + str(gene_information[4])
                                    # exon_num_liner += 1
                                # if int(circ[2]) <= int(gene_information[3]) and int(circ[3]) >= int(gene_information[4]) and exon_count < min_exon_num:
                                    # exon_count += 1
                                    # exon_liner[exon_num_liner] = str(gene_information[3]) + '-' + str(gene_information[4])
                                    # exon_num_liner += 1
                    # gene_gtf.close()
                    liner_transcript_name = []
                    exon_count1 = 0
                    exon_count2 = 0
                    for line in gennn:
                        gene_information = line.strip('\n')
                        gene_information = gene_information.split('\t')
                        if gene_name_list[-2] == gene_information[8].split('"')[1]:
                            if(gene_information[8].split('"')[3] not in liner_transcript_name):
                                liner_transcript_name.append(gene_information[8].split('"')[3])
                            if len(liner_transcript_name) == 1:
                                if int(gene_information[4]) < int(min_circ_left) or int(gene_information[3]) > int(max_circ_right):
                                    temp = str(gene_information[3]) + '-' + str(gene_information[4])
                                    if len(exon_liner1) == 0:
                                        exon_liner1.append(temp)
                                        exon_num_liner1 += 1
                                    elif int(gene_information[3]) > int(exon_liner1[-1].split('-')[1]):
                                        exon_liner1.append(temp)
                                        exon_num_liner1 += 1
                                if int(min_circ_left) <= int(gene_information[3]) and int(max_circ_right) >= int(gene_information[4]) and exon_count1 < min_exon_num:
                                    exon_count1 += 1
                                    temp = str(gene_information[3]) + '-' + str(gene_information[4])
                                    if len(exon_liner_variable1) == 0:
                                        exon_liner_variable1.append(temp)
                                    elif int(gene_information[3]) > int(exon_liner_variable1[-1].split('-')[1]):
                                        exon_liner_variable1.append(temp)
                            if len(liner_transcript_name) == 2:
                                if int(gene_information[4]) < int(min_circ_left) or int(gene_information[3]) > int(max_circ_right):
                                    temp = str(gene_information[3]) + '-' + str(gene_information[4])
                                    if len(exon_liner2) == 0:
                                        exon_liner2.append(temp)
                                        exon_num_liner2 += 1
                                    elif int(gene_information[3]) > int(exon_liner2[-1].split('-')[1]):
                                        exon_liner2.append(temp)
                                        exon_num_liner2 += 1
                                if int(min_circ_left) <= int(gene_information[3]) and int(max_circ_right) >= int(gene_information[4]) and exon_count2 < min_exon_num:
                                    exon_count2 += 1
                                    temp = str(gene_information[3]) + '-' + str(gene_information[4])
                                    if len(exon_liner_variable2) == 0:
                                        exon_liner_variable2.append(temp)
                                    elif int(gene_information[3]) > int(exon_liner_variable2[-1].split('-')[1]):
                                        exon_liner_variable2.append(temp)
                    gene_gtf.close()
                    exon_liner1_outside = exon_liner1

                    if exon_num_liner1 != 0:
                        write_liner_gtf = []
                        write_liner_gtf.append(chr_name_record[-2])  ##染色体
                        write_liner_gtf.append("ncbiRefSeq")  ##
                        write_liner_gtf.append("transcript")  ##
                        if min_exon_num != 0 and len(exon_liner_variable1) != 0:
                            if (int(exon_liner_variable1[0].split('-')[0]) < int(exon_liner1[0].split('-')[0])):
                                left_cor = exon_liner_variable1[0].split('-')[0]
                            else:
                                left_cor = exon_liner1[0].split('-')[0]  ##
                            if (int(exon_liner_variable1[-1].split('-')[1]) > int(exon_liner1[-1].split('-')[1])):
                                right_cor = exon_liner_variable1[-1].split('-')[1]  ##
                            else:
                                right_cor = exon_liner1[-1].split('-')[1]  ##
                        else:
                            left_cor = exon_liner1[0].split('-')[0]  ##
                            right_cor = exon_liner1[-1].split('-')[1]  ##
                        write_liner_gtf.append(left_cor)  ##
                        write_liner_gtf.append(right_cor)  ##
                        write_liner_gtf.append(".")  ##
                        write_liner_gtf.append(chain_record[-2])  ##
                        write_liner_gtf.append(".")
                        write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                        imforma_gtf = '\t'.join(write_liner_gtf)
                        out_liner_gtf.write(imforma_gtf)
                        out_liner_gtf.write('\n')
                        circ_information = []
                        exon_length = []
                        exon_start_loci = []
                        exon_num_choose = 0
                        write_overlap_count = 0
                        rist = range(0, exon_num_liner1)  ##全部
                        for ex in rist:
                            k = int(ex)
                            start_end = exon_liner1[k].split('-')
                            if min_exon_num != 0 and len(exon_liner_variable1) != 0:
                                if (int(start_end[0]) > int(exon_liner_variable1[0].split('-')[1]) or ((k + 1) == exon_num_liner1 and (int(start_end[1]) < int(exon_liner_variable1[0].split('-')[0])))) and write_overlap_count == 0:
                                    write_overlap_count += 1
                                    if ((k + 1) == exon_num_liner1 and (int(start_end[1]) < int(exon_liner_variable1[0].split('-')[0]))):##如果是线性外显子右侧位点在右侧可变位点左边，
                                        exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                                        exon_start_loci.append(str(int(start_end[0]) - int(left_cor)))  ##exon相对于线性RNA起始位点的偏移量
                                        exon_num_choose += 1
                                        write_liner_gtf = []
                                        write_liner_gtf.append(chr_name_record[-2])  ##染色体
                                        write_liner_gtf.append("ncbiRefSeq")  ##
                                        write_liner_gtf.append("exon")  ##
                                        write_liner_gtf.append(exon_liner1[k].split('-')[0])  ##
                                        write_liner_gtf.append(exon_liner1[k].split('-')[1])  ##
                                        write_liner_gtf.append(".")  ##
                                        write_liner_gtf.append(chain_record[-2])  ##
                                        write_liner_gtf.append(".")
                                        write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                                        imforma_gtf = '\t'.join(write_liner_gtf)
                                        out_liner_gtf.write(imforma_gtf)
                                        out_liner_gtf.write('\n')
                                        for i in range(len(exon_liner_variable1)):
                                            exon_length.append(str(int(exon_liner_variable1[i].split('-')[1]) + 1 - int(exon_liner_variable1[i].split('-')[0])))
                                            exon_start_loci.append(str(int(exon_liner_variable1[i].split('-')[0]) - int(left_cor)))
                                            exon_num_choose += 1
                                            write_liner_gtf = []
                                            write_liner_gtf.append(chr_name_record[-2])  ##染色体
                                            write_liner_gtf.append("ncbiRefSeq")  ##
                                            write_liner_gtf.append("exon")  ##
                                            write_liner_gtf.append(exon_liner_variable1[i].split('-')[0])  ##
                                            write_liner_gtf.append(exon_liner_variable1[i].split('-')[1])  ##
                                            write_liner_gtf.append(".")  ##
                                            write_liner_gtf.append(chain_record[-2])  ##
                                            write_liner_gtf.append(".")
                                            write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                                            imforma_gtf = '\t'.join(write_liner_gtf)
                                            out_liner_gtf.write(imforma_gtf)
                                            out_liner_gtf.write('\n')
                                        break
                                        
                                    for i in range(len(exon_liner_variable1)):
                                        exon_length.append(str(int(exon_liner_variable1[i].split('-')[1]) + 1 - int(exon_liner_variable1[i].split('-')[0])))
                                        exon_start_loci.append(str(int(exon_liner_variable1[i].split('-')[0]) - int(left_cor)))
                                        exon_num_choose += 1
                                        write_liner_gtf = []
                                        write_liner_gtf.append(chr_name_record[-2])  ##染色体
                                        write_liner_gtf.append("ncbiRefSeq")  ##
                                        write_liner_gtf.append("exon")  ##
                                        write_liner_gtf.append(exon_liner_variable1[i].split('-')[0])  ##
                                        write_liner_gtf.append(exon_liner_variable1[i].split('-')[1])  ##
                                        write_liner_gtf.append(".")  ##
                                        write_liner_gtf.append(chain_record[-2])  ##
                                        write_liner_gtf.append(".")
                                        write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                                        imforma_gtf = '\t'.join(write_liner_gtf)
                                        out_liner_gtf.write(imforma_gtf)
                                        out_liner_gtf.write('\n')
                            exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                            exon_start_loci.append(str(int(start_end[0]) - int(left_cor)))  ##exon相对于线性RNA起始位点的偏移量
                            exon_num_choose += 1
                            write_liner_gtf = []
                            write_liner_gtf.append(chr_name_record[-2])  ##染色体
                            write_liner_gtf.append("ncbiRefSeq")  ##
                            write_liner_gtf.append("exon")  ##
                            write_liner_gtf.append(exon_liner1[k].split('-')[0])  ##
                            write_liner_gtf.append(exon_liner1[k].split('-')[1])  ##
                            write_liner_gtf.append(".")  ##
                            write_liner_gtf.append(chain_record[-2])  ##
                            write_liner_gtf.append(".")
                            write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                            imforma_gtf = '\t'.join(write_liner_gtf)
                            out_liner_gtf.write(imforma_gtf)
                            out_liner_gtf.write('\n')
                        circ_information.append(chr_name_record[-2])  ##染色体？
                        circ_information.append(left_cor)  ##mRNA起始
                        circ_information.append(right_cor)  ##最后一个环内exon的右端位点
                        circ_information.append(gene_name_list[-2])## 基因名
                        circ_information.append(chain_record[-2]) ##正负链
                        circ_information.append(str(exon_num_choose))
                        circ_information.append(','.join(exon_length))
                        circ_information.append(','.join(exon_start_loci))
                        circ_information.append('liner')  ## isform属性
                        informa = '\t'.join(circ_information)
                        circ_infor.write(informa)
                        circ_infor.write('\n')
                    
                    if exon_num_liner2 != 0:
                        write_liner_gtf = []
                        write_liner_gtf.append(chr_name_record[-2])  ##染色体
                        write_liner_gtf.append("ncbiRefSeq")  ##
                        write_liner_gtf.append("transcript")  ##
                        if min_exon_num != 0 and len(exon_liner_variable2) != 0:
                            if (int(exon_liner_variable2[0].split('-')[0]) < int(exon_liner2[0].split('-')[0])):
                                left_cor = exon_liner_variable2[0].split('-')[0]
                            else:
                                left_cor = exon_liner2[0].split('-')[0]  ##
                            if (int(exon_liner_variable2[-1].split('-')[1]) > int(exon_liner2[-1].split('-')[1])):
                                right_cor = exon_liner_variable2[-1].split('-')[1]  ##
                            else:
                                right_cor = exon_liner2[-1].split('-')[1]  ##
                        else:
                            left_cor = exon_liner2[0].split('-')[0]  ##
                            right_cor = exon_liner2[-1].split('-')[1]  ##
                        write_liner_gtf.append(left_cor)  ##
                        write_liner_gtf.append(right_cor)  ##
                        write_liner_gtf.append(".")  ##
                        write_liner_gtf.append(chain_record[-2])  ##
                        write_liner_gtf.append(".")
                        write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                        imforma_gtf = '\t'.join(write_liner_gtf)
                        out_liner_gtf.write(imforma_gtf)
                        out_liner_gtf.write('\n')
                        circ_information = []
                        exon_length = []
                        exon_start_loci = []
                        exon_num_choose = 0
                        write_overlap_count = 0
                        rist = range(0, exon_num_liner2)  ##全部
                        for ex in rist:
                            k = int(ex)
                            start_end = exon_liner2[k].split('-')
                            if min_exon_num != 0 and len(exon_liner_variable2) != 0:
                                if (int(start_end[0]) > int(exon_liner_variable2[0].split('-')[1]) or ((k + 1) == exon_num_liner2 and (int(start_end[1]) < int(exon_liner_variable2[0].split('-')[0])))) and write_overlap_count == 0:
                                    write_overlap_count += 1
                                    if ((k + 1) == exon_num_liner2 and (int(start_end[1]) < int(exon_liner_variable2[0].split('-')[0]))):##如果是线性外显子右侧位点在右侧可变位点左边，
                                        exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                                        exon_start_loci.append(str(int(start_end[0]) - int(left_cor)))  ##exon相对于线性RNA起始位点的偏移量
                                        exon_num_choose += 1
                                        write_liner_gtf = []
                                        write_liner_gtf.append(chr_name_record[-2])  ##染色体
                                        write_liner_gtf.append("ncbiRefSeq")  ##
                                        write_liner_gtf.append("exon")  ##
                                        write_liner_gtf.append(exon_liner2[k].split('-')[0])  ##
                                        write_liner_gtf.append(exon_liner2[k].split('-')[1])  ##
                                        write_liner_gtf.append(".")  ##
                                        write_liner_gtf.append(chain_record[-2])  ##
                                        write_liner_gtf.append(".")
                                        write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                                        imforma_gtf = '\t'.join(write_liner_gtf)
                                        out_liner_gtf.write(imforma_gtf)
                                        out_liner_gtf.write('\n')
                                        for i in range(len(exon_liner_variable2)):
                                            exon_length.append(str(int(exon_liner_variable2[i].split('-')[1]) + 1 - int(exon_liner_variable2[i].split('-')[0])))
                                            exon_start_loci.append(str(int(exon_liner_variable2[i].split('-')[0]) - int(left_cor)))
                                            exon_num_choose += 1
                                            write_liner_gtf = []
                                            write_liner_gtf.append(chr_name_record[-2])  ##染色体
                                            write_liner_gtf.append("ncbiRefSeq")  ##
                                            write_liner_gtf.append("exon")  ##
                                            write_liner_gtf.append(exon_liner_variable2[i].split('-')[0])  ##
                                            write_liner_gtf.append(exon_liner_variable2[i].split('-')[1])  ##
                                            write_liner_gtf.append(".")  ##
                                            write_liner_gtf.append(chain_record[-2])  ##
                                            write_liner_gtf.append(".")
                                            write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                                            imforma_gtf = '\t'.join(write_liner_gtf)
                                            out_liner_gtf.write(imforma_gtf)
                                            out_liner_gtf.write('\n')
                                        break
                                        
                                    for i in range(len(exon_liner_variable2)):
                                        exon_length.append(str(int(exon_liner_variable2[i].split('-')[1]) + 1 - int(exon_liner_variable2[i].split('-')[0])))
                                        exon_start_loci.append(str(int(exon_liner_variable2[i].split('-')[0]) - int(left_cor)))
                                        exon_num_choose += 1
                                        write_liner_gtf = []
                                        write_liner_gtf.append(chr_name_record[-2])  ##染色体
                                        write_liner_gtf.append("ncbiRefSeq")  ##
                                        write_liner_gtf.append("exon")  ##
                                        write_liner_gtf.append(exon_liner_variable2[i].split('-')[0])  ##
                                        write_liner_gtf.append(exon_liner_variable2[i].split('-')[1])  ##
                                        write_liner_gtf.append(".")  ##
                                        write_liner_gtf.append(chain_record[-2])  ##
                                        write_liner_gtf.append(".")
                                        write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                                        imforma_gtf = '\t'.join(write_liner_gtf)
                                        out_liner_gtf.write(imforma_gtf)
                                        out_liner_gtf.write('\n')
                            exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                            exon_start_loci.append(str(int(start_end[0]) - int(left_cor)))  ##exon相对于线性RNA起始位点的偏移量
                            exon_num_choose += 1
                            write_liner_gtf = []
                            write_liner_gtf.append(chr_name_record[-2])  ##染色体
                            write_liner_gtf.append("ncbiRefSeq")  ##
                            write_liner_gtf.append("exon")  ##
                            write_liner_gtf.append(exon_liner2[k].split('-')[0])  ##
                            write_liner_gtf.append(exon_liner2[k].split('-')[1])  ##
                            write_liner_gtf.append(".")  ##
                            write_liner_gtf.append(chain_record[-2])  ##
                            write_liner_gtf.append(".")
                            write_liner_gtf.append("gene_id \"" + str(gene_name_list[-2]) + "\"; transcript_id \"" + str(gene_name_list[-2]) + "_1\"")
                            imforma_gtf = '\t'.join(write_liner_gtf)
                            out_liner_gtf.write(imforma_gtf)
                            out_liner_gtf.write('\n')
                        circ_information.append(chr_name_record[-2])  ##染色体？
                        circ_information.append(left_cor)  ##mRNA起始
                        circ_information.append(right_cor)  ##最后一个环内exon的右端位点
                        circ_information.append(gene_name_list[-2])## 基因名
                        circ_information.append(chain_record[-2]) ##正负链
                        circ_information.append(str(exon_num_choose))
                        circ_information.append(','.join(exon_length))
                        circ_information.append(','.join(exon_start_loci))
                        circ_information.append('liner')  ## isform属性
                        informa = '\t'.join(circ_information)
                        circ_infor.write(informa)
                        circ_infor.write('\n')
                
                
                
            if circ[9][:-1] not in gene_name_set:
                gene_name_set.append(circ[9][:-1])
                min_circ_left = int(circ[2])
                max_circ_right = int(circ[3])
            if circ[9][:-1] in gene_name_set:
                if int(circ[2]) < min_circ_left:
                    min_circ_left = int(circ[2])
                if int(circ[3]) > max_circ_right:
                    max_circ_right = int(circ[3])
            
            circ_information = []
            exon_length = []
            exon_start_loci = []
            exon_num_choose = 0
            if len(exon) > 3:
                iso_num = random.randint(1, 100)
                if 1 <= iso_num <= 50:
                    # rist = range(1, exon_num)  ##全部
                    # # for i in range(1,len(exon) + 1):
                        # # if exon[i] not in exon_liner1:
                            # # rist_list.remove(i)
                    # for ex in rist:
                        # k = int(ex)
                        # start_end = exon[k].split('-')
                        # exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        # exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))  ##exon相对于环形RNA起始位点的偏移量
                        # exon_num_choose += 1
                    # if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        # circ_information.append(circ[1])  ##染色体？
                        # circ_information.append(str(int(circ[2])))  ##circRNA起始
                        # circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        # circ_information.append(circ[9][:-1])## 基因名
                        # circ_information.append(circ[10]) ##正负链
                        # circ_information.append(str(exon_num_choose))
                        # circ_information.append(','.join(exon_length))
                        # circ_information.append(','.join(exon_start_loci))
                        # circ_information.append(str(float(circ[7])))
                        # circ_information.append('circ')  ## isform属性
                        # informa = '\t'.join(circ_information)
                        # circ_infor.write(informa)
                        # circ_infor.write('\n')
                        # if circ[0] not in circrna_name:
                            # circrna_name.append(circ[0])
                            # ciri_input.write(infor)
                            # ciri_input.write("\n")
                            # if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                # ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))
                    exon_num_choose = 0
                    circ_information = []
                    exon_length = []
                    exon_start_loci = []
                    rist_2 = range(1, exon_num)
                    rist_2_list = list(rist_2)
                    rist_2_list.remove(random.randint(2, int(exon_num) - 2))  # 中间删一个
                    for ex in rist_2_list:
                        k = int(ex)
                        start_end = exon[k].split('-')
                        exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                        exon_num_choose += 1
                    if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        circ_information.append(circ[1])  ##染色体？
                        circ_information.append(str(int(circ[2])))  ##circRNA起始
                        circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        circ_information.append(circ[9][:-1])
                        circ_information.append(circ[10])
                        circ_information.append(str(exon_num_choose))
                        circ_information.append(','.join(exon_length))
                        circ_information.append(','.join(exon_start_loci))
                        circ_information.append(str(float(circ[7])))
                        circ_information.append('circ')  ## isform属性
                        informa = '\t'.join(circ_information)
                        circ_infor.write(informa)
                        circ_infor.write('\n')
                        if circ[0] not in circrna_name:
                            circrna_name.append(circ[0])
                            ciri_input.write(infor)
                            ciri_input.write("\n")
                            if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))
                    # exon_num_choose = 0
                    # circ_information = []
                    # exon_length = []
                    # exon_start_loci = []
                    # rist_3 = range(1, exon_num)
                    # rrlist = rist_3[1:int(exon_num - 2)]
                    # dele_exon = random.sample(rrlist, 2)
                    # rist_3_list = list(rist_3)
                    # rist_3_list.remove(dele_exon[0])
                    # rist_3_list.remove(dele_exon[1])  ##左边删两个
                    # for i in range(1,len(exon) + 1):
                        # if exon[i] not in exon_liner1_outside and i in rist_3_list and i != 1 and i!= len(exon):
                            # rist_3_list.remove(i)
                    # for ex in rist_3_list:
                        # k = int(ex)
                        # start_end = exon[k].split('-')
                        # exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        # exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                        # exon_num_choose += 1
                    # if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        # circ_information.append(circ[1])  ##染色体？
                        # circ_information.append(str(int(circ[2])))  ##circRNA起始
                        # circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        # circ_information.append(circ[9][:-1])
                        # circ_information.append(circ[10])
                        # circ_information.append(str(exon_num_choose))
                        # circ_information.append(','.join(exon_length))
                        # circ_information.append(','.join(exon_start_loci))
                        # circ_information.append(str(float(circ[7])))
                        # circ_information.append('circ')  ## isform属性
                        # informa = '\t'.join(circ_information)
                        # circ_infor.write(informa)
                        # circ_infor.write('\n')
                        # if circ[0] not in circrna_name:
                            # circrna_name.append(circ[0])
                            # ciri_input.write(infor)
                            # ciri_input.write("\n")
                            # if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                # ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))

                elif 50 < iso_num <= 100:
                    rist = range(1, exon_num)
                    rist_list = list(rist)
                    rist_list.remove(random.randint(2, int(exon_num) - 2))  ##中间删
                    for ex in rist_list:
                        k = int(ex)
                        start_end = exon[k].split('-')
                        exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                        exon_num_choose += 1
                    if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        circ_information.append(circ[1])  ##染色体？
                        circ_information.append(str(int(circ[2])))  ##circRNA起始
                        circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        circ_information.append(circ[9][:-1])
                        circ_information.append(circ[10])
                        circ_information.append(str(exon_num_choose))
                        circ_information.append(','.join(exon_length))
                        circ_information.append(','.join(exon_start_loci))
                        circ_information.append(str(float(circ[7])))
                        circ_information.append('circ')  ## isform属性
                        informa = '\t'.join(circ_information)
                        circ_infor.write(informa)
                        circ_infor.write('\n')
                        if circ[0] not in circrna_name:
                            circrna_name.append(circ[0])
                            ciri_input.write(infor)
                            ciri_input.write("\n")
                            if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))
                    exon_num_choose = 0
                    circ_information = []
                    exon_length = []
                    exon_start_loci = []
                    rist_2 = range(1, exon_num)
                    rrlist = rist_2[1:int(exon_num - 2)]
                    dele_exon = random.sample(rrlist, 2)
                    rist_2_list = list(rist_2)
                    rist_2_list.remove(dele_exon[0])
                    rist_2_list.remove(dele_exon[1])  ##左边删
                    for ex in rist_2_list:
                        k = int(ex)
                        start_end = exon[k].split('-')
                        exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                        exon_num_choose += 1
                    if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        circ_information.append(circ[1])  ##染色体？
                        circ_information.append(str(int(circ[2])))  ##circRNA起始
                        circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        circ_information.append(circ[9][:-1])
                        circ_information.append(circ[10])
                        circ_information.append(str(exon_num_choose))
                        circ_information.append(','.join(exon_length))
                        circ_information.append(','.join(exon_start_loci))
                        circ_information.append(str(float(circ[7])))
                        circ_information.append('circ')  ## isform属性
                        informa = '\t'.join(circ_information)
                        circ_infor.write(informa)
                        circ_infor.write('\n')
                        if circ[0] not in circrna_name:
                            circrna_name.append(circ[0])
                            ciri_input.write(infor)
                            ciri_input.write("\n")
                            if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))


                # elif 50 < iso_num <= 100:
                    # rist = range(1, exon_num)
                    # rist_list = list(rist)
                    # for i in range(1,len(exon) + 1):
                        # if exon[i] not in exon_liner1 and i in rist_list and i != 1 and i!= len(exon):
                            # rist_list.remove(i)
                    # for ex in rist_list:
                        # k = int(ex)
                        # start_end = exon[k].split('-')
                        # exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        # exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                        # exon_num_choose += 1
                    # if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        # circ_information.append(circ[1])  ##染色体？
                        # circ_information.append(str(int(circ[2])))  ##circRNA起始
                        # circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        # circ_information.append(circ[9][:-1])
                        # circ_information.append(circ[10])
                        # circ_information.append(str(exon_num_choose))
                        # circ_information.append(','.join(exon_length))
                        # circ_information.append(','.join(exon_start_loci))
                        # circ_information.append(str(float(circ[7])))
                        # circ_information.append('circ')  ## isform属性
                        # informa = '\t'.join(circ_information)
                        # circ_infor.write(informa)
                        # circ_infor.write('\n')
                        # if circ[0] not in circrna_name:
                            # circrna_name.append(circ[0])
                            # ciri_input.write(infor)
                            # ciri_input.write("\n")
                            # if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                # ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))
                    # exon_num_choose = 0
                    # circ_information = []
                    # exon_length = []
                    # exon_start_loci = []
                    # rist_2 = range(1, exon_num)
                    # rist_2_list = list(rist_2)
                    # rist_2_list.remove(random.randint(2, int(exon_num) - 2))  # 中间删一个
                    # for i in range(1,len(exon) + 1):
                        # if exon[i] not in exon_liner1 and i in rist_2_list and i != 1 and i!= len(exon):
                            # rist_2_list.remove(i)
                    # for ex in rist_2_list:
                        # k = int(ex)
                        # start_end = exon[k].split('-')
                        # exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                        # exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                        # exon_num_choose += 1
                    # if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                        # circ_information.append(circ[1])  ##染色体？
                        # circ_information.append(str(int(circ[2])))  ##circRNA起始
                        # circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                        # circ_information.append(circ[9][:-1])
                        # circ_information.append(circ[10])
                        # circ_information.append(str(exon_num_choose))
                        # circ_information.append(','.join(exon_length))
                        # circ_information.append(','.join(exon_start_loci))
                        # circ_information.append(str(float(circ[7])))
                        # circ_information.append('circ')  ## isform属性
                        # informa = '\t'.join(circ_information)
                        # circ_infor.write(informa)
                        # circ_infor.write('\n')
                        # if circ[0] not in circrna_name:
                            # circrna_name.append(circ[0])
                            # ciri_input.write(infor)
                            # ciri_input.write("\n")
                            # if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                # ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))

            
            elif len(exon) == 3:
                rist = range(1, exon_num)
                rist_list = list(rist)
                rist_list.remove(2)  # 中间删一个
                for ex in rist_list:
                    k = int(ex)
                    start_end = exon[k].split('-')
                    exon_length.append(str(int(start_end[1]) + 1 - int(start_end[0])))
                    exon_start_loci.append(str(int(start_end[0]) - int(circ[2])))
                    exon_num_choose += 1
                if int(exon_start_loci[0]) == 0 and (int(exon[int(exon_num) - 1].split('-')[1]) == int(circ[3])):
                    circ_information.append(circ[1])  ##染色体？
                    circ_information.append(str(int(circ[2])))  ##circRNA起始
                    circ_information.append(exon[int(exon_num) - 1].split('-')[1])  ##最后一个环内exon的右端位点
                    circ_information.append(circ[9][:-1])
                    circ_information.append(circ[10])
                    circ_information.append(str(exon_num_choose))
                    circ_information.append(','.join(exon_length))
                    circ_information.append(','.join(exon_start_loci))
                    circ_information.append(str(float(circ[7])))
                    circ_information.append('circ')  ## isform属性
                    informa = '\t'.join(circ_information)
                    circ_infor.write(informa)
                    circ_infor.write('\n')
                    if circ[0] not in circrna_name:
                        circrna_name.append(circ[0])
                        ciri_input.write(infor)
                        ciri_input.write("\n")
                        if str(circ[2]) + '-' + str(circ[3]) not in ciri_dict_name:
                                ciri_dict_name.append(str(circ[2]) + '-' + str(circ[3]))
circ_dict_num = numpy.zeros(len(ciri_dict_name))
circ_dict = dict(zip(ciri_dict_name,circ_dict_num))
reference.close()
circ_infor.close()
ciri_input.close()
out_liner_gtf.close()
os.system(' cat circRNA_information_1.txt | sort -u > circRNA_information_before.txt ')
circ_resu = open("circRNA_information.txt", "a")
double_check = []
with open("circRNA_information_before.txt", "r") as circ_before:  ##筛选同一gene的isoform
    for line in circ_before:
        start_gene_name = []
        end_gene_name = []
        jj = line.strip('\n')
        cc_list = jj.split()
        if cc_list[-1] == 'liner':
            if cc_list[3] not in double_check:
                double_check.append(cc_list[3])
                expression_liner = random.uniform(50, 100)   ###合理吗
                cc_list.append(str(expression_liner))  ##表达量？
                resulll = '\t'.join(cc_list)
                circ_resu.write(resulll)
                circ_resu.write('\n')
                expression_liner_big = expression_liner
                continue
            if cc_list[3] in double_check:
                expression_liner = random.uniform(5, 10)   ###合理吗
                cc_list.append(str(expression_liner))  ##表达量？
                resulll = '\t'.join(cc_list)
                circ_resu.write(resulll)
                circ_resu.write('\n')
        if cc_list[-1] == 'circ':
            with open(tempfile_path + str(cc_list[0]) + '.gtf', 'r') as ch:
                for lines_1 in ch:
                    lines_2 = lines_1.strip('\n')
                    lines_2 = lines_2.split('\t')
                    trans_id = lines_2[8].split('"')
                    if lines_2[3] == cc_list[1]:
                        if trans_id[3] not in start_gene_name:
                            start_gene_name.append(trans_id[3])
                    if lines_2[4] == cc_list[2]:
                        if trans_id[3] not in end_gene_name:
                            end_gene_name.append(trans_id[3])
                gene_name_same = list(set(start_gene_name).intersection(set(end_gene_name)))  ##gene_name取交集
                if gene_name_same != []:
                    if float(cc_list[-2]) >= 0.1:
                        expression_circ = expression_liner_big * float(cc_list[-2])
                        # if expression_circ <= 1:
                        #     expression_circ = 1
                        # # if expression_circ <= 0.1:
                        # #     expression_left = 0
                        # #     expression_right = expression_circ + 0.1
                        # # else:
                        # #     expression_left = expression_circ - 0.1
                        # #     expression_right = expression_circ + 0.1
                        expression = random.uniform(expression_circ - 1, expression_circ + 1)
                        cc_list.append(str(expression))  ##表达量？
                        resulll = '\t'.join(cc_list)
                        circ_resu.write(resulll)
                        circ_resu.write('\n')
circ_resu.close()

chromosome={}

chromosome_path_list=os.listdir(str(inputChromosomePath))##fa文件包含的文件名

for filename in chromosome_path_list:
	chromosome_path = os.path.join(str(inputChromosomePath),filename)
	file = open(chromosome_path,'r')
	chr_seq = file.read()
	chr_seq = chr_seq.replace("\n","")
	rl = r'chr(\w*).fa'
	chr_num = re.findall(rl,filename)
	chromosome[chr_num[0]]=chr_seq
	file.close()
start_simulate_time = time.time()

reads_1 = open("reads_1.fastq", "a")
reads_2 = open("reads_2.fastq", "a")
circ_list = open("circRNA_information_before.txt", "a")
circ_isoform_infor = open("circRNA_isoform_list.txt", "a")
junction_record = open("junction_record.txt", "a")
read_num = 0
length = []


with open("circRNA_information.txt", "r") as reference_circ:
    for infor_2 in reference_circ:
        infor_2 = infor_2.strip('\n')
        circ_2 = infor_2.split()
        if circ_2[-2] == 'liner':
            liner_information = []
            liner_exon_len = circ_2[6].split(',')
            liner_exon_start = circ_2[7].split(',')
            exon_num_1 = circ_2[5]
            junction_num = 0
            liner_information_2 = []
            liner_seq_1 = ""
            liner_chr_num = circ_2[0].split("chr")
            liner_alseq = chromosome[liner_chr_num[1]]
            for exon_2 in range(1, int(exon_num_1) + 1):
                k = int(exon_2)
                liner_seq_1 = liner_seq_1 + liner_alseq[(int(circ_2[1]) + int(liner_exon_start[k - 1]) + len(circ_2[0])):(int(circ_2[1]) + int(liner_exon_start[k - 1]) + int(liner_exon_len[k - 1]) + len(circ_2[0]))]  ##  为什么要加len(circ_2[0])？ 这是是把所有exon上的序列合在一起
            liner_seq = liner_seq_1
            seq_length = len(liner_seq)
            liner_time = int(float(circ_2[-1])) * (seq_length) * int(inputTotalNumber) / 1000000000  ##10亿，inputTotalNumber = 20000000。表达量 * isoform长度 * 测序深度 = 取样数目
            count = 1
            while count <= liner_time:
                if seq_length >= inputReadLength:
                    start_loci = random.randint(0, seq_length - inputReadLength - 1)
                else:
                    break
                if length == []:
                    length_array = 60 * numpy.random.randn(1, int((int(inputTotalNumber) / 10))) + (int(inputInsertLength) + (2 * int(inputReadLength)))  ##？
                    length = list(length_array[0])
                    length_array = []
                reads1_seq = liner_seq[start_loci:(start_loci + int(inputReadLength))]
                while int(length[0]) < int(inputReadLength):
                    length.pop(0)
                start_loci_2 = int(start_loci + int(length[0]) - int(inputReadLength))
                if (start_loci_2 + int(inputReadLength)) > seq_length:
                    reads2_seq = liner_seq[(int(seq_length)-int(inputReadLength)):int(seq_length)]
                else:
                    reads2_seq = liner_seq[start_loci_2:(start_loci_2 + int(inputReadLength))]
                reads2_seq = reads2_seq[::-1]  ##取反
                reads2_seq = reads2_seq.replace("A", "X")
                reads2_seq = reads2_seq.replace("a", "X")
                reads2_seq = reads2_seq.replace("T", "A")
                reads2_seq = reads2_seq.replace("t", "A")
                reads2_seq = reads2_seq.replace("X", "T")
                reads2_seq = reads2_seq.replace("C", "Y")
                reads2_seq = reads2_seq.replace("c", "Y")
                reads2_seq = reads2_seq.replace("G", "C")
                reads2_seq = reads2_seq.replace("g", "C")
                reads2_seq = reads2_seq.replace("Y", "G")  ###互补
                if circ_2[4] == "+":
                    read_num += 1
                    reads_1.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads1_seq) + '\n')
                    reads_1.write("+"+'\n')
                    reads_1.write("!" * len(reads1_seq))
                    reads_1.write('\n')
                    reads_2.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads2_seq) + '\n')
                    reads_2.write("+"+'\n')
                    reads_2.write("!" * len(reads2_seq))
                    reads_2.write('\n')
                    length.pop(0)
                else:
                    read_num += 1
                    reads_1.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads2_seq) + '\n')
                    reads_1.write("+"+'\n')
                    reads_1.write("!" * len(reads2_seq))
                    reads_1.write('\n')
                    reads_2.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads1_seq) + '\n')
                    reads_2.write("+"+'\n')
                    reads_2.write("!" * len(reads1_seq))
                    reads_2.write('\n')
                    length.pop(0)
                count += 1
            liner_information.append(circ_2[0])
            liner_information.append(str(int(circ_2[1])))
            liner_information.append(circ_2[2])
            liner_information.append(circ_2[4])
            liner_information.append(circ_2[3])
            liner_information.append(str(seq_length))
            liner_information.append(str(junction_num))
            informa_2 = '\t'.join(liner_information)
            circ_list.write(informa_2)
            circ_list.write('\n')
            circ_2.append(str(junction_num))
            circ_2.append(str(liner_time))
            iso_infor = '\t'.join(circ_2)
            circ_isoform_infor.write(iso_infor)
            circ_isoform_infor.write('\n')

        if circ_2[-2] == 'circ':
            junction_record.write(circ_2[0] + ":" + circ_2[1] + "|" + circ_2[2] + '\t')
            circ_exon_len = circ_2[6].split(',')
            circ_exon_start = circ_2[7].split(',')
            exon_num_2 = circ_2[5]
            junction_num = 0
            circ_information_2 = []
            circ_seq_1 = ""
            circ_chr_num = circ_2[0].split("chr")
            alseq = chromosome[circ_chr_num[1]]
            for exon_2 in range(1, int(exon_num_2) + 1):
                k = int(exon_2)
                circ_seq_1 = circ_seq_1 + alseq[(int(circ_2[1]) + int(circ_exon_start[k - 1]) + len(circ_2[0])):(int(circ_2[1]) + int(circ_exon_start[k - 1]) + int(circ_exon_len[k - 1]) + len(circ_2[0]))]  ##  为什么要加len(circ_2[0])？ 这是是把所有exon上的序列合在一起
            circ_seq = circ_seq_1
            seq_length = len(circ_seq)
            circ_time = int(float(circ_2[-1])) * (seq_length) * int(inputTotalNumber) / 1000000000  ##10亿，inputTotalNumber = 20000000。表达量 * isoform长度 * 测序深度 = 取样数目
            count = 1
            while count <= circ_time:
                junction_mark = 0
                start_loci = random.randint(0, seq_length - 1)
                if length == []:
                    length_array = 60 * numpy.random.randn(1, int((int(inputTotalNumber) / 10))) + (int(inputInsertLength) + (2 * int(inputReadLength)))  ##？
                    length = list(length_array[0])
                    length_array = []
                if start_loci + int(inputReadLength) <= seq_length:
                    reads1_seq = circ_seq[start_loci:(start_loci + int(inputReadLength))]
                elif seq_length < start_loci + int(inputReadLength) < 2 * seq_length:
                    reads1_seq = circ_seq[start_loci:int(seq_length)] + circ_seq[0:(start_loci + int(inputReadLength) - int(seq_length))]
                    junction_num += 1
                    junction_mark = 1
                elif start_loci + int(inputReadLength) > 2 * seq_length:
                    reads1_seq = circ_seq[start_loci:int(seq_length)] + circ_seq + circ_seq[0:(start_loci + int(inputReadLength) - 2 * int(seq_length))]
                    junction_num += 1
                    junction_mark = 1
                while int(length[0]) < int(inputReadLength):
                    length.pop(0)
                start_loci_2 = int(start_loci + int(length[0]) - int(inputReadLength))
                while start_loci_2 > seq_length:
                    start_loci_2 = start_loci_2 - seq_length
                if start_loci_2 + int(inputReadLength) <= seq_length:
                    reads2_seq = circ_seq[start_loci_2:(start_loci_2 + int(inputReadLength))]
                elif seq_length < start_loci_2 + int(inputReadLength) <= 2 * seq_length:
                    reads2_seq = circ_seq[start_loci_2:int(seq_length)] + circ_seq[0:(start_loci_2 + int(inputReadLength) - int(seq_length))]
                    junction_num += 1
                    junction_mark = 1
                elif start_loci_2 + int(inputReadLength) > 2 * seq_length:
                    reads2_seq = circ_seq[start_loci_2:int(seq_length)] + circ_seq + circ_seq[0:(start_loci_2 + int(inputReadLength) - 2 * int(seq_length))]
                    junction_num += 1
                    junction_mark = 1
                reads2_seq = reads2_seq[::-1]  ##取反
                reads2_seq = reads2_seq.replace("A", "X")
                reads2_seq = reads2_seq.replace("a", "X")
                reads2_seq = reads2_seq.replace("T", "A")
                reads2_seq = reads2_seq.replace("t", "A")
                reads2_seq = reads2_seq.replace("X", "T")
                reads2_seq = reads2_seq.replace("C", "Y")
                reads2_seq = reads2_seq.replace("c", "Y")
                reads2_seq = reads2_seq.replace("G", "C")
                reads2_seq = reads2_seq.replace("g", "C")
                reads2_seq = reads2_seq.replace("Y", "G")  ###互补
                if circ_2[4] == "+":
                    read_num += 1
                    if junction_mark == 1:
                        junction_record.write("@>simulate:" + str(read_num) + ",")
                    reads_1.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads1_seq) + '\n')
                    reads_1.write("+"+'\n')
                    reads_1.write("!" * len(reads1_seq))
                    reads_1.write('\n')
                    reads_2.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads2_seq) + '\n')
                    reads_2.write("+"+'\n')
                    reads_2.write("!" * len(reads2_seq))
                    reads_2.write('\n')
                    length.pop(0)
                else:
                    read_num += 1
                    if junction_mark == 1:
                        junction_record.write("@>simulate:" + str(read_num) + ",")
                    reads_1.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads2_seq) + '\n')
                    reads_1.write("+"+'\n')
                    reads_1.write("!" * len(reads2_seq))
                    reads_1.write('\n')
                    reads_2.write("@>simulate:" + str(read_num) + '\t' + "length=" + str(inputReadLength) + '\n' + str(reads1_seq) + '\n')
                    reads_2.write("+"+'\n')
                    reads_2.write("!" * len(reads1_seq))
                    reads_2.write('\n')
                    length.pop(0)

                count += 1
            junction_record.write("\n")
            circ_information_2.append(circ_2[0])
            circ_information_2.append(str(int(circ_2[1])))
            circ_information_2.append(circ_2[2])
            circ_information_2.append(circ_2[4])
            circ_information_2.append(circ_2[3])
            circ_information_2.append(str(seq_length))
            circ_information_2.append(str(junction_num))
            informa_2 = '\t'.join(circ_information_2)
            circ_list.write(informa_2)
            circ_list.write('\n')
            circ_2.append(str(junction_num))
            circ_2.append(str(circ_time))
            iso_infor = '\t'.join(circ_2)
            circ_isoform_infor.write(iso_infor)
            circ_isoform_infor.write('\n')
            if(str(circ_2[1]) + '-' + str(circ_2[2]) in ciri_dict_name):
                circ_dict[str(circ_2[1]) + '-' + str(circ_2[2])] = int(circ_dict[str(circ_2[1]) + '-' + str(circ_2[2])] + int(junction_num/2))
reads_1.close()
reads_2.close()
circ_list.close()
circ_isoform_infor.close()
junction_record.close()

i = 0
in_txt = open("ciri_input.txt", "r")
out_txt = open("ciri_input_new.txt", "a")
for line in in_txt:
    i = i+1
    line = line.strip('\n')
    if i < 2:
        out_txt.write(line + "\n")
        continue
    else:
        array = line.split()
        if (str(array[2]) + '-' + str(array[3])) in ciri_dict_name and circ_dict[str(array[2]) + '-' + str(array[3])] != 0.0:
            circ_info = []
            circ_info = array
            circ_info[4] = str(circ_dict[str(array[2]) + '-' + str(array[3])])
            informa = '\t'.join(circ_info)
            out_txt.write(informa + "\n")
in_txt.close()
out_txt.close()

cc = open("circRNA_list_before.txt","r")

result = open("circRNA_list_result.txt","a")

circ_info={}

for i in cc :
	j=i.strip('\n')
	j=j.split('\t')
	kk=j[0]+"*"+j[1]+"*"+j[2]+"*"+j[3]+"*"+j[4]
	if kk not in circ_info:
		circ_info[kk]=j[6]
	else:
		circ_info[kk]=str(int(circ_info[kk])+int(j[6]))

for key in circ_info:
	readsss_num=circ_info[key]
	info=key
	circ_id=info.split('*')
	readsss_len=int(circ_id[2])-int(circ_id[1])
	circRNA_id='\t'.join(circ_id)
	result.write(circRNA_id)
	result.write('\t')
	result.write(str(readsss_len))
	result.write('\t')
	result.write(str(readsss_num))
	result.write('\n')

cc.close()
result.close()

if os.path.exists('circRNA_list_before.txt'):
	os.remove('circRNA_list_before.txt')
if os.path.exists('circRNA_information_1.txt'):
	os.remove('circRNA_information_1.txt')
if os.path.exists('circRNA_information_before.txt'):
	os.remove('circRNA_information_before.txt')
if os.path.exists('simulate_temp/'):
	shutil.rmtree("simulate_temp/")

