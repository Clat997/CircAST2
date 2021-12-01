#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#by hongwen 2021-5-25

'''determine whether the version of user's python comply with the requirements of  this procedure'''
import sys
if sys.version_info[0] != 3 or sys.version_info[1] != 8:
    print(sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " CircAST2.0 needs python3.8!\n")
    sys.exit()


import time
import os
import optparse
import re
import shutil
import numpy as np
import copy
import collections
from difflib import SequenceMatcher
import difflib
import multiprocessing as mul
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from functools import partial
#==============================================================================
# functions definition
#==============================================================================
def sorted_or_not(infile):       
    '''determine if the input SAM file is sorted by name'''
    in_file = open(infile)
    flag = 0
    frag_name = []
    for line in in_file:
        line = line.strip('\n') 
        array = line.split()
        if array[0][0] == "@":
            continue
        else:
            if array[6] == "=":
                frag_name.append(array[0])
            if len(frag_name) == 6:
                break
    in_file.close()
    if frag_name[0] == frag_name[1] and frag_name[2] == frag_name[3] and frag_name[4] == frag_name[5]:
        flag = 1        
    return flag
def findname1(chr,a):     
    gname1 = []
    in_table = open(tempfile_path + chr + ".gtf","r")
    for line in in_table:
        line = line.strip('\n')
        array = line.split()
        if array[3] == a:
            t = array.index("gene_id")
            temp_name = array[t + 1][1:-2]
            if temp_name not in gname1:
                gname1.append(temp_name)
    in_table.close()
    return gname1 
def findname2(chr,b):  
    gname2 = []
    in_table = open(tempfile_path + chr + ".gtf","r")
    for line in in_table:
        line = line.strip('\n')
        array = line.split()
        if array[4] == b:
            t = array.index("gene_id")
            temp_name = array[t + 1][1:-2]
            if temp_name not in gname2:
                gname2.append(temp_name)
    in_table.close()
    return gname2 
def reverse(reads2_seq):
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
    return reads2_seq
def convert_flag(array):
    array[1] = str(int(array[1]) + 2)
    t = '\t'
    line = t.join(array)
    return line

    gname2 = []
    in_table = open(tempfile_path + chr + ".gtf", "r")
    for line in in_table:
        line = line.strip('\n')
        array = line.split()
        if array[4] == b:
            t = array.index("gene_id")
            temp_name = array[t + 1][1:-2]
            if temp_name not in gname2:
                gname2.append(temp_name)
    in_table.close()
    return gname2
def calculate(gene,gene_name,chr_no,match_num_total,head_num,inputfrag_len,tempfile_path):
    gene_id = gene
    i = gene_name.index(gene_id)
    chr_name = chr_no[i]
    bound = 3
    def max_match(A):
        dim = A.shape
        m = dim[0]
        n = dim[1]
        M = np.zeros((m,n))
        for i in range(0,m):
            for j in range(0,n):
                if A[i,j] == 1:
                    M[i,j] = 1
                    break
            if np.sum(M) > 0:
                break

        while 1:
            x = []
            for i in range(0,m):
                x.append(0)
            y = []
            for j in range(0,n):
                y.append(0)

            for i in range(0,m):
                pd = 1
                for j in range(0,n):
                    if M[i,j] == 1:
                        pd = 0
                if pd == 1:
                    x[i] = - n 

            pd = 0
            while 1:
                xi = -1
                for i in range(0,m):
                    if x[i] < 0:
                        xi = i
                        break
                if xi == -1:
                    pd = 1
                    break
                x[xi] = x[xi] * (-1)
                yy = []
                for j in range(0,n):
                    if A[xi,j] == 1 and y[j] == 0:
                        y[j] = xi
                        yy.append(j)

                if len(yy) >= 1:
                    for j in range(0,len(yy)):
                        pdd = 1
                        for i in range(0,m):
                            if M[i,yy[j]] == 1:
                                x[i] = -yy[j]
                                pdd = 0
                                break
                        if pdd == 1:
                            break
                    if pdd == 1:
                        j = yy[j]
                        P1 = []
                        P2 = []
                        while 1:
                            P2.append(j)
                            P1.append(y[j])
                            j = abs(x[y[j]])
                            if j == n:
                                break
                        for i in range(0,len(P1)):
                            if M[P1[i],P2[i]] == 1:
                                M[P1[i],P2[i]] = 0
                            else:
                                M[P1[i],P2[i]] = 1
                        break
            if pd == 1:
                break
        return M
    def findPath(Graph,partialPath,destination):  
        pathLength = len(partialPath)
        lastNode = partialPath[-1]
        possiablePaths = []
        nextNodes = []
        for i in range(lastNode,destination + 1):
            if Graph[lastNode,i] == 1:
                nextNodes.append(i)
        for i in range(0,len(nextNodes)):
            if destination == nextNodes[i]:
                tmpPath = partialPath + [destination]
                possiablePaths.append(tmpPath)
                nextNodes[i] = 0
        for i in range(len(nextNodes) - 1,-1,-1):
            if nextNodes[i] == 0:
                nextNodes.pop(i)
        for i in range(0,len(nextNodes)):
            tmpPath = partialPath + [nextNodes[i]]
            tmpPsbPaths = findPath(Graph, tmpPath, destination)     
            possiablePaths = possiablePaths + tmpPsbPaths
        return possiablePaths
    def quantification(matrix,isoform_len):
        k = 0
        in_sam = open(tempfile_path + gene_id + "_A.sam","r")    
        while k < head_num:
            line1 = in_sam.readline()
            k += 1
        k = 0
        while k != sam_row_num_A:
            line1 = in_sam.readline()
            line1 = line1.strip('\n')
            array1 = line1.split()
            data11 = re.findall("\d+",array1[5])
            data12 = re.findall("\D",array1[5])
            line2 = in_sam.readline()
            line2 = line2.strip('\n')
            array2 = line2.split()
            data21 = re.findall("\d+",array2[5])
            data22 = re.findall("\D",array2[5])
            count_array = np.zeros(len(matrix))
            frag_len_isoform = np.zeros(len(matrix))
            if data12.count("N") == 0 and data22.count("N") == 0:    
                match1 = 0
                for i in range(0,len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match1 = match1 + int(data11[i])
                match2 = 0
                for i in range(0,len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match2 = match2 + int(data21[i])
                t11 = int(array1[3])
                t12 = t11 + match1 - 1
                t21 = int(array2[3])
                t22 = t21 + match2 - 1
                for i in range(0,len(matrix)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0,len(matrix[i])):
                        if t11 >= int(matrix[i][j][0]) - 2 and t12 <= int(matrix[i][j][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix[i][j][1]) - t11 + 1
                            break
                    if t > -1: 
                        if t21 >= int(matrix[i][t][0]) - 2 and t22 <= int(matrix[i][t][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = t21 - t11 + match2
                        else:
                            for j in range(t + 1,len(matrix[i])):
                                if t21 >= int(matrix[i][j][0]) - 2 and t22 <= int(matrix[i][j][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + t22 - int(matrix[i][j][0]) + 1
                                    break
                                else:
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) -int(matrix[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000
            elif data12.count("N") == 0 and data22.count("N") == 1:     
                match1 = 0
                for i in range(0,len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match1 = match1 + int(data11[i])                        
                t11 = int(array1[3])
                t12 = t11 + match1 - 1                      
                N_id2 = data22.index("N")
                match21 = 0
                for i in range(0,N_id2):
                    if data22[i] == "M" or data22[i] == "D":
                        match21 = match21 + int(data21[i])                      
                match22 = 0
                for i in range(N_id2 + 1,len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match22 = match22 + int(data21[i])
                t21 = int(array2[3])
                t22 = t21 + match21 - 1
                t23 = t22 + int(data21[N_id2]) + 1
                t24 = t23 + match22 - 1     
                for i in range(0,len(matrix)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0,len(matrix[i]) -1):                    
                        if t11 >= int(matrix[i][j][0]) - 2 and t12 <= int(matrix[i][j][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix[i][j][1]) - t11 + 1
                            break
                    if t > -1:
                        if t21 >= int(matrix[i][t][0]) - 2 and abs(t22 - int(matrix[i][t][1])) <= 2 and abs(t23 - int(matrix[i][t + 1][0])) <= 2 and t24 <= int(matrix[i][t + 1][1]) + 2:
                            match[1]= 1
                            frag_len_isoform[i] = int(matrix[i][t][1]) - t11 + 1 + match22
                        else:
                            for j in range(t + 1,len(matrix[i]) - 1):
                                if t21 >= int(matrix[i][j][0]) -2 and abs(t22 - int(matrix[i][j][1])) <= 2 and abs(t23 - int(matrix[i][j + 1][0])) <= 2 and t24 <= int(matrix[i][j + 1][1]) + 2:
                                    match[1]= 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) - int(matrix[i][j][0]) + 1 + match22
                                    break
                                else:
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) -int(matrix[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000

            elif data12.count("N") == 1 and data22.count("N") == 0:   
                N_id1 = data12.index("N")
                match11 = 0
                for i in range(0,N_id1):
                    if data12[i] == "M" or data12[i] == "D":
                        match11 = match11 + int(data11[i])                      
                match12 = 0
                for i in range(N_id1 + 1,len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match12 = match12 + int(data11[i])
                t11 = int(array1[3])
                t12 = t11 + match11 - 1
                t13 = t12 + int(data11[N_id1]) + 1
                t14 = t13 + match12 - 1
                match2 = 0
                for i in range(0,len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match2 = match2 + int(data21[i])                        
                t21 = int(array2[3])
                t22 = t21 + match2 - 1

                for i in range(0,len(matrix)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0,len(matrix[i]) - 1):                   
                        if t11 >= int(matrix[i][j][0]) - 2 and abs(t12 - int(matrix[i][j][1])) <= 2 and abs(t13 - int(matrix[i][j + 1][0])) <= 2 and t14 <= int(matrix[i][j + 1][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix[i][j][1]) - t11 + 1
                            break
                    if t > -1:
                        for j in range(t + 1,len(matrix[i])):
                            if t21 >= int(matrix[i][j][0]) - 2 and t22 <= int(matrix[i][j][1]) + 2:
                                match[1] = 1
                                frag_len_isoform[i] = frag_len_isoform[i] + t21 - int(matrix[i][j][0]) + match2 
                                break
                            else:
                                frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) -int(matrix[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000

            elif data12.count("N") == 1 and data22.count("N") == 1:     
                N_id1 = data12.index("N")
                match11 = 0
                for i in range(0,N_id1):
                    if data12[i] == "M" or data12[i] == "D":
                        match11 = match11 + int(data11[i])                      
                match12 = 0
                for i in range(N_id1 + 1,len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match12 = match12 + int(data11[i])
                t11 = int(array1[3])
                t12 = t11 + match11 - 1
                t13 = t12 + int(data11[N_id1]) + 1
                t14 = t13 + match12 - 1
                N_id2 = data22.index("N")
                match21 = 0
                for i in range(0,N_id2):
                    if data22[i] == "M" or data22[i] == "D":
                        match21 = match21 + int(data21[i])                      
                match22 = 0
                for i in range(N_id2 + 1,len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match22 = match22 + int(data21[i])
                t21 = int(array2[3])
                t22 = t21 + match21 - 1
                t23 = t22 + int(data21[N_id2]) + 1
                t24 = t23 + match22 - 1     
        
                for i in range(0,len(matrix)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0,len(matrix[i]) - 2):
                        if t11 >= int(matrix[i][j][0]) - 2 and abs(t12 - int(matrix[i][j][1])) <= 2 and abs(t13 - int(matrix[i][j + 1][0])) <= 2 and t14 <= int(matrix[i][j + 1][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix[i][j][1]) - t11 + 1 + int(matrix[i][j + 1][1]) - int(matrix[i][j + 1][0]) + 1
                            break
                    if t > -1 and t + 1 < len(matrix[i]):
                        t = t + 1
                        if t21 >= int(matrix[i][t - 1][0]) - 2 and abs(t22 - int(matrix[i][t - 1][1])) <= 2 and abs(t23 - int(matrix[i][t][0])) <= 2 and t24 <= int(matrix[i][t][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = int(matrix[i][t - 1][1]) - t11 + 1 + match22 
                        elif t21 >= int(matrix[i][t][0]) - 2 and abs(t22 - int(matrix[i][t][1])) <= 2 and abs(t23 - int(matrix[i][t + 1][0])) <= 2 and t24 <= int(matrix[i][t + 1][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = frag_len_isoform[i] + match22
                        else:
                            for j in range(t + 1,len(matrix[i]) - 1):
                                if t21 >= int(matrix[i][j][0]) and abs(t22 - int(matrix[i][j][1])) <= 2 and abs(t23 - int(matrix[i][j + 1][0])) <= 2 and t24 <= int(matrix[i][j + 1][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) - int(matrix[i][j][0]) + match22
                                    break
                                else:
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) -int(matrix[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000
            else:
                for i in range(0,len(matrix)):
                    frag_len_isoform[i] = 1000000
            if sum(count_array) > 0:
                match_count = open(tempfile_path + gene_id + "_match.txt","a")       #a属性，每次打开后接着写，所以用完应该删除
                for i in range(0,len(matrix)):
                    match_count.write(str(count_array[i]) + "\t")
                match_count.write("\n")
                match_count.close()
                fraglen = open(tempfile_path + gene_id + "_fraglen.txt","a") 
                for i in range(0,len(matrix)):
                    len_adjust = frag_len_isoform[i] + isoform_len[i]
                    if abs(len_adjust - inputfrag_len) < abs(frag_len_isoform[i] - inputfrag_len):
                        frag_len_isoform[i] = len_adjust
                    fraglen.write(str(frag_len_isoform[i]) + "\t")
                fraglen.write("\n")
                fraglen.close()
            k = k + 2
        in_sam.close()
        
        if sam_row_num_B > 0:
            k = 0
            in_sam = open(tempfile_path + gene_id + "_B.sam","r")    
            while k < head_num:
                line1 = in_sam.readline()
                k += 1
            k = 0
            while k != sam_row_num_B:
                line1 = in_sam.readline()
                line1 = line1.strip('\n')
                array1 = line1.split()
                data11 = re.findall("\d+",array1[5])
                data12 = re.findall("\D",array1[5])
                line2 = in_sam.readline()
                line2 = line2.strip('\n')
                array2 = line2.split()
                data21 = re.findall("\d+",array2[5])
                data22 = re.findall("\D",array2[5])
                count_array = np.zeros(len(matrix))
                frag_len_isoform = np.zeros(len(matrix))
                if data12.count("N") == 0 and data22.count("N") == 0:
                    match1 = 0
                    for i in range(0,len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match1 = match1 + int(data11[i])
                    match2 = 0
                    for i in range(0,len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match2 = match2 + int(data21[i])
                    t11 = int(array1[3])
                    t12 = t11 + match1 - 1
                    t21 = int(array2[3])
                    t22 = t21 + match2 - 1
                    for i in range(0,len(matrix)):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0,len(matrix[i])):
                            if t11 >= int(matrix[i][j][0]) - 2 and t12 <= int(matrix[i][j][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix[i][j][0])
                                for m in range(0,t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1
                                break
                        if t > -1: 
                            for j in range(t,len(matrix[i])):
                                if t21 >= int(matrix[i][j][0]) - 2 and t22 <= int(matrix[i][j][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) - t22
                                    for m in range(j + 1,len(matrix[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match1 + match2
                                    break                           
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000
                            
                elif data12.count("N") == 0 and data22.count("N") == 1:
                    match1 = 0
                    for i in range(0,len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match1 = match1 + int(data11[i])
                    t11 = int(array1[3])
                    t12 = t11 + match1 - 1
                    N_id2 = data22.index("N")
                    match21 = 0
                    for i in range(0,N_id2):
                        if data22[i] == "M" or data22[i] == "D":
                            match21 = match21 + int(data21[i])
                    match22 = 0
                    for i in range(N_id2 + 1,len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match22 = match22 + int(data21[i])
                    t21 = int(array2[3])
                    t22 = t21 + match21 - 1
                    t23 = t22 + int(data21[N_id2]) + 1
                    t24 = t23 + match22 - 1     
                    for i in range(0,len(matrix)):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0,len(matrix[i]) -1):                    
                            if t11 >= int(matrix[i][j][0]) - 2 and t12 <= int(matrix[i][j][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix[i][j][0])
                                for m in range(0,t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1
                                break
                        if t > -1:
                            for j in range(t,len(matrix[i]) - 1):
                                if t21 >= int(matrix[i][j][0]) - 2 and abs(t22 - int(matrix[i][j][1])) <= 2 and abs(t23 - int(matrix[i][j + 1][0])) <= 2 and t24 <= int(matrix[i][j + 1][1]) + 2:
                                    match[1]= 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j + 1][1]) - t24
                                    for m in range(j + 2,len(matrix[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match1 + match21 + match22
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000

                elif data12.count("N") == 1 and data22.count("N") == 0:    
                    N_id1 = data12.index("N")
                    match11 = 0
                    for i in range(0,N_id1):
                        if data12[i] == "M" or data12[i] == "D":
                            match11 = match11 + int(data11[i])
                    match12 = 0
                    for i in range(N_id1 + 1,len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match12 = match12 + int(data11[i])
                    t11 = int(array1[3])
                    t12 = t11 + match11 - 1
                    t13 = t12 + int(data11[N_id1]) + 1
                    t14 = t13 + match12 - 1
                    match2 = 0
                    for i in range(0,len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match2 = match2 + int(data21[i])
                    t21 = int(array2[3])
                    t22 = t21 + match2 - 1

                    for i in range(0,len(matrix)):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0,len(matrix[i]) - 1):                   
                            if t11 >= int(matrix[i][j][0]) - 2 and abs(t12 - int(matrix[i][j][1])) <= 2 and abs(t13 - int(matrix[i][j + 1][0])) <= 2 and t14 <= int(matrix[i][j + 1][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix[i][j][0])
                                for m in range(0,t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1 
                                break
                        if t > -1:
                            for j in range(t + 1,len(matrix[i])):
                                if t21 >= int(matrix[i][j][0]) - 2 and t22 <= int(matrix[i][j][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j][1]) - t22
                                    for m in range(j + 1,len(matrix[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match11 + match12 + match2
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000

                elif data12.count("N") == 1 and data22.count("N") == 1:     
                    N_id1 = data12.index("N")
                    match11 = 0
                    for i in range(0,N_id1):
                        if data12[i] == "M" or data12[i] == "D":
                            match11 = match11 + int(data11[i])                      
                    match12 = 0
                    for i in range(N_id1 + 1,len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match12 = match12 + int(data11[i])
                    t11 = int(array1[3])
                    t12 = t11 + match11 - 1
                    t13 = t12 + int(data11[N_id1]) + 1
                    t14 = t13 + match12 - 1
                    N_id2 = data22.index("N")
                    match21 = 0
                    for i in range(0,N_id2):
                        if data22[i] == "M" or data22[i] == "D":
                            match21 = match21 + int(data21[i])
                    match22 = 0
                    for i in range(N_id2 + 1,len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match22 = match22 + int(data21[i])
                    t21 = int(array2[3])
                    t22 = t21 + match21 - 1
                    t23 = t22 + int(data21[N_id2]) + 1
                    t24 = t23 + match22 - 1
            
                    for i in range(0,len(matrix)):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0,len(matrix[i]) - 2):
                            if t11 >= int(matrix[i][j][0]) - 2 and abs(t12 - int(matrix[i][j][1])) <= 2 and abs(t13 - int(matrix[i][j + 1][0])) <= 2 and t14 <= int(matrix[i][j + 1][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix[i][j][0])
                                for m in range(0,t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1 
                                break
                        if t > -1 and t + 1 < len(matrix[i]):
                            for j in range(t,len(matrix[i]) - 1):
                                if t21 >= int(matrix[i][j][0]) and abs(t22 - int(matrix[i][j][1])) <= 2 and abs(t23 - int(matrix[i][j + 1][0])) <= 2 and t24 <= int(matrix[i][j + 1][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][j + 1][1]) - t24
                                    for m in range(j + 2,len(matrix[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix[i][m][1]) - int(matrix[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match11 + match12 + match21 + match22
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000
                else:
                    for i in range(0,len(matrix)):
                        frag_len_isoform[i] = 1000000
                if sum(count_array) > 0:
                    match_count = open(tempfile_path + gene_id + "_match.txt","a")  
                    for i in range(0,len(matrix)):
                        match_count.write(str(count_array[i]) + "\t")
                    match_count.write("\n")
                    match_count.close()
                    fraglen = open(tempfile_path + gene_id + "_fraglen.txt","a")                        
                    for i in range(0,len(matrix)):
                        len_adjust = frag_len_isoform[i] + isoform_len[i]
                        if abs(len_adjust - inputfrag_len) < abs(frag_len_isoform[i] - inputfrag_len):
                            frag_len_isoform[i] = len_adjust
                        fraglen.write(str(frag_len_isoform[i]) + "\t")                                              
                    fraglen.write("\n")
                    fraglen.close()
                k = k + 2               
            in_sam.close()

        FPKM = np.zeros(len(matrix))
        if os.path.exists(tempfile_path + gene_id + "_fraglen.txt"):        
            y = np.loadtxt(tempfile_path + gene_id + "_fraglen.txt")   # 
            dim = np.shape(y)
            if len(dim) == 1: 
                if len(matrix) == 1: 
                    for i in range(0,len(matrix)):
                        for k in range(0,circ_spot_num):
                            if circ_spot[k][0] == matrix[i][0][0] and circ_spot[k][1] == matrix[i][-1][-1]:
                                FPKM[i] = 10 ** 9 * (circ_num[k] + dim[0]) / (isoform_len[i] * match_num_total * 1.0)
                                FPKM[i] = float('%.4f' % FPKM[i])

            elif len(dim) != 0:         
                frag_len = []
                for i in range(0,dim[0]):
                    for j in range(0,len(matrix)):
                        if y[i,j] != 1000000:
                            frag_len.append(y[i,j])
                miu = np.mean(frag_len)   
                sigma_square = np.var(frag_len)  
                sigma = np.sqrt(sigma_square)   
                if sigma != 0:
                    for i in range(0,circ_spot_num):
                        count_array = np.zeros(len(matrix))
                        frag_len_isoform = np.zeros(len(matrix))
                        for k in range(0,len(matrix)):
                            if circ_spot[i][0] == matrix[k][0][0] and circ_spot[i][1] == matrix[k][-1][-1]:
                                count_array[k] = 1
                                frag_len_isoform[k] = int(miu)
                            else:
                                frag_len_isoform[k] = 1000000   
                        if sum(count_array) > 0:                                        
                            match_count = open(tempfile_path + gene_id + "_match.txt","a")      
                            for j in range(0,circ_num[i]):  
                                for t in range(0,len(matrix)):
                                    match_count.write(str(count_array[t]) + "\t")
                                match_count.write("\n")
                            match_count.close()
                            fraglen = open(tempfile_path + gene_id + "_fraglen.txt","a")
                            for j in range(0,circ_num[i]):  
                                for t in range(0,len(matrix)):
                                    fraglen.write(str(frag_len_isoform[t]) + "\t")
                                fraglen.write("\n")
                            fraglen.close()
                    new_frag = np.loadtxt(tempfile_path + gene_id + "_fraglen.txt")
                    match_num = len(new_frag)

                    lt = np.array(isoform_len)
                    lu = np.zeros(len(matrix))
                    lu_iter = np.zeros(len(matrix))
                    error = np.zeros(len(matrix))
                    F = np.zeros((match_num,len(matrix)))
                    for j in range(0,len(matrix)):        
                        lu[j] =  1.0 / len(matrix)

                    for i in range(0,match_num):
                        for j in range(0,len(matrix)):
                            if new_frag[i,j] == 1000000:
                                F[i,j] = 0
                            else:                       
                                F[i,j] = np.exp((new_frag[i,j] - miu) ** 2 *(-0.5) / sigma_square) / (np.sqrt(2 * np.pi) * sigma)

                    flag = 1
                    flag1 = 0
                    z = np.zeros((match_num,len(matrix)))
                    iterate = 0
                    lul = np.zeros(len(matrix))
                    alpha = np.zeros(len(matrix))
                    alpha_iter = np.zeros(len(matrix))
                    al_lt = np.zeros(len(matrix))
                    z_fenzi = np.zeros((match_num,len(matrix)))
                    while flag == 1:
                        t = 0
                        for j in range(0,len(matrix)):
                            lul[j] = lu[j] * lt[j]
                        t = sum(lul)
                        for j in range(0,len(matrix)):
                            alpha[j] = lul[j] / (t * 1.0)
                        for i in range(0,match_num):
                            sum1 = 0
                            for j in range(0,len(matrix)):          
                                z_fenzi[i,j] = alpha[j] * F[i,j] / (lt[j] - bound * 2 + 2.0)
                                sum1 = sum1 + z_fenzi[i,j]
                            if sum1 != 0:
                                for j in range(0,len(matrix)):
                                    z[i,j] = z_fenzi[i,j] * 1.0 / sum1
                            else:
                                flag = 0
                                break
                        if flag == 1:
                            n = sum(z)
                            for j in range(0,len(matrix)):
                                alpha_iter[j] = n[j] * 1.0 / match_num
                            for j in range(0,len(matrix)):
                                al_lt[j] = alpha_iter[j] * 1.0 / lt[j]
                            k = sum(al_lt)
                            for j in range(0,len(matrix)):
                                lu_iter[j] = al_lt[j] * 1.0 / k
                                error[j] = abs(lu[j]-lu_iter[j])
                            iterate = iterate + 1
                            if sum(error) < 0.0001:
                                flag = 0
                                flag1 = 1
                            else:
                                for j in range(0,len(matrix)):
                                    lu[j] = lu_iter[j]
                    if flag1 == 1:
                        for j in range(0,len(matrix)):
                            FPKM[j] = 10 ** 9 * alpha_iter[j] * match_num / (lt[j] * match_num_total * 1.0)
                            FPKM[j] = float('%.4f' % FPKM[j])
            os.remove(tempfile_path + gene_id + "_match.txt")
            os.remove(tempfile_path + gene_id + "_fraglen.txt")
        return FPKM

    #######在gtf文件中找属于该染色体，该基因位点的exon条目，写入gene_id.gtf，记录基因的起始和中止坐标
    gene_coordinate_min = 0
    gene_coordinate_max = 0
    row = 0
    in_gtf = open(tempfile_path + chr_name + ".gtf")
    out_gtf = open(tempfile_path + gene_id + ".gtf", "w")
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        t = array.index("gene_id")
        temp_name = array[t + 1][1:-2]
        if array[0] == chr_name and array[2] == "exon" and temp_name == gene_id:
            out_gtf.write(line + "\n")
            row += 1
            if row == 1:
                gene_coordinate_min = int(array[3])
            if int(array[4]) > gene_coordinate_max:
                gene_coordinate_max = int(array[4])
    in_gtf.close()
    out_gtf.close()
    ##处理selected_circAST，先组装
    sam_row_num_A = 0
    i = 0
    in_sam = open(tempfile_path + "selected_" + chr_name + "_A.sam","r")
    out_sam = open(tempfile_path + gene_id + "_A.sam","a")
    for line in in_sam:
        line = line.strip('\n')
        if i < head_num:
            out_sam.write(line + "\n")
        else:
            k = i - head_num
            line = line.strip('\n')
            array = line.split()
            if k % 2 == 1:
                t1 = int(array[7])
                data = re.findall("\d+\D",array[5])
                s = 0
                for j in data:
                    if j[-1] in CIGAR_set1:
                        s = s + int(j[0:-1])
                t2 = int(array[3]) + s - 1
                if t1 >= gene_coordinate_min - 2 and t2 <= gene_coordinate_max + 2:
                    out_sam.write(line1 + "\n")
                    out_sam.write(line + "\n")
                    sam_row_num_A += 2
            line1 = line
        i = i + 1
    in_sam.close()
    out_sam.close()

    shutil.copyfile(tempfile_path + gene_id + "_A.sam", tempfile_path + gene_id + ".sam")               
    sam_row_num_B = 0
    i = 0
    in_sam = open(tempfile_path + "selected_" + chr_name + "_B.sam","r") 
    out_sam1 = open(tempfile_path + gene_id + "_B.sam","w")
    out_sam2 = open(tempfile_path + gene_id + ".sam","a")
    for line in in_sam:
        line = line.strip('\n')
        if i < head_num:
            out_sam1.write(line + "\n")
        else:
            k = i - head_num
            line = line.strip('\n')
            array = line.split()
            if k % 2 == 1:
                t1 = int(array[7])
                data = re.findall("\d+\D",array[5])
                s = 0
                for j in data:
                    if j[-1] in CIGAR_set1:
                        s = s + int(j[0:-1])
                t2 = int(array[3]) + s - 1
                if t1 >= gene_coordinate_min - 2 and t2 <= gene_coordinate_max + 2:
                    out_sam1.write(line1 + "\n")
                    out_sam1.write(line + "\n")
                    sam_row_num_B += 2
                    out_sam2.write(line1 + "\n")
                    out_sam2.write(line + "\n")
            line1 = line
        i = i + 1
    in_sam.close()
    out_sam1.close()
    out_sam2.close()
    
    ##组装
    circ_spot = []
    circ_num = []
    strand = []
    in_txt = open(tempfile_path + "circRNA_list_selected.txt","r")
    for line in in_txt:
        line = line.strip('\n')
        array = line.split()
        if array[0] == chr_name and array[4] == gene_id:
            temp =  [0] * 2
            temp[0] = int(array[1]) + 1 
            temp[1] = int(array[2])
            circ_spot.append(temp)
            circ_num.append(int(float(array[6])))
            strand.append(array[3])
    in_txt.close()
    
    circ_spot_num = len(circ_spot)
    sam_num = [0] * circ_spot_num     
    for t in range(0,circ_spot_num):
        i = 0
        sam_num[t] = 0
        in_sam = open(tempfile_path + gene_id + ".sam","r")
        out_sam = open(tempfile_path + gene_id + "_" + str(t) + ".sam","w")   
        for line in in_sam:         
            if i > head_num - 1:
                k = i - head_num
                line = line.strip('\n')
                array = line.split()
                if k % 2 == 1:      
                    t1 = int(array[7])
                    data = re.findall("\d+\D",array[5])
                    s = 0               
                    for j in data:
                        if j[-1] in CIGAR_set1:
                            s = s + int(j[0:-1])
                    t2 = int(array[3]) + s - 1
                    if t1 >= circ_spot[t][0] - 2 and t2 <= circ_spot[t][1] + 2:
                        out_sam.write(line1 + "\n")
                        out_sam.write(line + "\n")
                        sam_num[t] += 2
                line1 = line
            i = i + 1
        in_sam.close()
        out_sam.close()

    #construct DAGs
    exon_bin_0 = [None] * circ_spot_num
    exon_bin = [None] * circ_spot_num
    transcript_bin_idx = [None] * circ_spot_num
    transcript_bin = [None] * circ_spot_num
    subpath_constraints = [None] * circ_spot_num   
    edge = [None] * circ_spot_num
    new_edge = [None] * circ_spot_num
    isoform = [None] * circ_spot_num
    path_exon = [None] * circ_spot_num
    isoform = [None] * circ_spot_num
    isoform_num = [0] * circ_spot_num
    for t in range(0,circ_spot_num):
        exon_bin_0[t] = []
        transcript_bin[t] = []
        isoform[t] = []
        start = circ_spot[t][0]   
        end = circ_spot[t][1]    
        in_gtf = open(tempfile_path + gene_id + ".gtf","r")     
        for line in in_gtf:         
            line = line.strip('\n')
            array = line.split()
            if int(array[3]) >= start and int(array[4]) <= end:
                temp =  [0] * 2
                temp[0] =  int(array[3])
                temp[1] =  int(array[4])
                if temp not in exon_bin_0[t]:
                    exon_bin_0[t].append(temp)
        in_gtf.close()
        exon_bin_0[t].sort(key=lambda x:(x[0],x[1]))     
        exon_lb = []
        exon_rb = []
        for i in range(0,len(exon_bin_0[t])):
            exon_lb.append(exon_bin_0[t][i][0])
            exon_rb.append(exon_bin_0[t][i][1])
        if start in exon_lb and end in exon_rb: 
            while 1:
                if exon_bin_0[t][-1][-1] != end:
                    exon_bin_0[t].pop(-1)
                else:
                    break

            exon_bin[t] = copy.deepcopy(exon_bin_0[t])
                        

            if len(exon_bin[t]) == 1: 
                isoform[t] = [exon_bin[t]]
                isoform_num[t] = 1
            if len(exon_bin[t]) == 2:
                if exon_bin[t][0][0] == start and exon_bin[t][0][1] < exon_bin[t][1][0] and exon_bin[t][1][1] == end:
                    isoform[t] = [exon_bin[t]]
                    isoform_num[t] = 1
                elif exon_bin[t][0][0] == start and exon_bin[t][0][1] == end:
                    exon_bin[t].pop(1)
                    isoform[t] = [exon_bin[t]]
                    isoform_num[t] = 1
                elif exon_bin[t][1][0] == start and exon_bin[t][1][1] == end:
                    exon_bin[t].pop(0)
                    isoform[t] = [exon_bin[t]]
                    isoform_num[t] = 1

            if len(exon_bin[t]) >= 3 and exon_bin[t][0][0] < exon_bin[t][1][0] and exon_bin[t][-2][1] < exon_bin[t][-1][1]:
                node_num = len(exon_bin[t])
                ASG = np.zeros((node_num,node_num))
                subpath_constraints[t] = []
                in_sam = open(tempfile_path + gene_id + "_" + str(t) + ".sam","r")
                for line in in_sam:
                    line = line.strip('\n')
                    array = line.split()
                    data1 = re.findall("\d+",array[5])
                    data2 = re.findall("\D",array[5])
                    if data2.count("N") == 1:                       
                        N_id = data2.index("N")
                        match1 = 0
                        for i in range(0,N_id):
                            if data2[i] == "M" or data2[i] == "D":
                                match1 = match1 + int(data1[i])                     
                        match2 = 0
                        for i in range(N_id + 1,len(data1)):
                            if data2[i] == "M" or data2[i] == "D":
                                match2 = match2 + int(data1[i])
                        if match1 >= bound and match2 >= bound:    
                            t1 = int(array[3])
                            t2 = t1 + match1 - 1
                            t3 = t2 + int(data1[N_id]) + 1
                            t4 = t3 + match2 - 1
                            for i in range(0,node_num - 1):
                                if t1 >= exon_bin[t][i][0] - 2 and abs(t2 - exon_bin[t][i][1]) <= 2 :
                                    for j in range(i + 1,node_num):
                                        if abs(t3 - exon_bin[t][j][0]) <=2 and t4 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1                       
                    elif data2.count("N") == 2:     
                        N_id1 = data2.index("N")
                        data2[N_id1] = "X"
                        N_id2 = data2.index("N")
                        data2[N_id1] = "N"
                        match1 = 0
                        for i in range(0,N_id1):
                            if data2[i] == "M" or data2[i] == "D":
                                match1 = match1 + int(data1[i]) 
                        match2 = 0
                        for i in range(N_id1 + 1,N_id2):
                            if data2[i] == "M" or data2[i] == "D":
                                match2 = match2 + int(data1[i])
                        match3 = 0
                        for i in range(N_id2 + 1,len(data1)):
                            if data2[i] == "M" or data2[i] == "D":
                                match3 = match3 + int(data1[i])
                        t1 = int(array[3])
                        t2 = t1 + match1 - 1
                        t3 = t2 + int(data1[N_id1]) + 1
                        t4 = t3 + match2 - 1
                        t5 = t4 + int(data1[N_id2]) + 1
                        t6 = t5 + match3 - 1
                        if match1 < bound and match3 >= bound:      
                            for i in range(1,node_num - 1):
                                if abs(t3 - exon_bin[t][i][0]) <= 2 and abs(t4 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num):
                                        if abs(t5 - exon_bin[t][j][0]) <= 2 and t6 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1                           
                        elif match1 >= bound and match3 < bound:  
                            for i in range(0,node_num - 2):
                                if t1 >= exon_bin[t][i][0] - 2 and abs(t2 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t3 - exon_bin[t][j][0]) <= 2 and abs(t4 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                        elif match1 >= bound and match3 >= bound:  
                            for i in range(0,node_num - 2):
                                if t1 >= exon_bin[t][i][0] - 2 and abs(t2 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t3 - exon_bin[t][j][0]) <= 2 and abs(t4 - exon_bin[t][j][1]) <= 2:
                                            for k in range(j + 1,node_num):
                                                if abs(t5 - exon_bin[t][k][0]) <= 2 and t6 <= exon_bin[t][k][1] + 2:
                                                    ASG[i,j] += 1
                                                    ASG[j,k] += 1
                                                    temp =  [None] * 3
                                                    temp[0] =  exon_bin[t][i]
                                                    temp[1] =  exon_bin[t][j]
                                                    temp[2] =  exon_bin[t][k]
                                                    if temp not in subpath_constraints[t]:
                                                        subpath_constraints[t].append(temp)
                in_sam.close()
                m = 0
                in_sam = open(tempfile_path + gene_id + "_" + str(t) + ".sam","r")
                while m != sam_num[t]:
                    line1 = in_sam.readline()
                    line1 = line1.strip('\n')
                    array1 = line1.split()
                    data11 = re.findall("\d+",array1[5])
                    data12 = re.findall("\D",array1[5])
                    line2 = in_sam.readline()
                    line2 = line2.strip('\n')
                    array2 = line2.split()
                    data21 = re.findall("\d+",array2[5])
                    data22 = re.findall("\D",array2[5])
                    if data12.count("N") == 0 and data22.count("N") == 0:    
                        match1 = 0
                        for i in range(0,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match1 = match1 + int(data11[i])                        
                        match2 = 0
                        for i in range(0,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match2 = match2 + int(data21[i])                        
                        t11 = int(array1[3])
                        t12 = t11 + match1 - 1
                        t21 = int(array2[3])
                        t22 = t21 + match2 - 1
                        for i in range(0,node_num - 1):
                            if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2:
                                flag = 0
                                for j in range(i + 1,node_num):
                                    if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                        flag = 1
                                        break
                                if flag == 1:
                                    if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                    if j < node_num - 1:
                                        if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                            j = j + 1
                                            if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                                ASG[i,j] += 1
                                    
                    elif data12.count("N") == 1 and data22.count("N") == 1:         
                        N_id1 = data12.index("N")
                        match11 = 0
                        for i in range(0,N_id1):
                            if data12[i] == "M" or data12[i] == "D":
                                match11 = match11 + int(data11[i])                      
                        match12 = 0
                        for i in range(N_id1 + 1,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match12 = match12 + int(data11[i])
                        t11 = int(array1[3])
                        t12 = t11 + match11 - 1
                        t13 = t12 + int(data11[N_id1]) + 1
                        t14 = t13 + match12 - 1
                        N_id2 = data22.index("N")
                        match21 = 0
                        for i in range(0,N_id2):
                            if data22[i] == "M" or data22[i] == "D":
                                match21 = match21 + int(data21[i])                      
                        match22 = 0
                        for i in range(N_id2 + 1,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match22 = match22 + int(data21[i])
                        t21 = int(array2[3])
                        t22 = t21 + match21 - 1
                        t23 = t22 + int(data21[N_id2]) + 1
                        t24 = t23 + match22 - 1                     

                        if match11 >= bound and match22 >= bound:   
                            for i in range(0,node_num - 2):     
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            for k in range(j + 1,node_num):
                                                if abs(t23 - exon_bin[t][k][0]) <= 2 and t24 <= exon_bin[t][k][1] + 2:
                                                    temp =  [None] * 3
                                                    temp[0] =  exon_bin[t][i]
                                                    temp[1] =  exon_bin[t][j]
                                                    temp[2] =  exon_bin[t][k]
                                                    if temp not in subpath_constraints[t]:
                                                        subpath_constraints[t].append(temp)
                                                    if match12 < bound:
                                                        ASG[i,j] += 1
                                                    if match21 < bound:
                                                        ASG[j,k] += 1
                            if match12 >= bound and match21 >= bound: 
                                for i in range(0,node_num - 3):
                                    if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                        for j in range(i + 1,node_num - 2):
                                            if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2:
                                                flag = 0
                                                for k in range(j + 1,node_num - 1):
                                                    if exon_bin[t][k][0] >= exon_bin[t][j][1]:
                                                        flag = 1
                                                        break
                                                if flag == 1:                                           
                                                    if t21 >= exon_bin[t][k][0] - 2 and abs(t22 - exon_bin[t][k][1]) <= 2:
                                                        for l in range(k + 1,node_num):
                                                            if abs(t23 - exon_bin[t][l][0]) <= 2 and t24 <= exon_bin[t][l][1] + 2:
                                                                ASG[j,k] += 1
                                                                temp =  [None] * 4
                                                                temp[0] =  exon_bin[t][i]
                                                                temp[1] =  exon_bin[t][j]
                                                                temp[2] =  exon_bin[t][k]
                                                                temp[3] =  exon_bin[t][l]
                                                                if temp not in subpath_constraints[t]:
                                                                    subpath_constraints[t].append(temp)
                                                    if k < node_num - 2:
                                                        if exon_bin[t][k][0] == exon_bin[t][k + 1][0]:
                                                            k = k + 1
                                                            if t21 >= exon_bin[t][k][0] - 2 and abs(t22 - exon_bin[t][k][1]) <= 2:
                                                                for l in range(k + 1,node_num):
                                                                    if abs(t23 - exon_bin[t][l][0]) <= 2 and t24 <= exon_bin[t][l][1] + 2:
                                                                        ASG[j,k] += 1
                                                                        temp =  [None] * 4
                                                                        temp[0] =  exon_bin[t][i]
                                                                        temp[1] =  exon_bin[t][j]
                                                                        temp[2] =  exon_bin[t][k]
                                                                        temp[3] =  exon_bin[t][l]
                                                                        if temp not in subpath_constraints[t]:
                                                                            subpath_constraints[t].append(temp)         
                                                                            
                        if match12 < bound and match22 < bound:   
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    ASG[i,j] += 1                                           
                                    
                        if match11 < bound and match21 < bound:   
                            for i in range(1,node_num - 1):
                                if abs(t13 - exon_bin[t][i][0]) <= 2 and t14 <= exon_bin[t][i][1] + 2 and t21 >= exon_bin[t][i][0] - 2 and abs(t22 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1

                        if match11 < bound and match22 < bound:   
                            for i in range(1,node_num - 2):
                                if abs(t13 - exon_bin[t][i][0]) <= 2 and t14 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 <= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 <= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    ASG[i,j] += 1
                                        

                    elif data12.count("N") == 1 and data22.count("N") == 0:   
                        N_id1 = data12.index("N")
                        match11 = 0
                        for i in range(0,N_id1):
                            if data12[i] == "M" or data12[i] == "D":
                                match11 = match11 + int(data11[i])                      
                        match12 = 0
                        for i in range(N_id1 + 1,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match12 = match12 + int(data11[i])
                        t11 = int(array1[3])
                        t12 = t11 + match11 - 1
                        t13 = t12 + int(data11[N_id1]) + 1
                        t14 = t13 + match12 - 1
                        match2 = 0
                        for i in range(0,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match2 = match2 + int(data21[i])                        
                        t21 = int(array2[3])
                        t22 = t21 + match2 - 1
                        if match11 >= bound and match12 >= bound:
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2:
                                            flag = 0
                                            for k in range(j + 1,node_num):
                                                if exon_bin[t][k][0] >= exon_bin[t][j][1]:
                                                    flag = 1
                                                    break
                                            if flag == 1:
                                                if t21 >= exon_bin[t][k][0] - 2 and t22 <= exon_bin[t][k][1] + 2:
                                                        ASG[j,k] += 1   
                                                        temp =  [None] * 3
                                                        temp[0] =  exon_bin[t][i]
                                                        temp[1] =  exon_bin[t][j]
                                                        temp[2] =  exon_bin[t][k]
                                                        if temp not in subpath_constraints[t]:
                                                            subpath_constraints[t].append(temp)
                                                if k < node_num - 1:
                                                    if exon_bin[t][k][0] == exon_bin[t][k + 1][0]:
                                                        k = k + 1
                                                        if t21 >= exon_bin[t][k][0] - 2 and t22 <= exon_bin[t][k][1] + 2:
                                                            ASG[j,k] += 1   
                                                            temp =  [None] * 3
                                                            temp[0] =  exon_bin[t][i]
                                                            temp[1] =  exon_bin[t][j]
                                                            temp[2] =  exon_bin[t][k]
                                                            if temp not in subpath_constraints[t]:
                                                                subpath_constraints[t].append(temp)
                                                                                                
                        elif match11 >= bound and match12 < bound:
                            for i in range(0,node_num - 1):
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and  t22 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and  t22 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1
                                                                                
                        elif match11 < bound and match12 >= bound:
                            for i in range(1,node_num - 1):
                                if abs(t13 - exon_bin[t][i][0]) <= 2 and t14 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1   
                                                    
                    elif data12.count("N") == 0 and data22.count("N") == 1:     
                        match1 = 0
                        for i in range(0,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match1 = match1 + int(data11[i])                        
                        t11 = int(array1[3])
                        t12 = t11 + match1 - 1                      
                        N_id2 = data22.index("N")
                        match21 = 0
                        for i in range(0,N_id2):
                            if data22[i] == "M" or data22[i] == "D":
                                match21 = match21 + int(data21[i])                      
                        match22 = 0
                        for i in range(N_id2 + 1,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match22 = match22 + int(data21[i])
                        t21 = int(array2[3])
                        t22 = t21 + match21 - 1
                        t23 = t22 + int(data21[N_id2]) + 1
                        t24 = t23 + match22 - 1     

                        if match21 >= bound and match22 >= bound:
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            for k in range(j + 1,node_num):
                                                if abs(t23 - exon_bin[t][k][0]) <= 2 and t24 <= exon_bin[t][k][1] + 2:
                                                    ASG[i,j] += 1
                                                    temp =  [None] * 3
                                                    temp[0] =  exon_bin[t][i]
                                                    temp[1] =  exon_bin[t][j]
                                                    temp[2] =  exon_bin[t][k]
                                                    if temp not in subpath_constraints[t]:
                                                        subpath_constraints[t].append(temp)
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    for k in range(j + 1,node_num):
                                                        if abs(t23 - exon_bin[t][k][0]) <= 2 and t24 <= exon_bin[t][k][1] + 2:
                                                            ASG[i,j] += 1
                                                            temp =  [None] * 3
                                                            temp[0] =  exon_bin[t][i]
                                                            temp[1] =  exon_bin[t][j]
                                                            temp[2] =  exon_bin[t][k]
                                                            if temp not in subpath_constraints[t]:
                                                                subpath_constraints[t].append(temp)                                         
                                                    
                        elif match21 >= bound and match22 < bound:
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    ASG[i,j] += 1
                                                    
                        elif match21 < bound and match22 >= bound:
                            for i in range(0,node_num - 1):
                                if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2 and t21 >= exon_bin[t][i][0] - 2 and abs(t22 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1                                    
                                            break
                                    if flag == 1:
                                        if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1
                                    
                                    
                    m += 2
                in_sam.close()
                for i in range(0,node_num):
                    for j in range(0,node_num):
                        if ASG[i,j] < 3:
                            ASG[i,j] = 0
                        else:
                            ASG[i,j] = 1          
                
                
                while 1:
                    isolate = []
                    for i in range(1,node_num - 1):
                        if sum(ASG[:,i]) == 0 or sum(ASG[i,:]) == 0:
                            isolate.append(i)
                    if isolate != []:
                        isolate.reverse()
                        for i in isolate:
                            exon_bin[t].pop(i)      
                        ASG = np.delete(ASG,isolate,0)
                        ASG = np.delete(ASG,isolate,1)
                        node_num = len(exon_bin[t])
                    else:
                        break
                
                if sum(ASG[0,:]) > 0 and sum(ASG[:,-1]) > 0:                        
                    for i in range(len(subpath_constraints[t]) - 1,-1,-1):
                        for j in range(0,len(subpath_constraints[t][i])):
                            if subpath_constraints[t][i][j] not in exon_bin[t]:
                                subpath_constraints[t].pop(i)
                                break
                    node_num_not_isolated1 = len(exon_bin[t])  
                
                    edge[t] = []
                    for i in range(0,node_num_not_isolated1):
                        for j in range(0,node_num_not_isolated1):
                            if ASG[i,j] == 1:
                                temp = [0] * 2
                                temp[0] = exon_bin[t][i]
                                temp[1] = exon_bin[t][j]
                                edge[t].append(temp)


                    subpath_constraints[t].sort(key=lambda x:(x[0],x[1],x[2]))  
                    for i in range(len(subpath_constraints[t]) - 1,0,-1):
                        flag = 0
                        for j in range(i - 1,-1,-1):
                            if subpath_constraints[t][i][0] == subpath_constraints[t][j][-1]:       
                                t1 = exon_bin[t].index(subpath_constraints[t][j][-1])
                                t2 = exon_bin[t].index(subpath_constraints[t][j][-2])
                                s = 0
                                for k in range(0,t1):                     
                                    s = s + ASG[k,t1]
                                if s == 1 and ASG[t2,t1] == 1:
                                    flag = 1
                                    for l in range(1,len(subpath_constraints[t][i])):
                                        subpath_constraints[t][j].append(subpath_constraints[t][i][l])
                            elif subpath_constraints[t][i][0] == subpath_constraints[t][j][-2] and subpath_constraints[t][i][1] == subpath_constraints[t][j][-1]: 
                                t1 = exon_bin[t].index(subpath_constraints[t][j][-2])
                                t2 = exon_bin[t].index(subpath_constraints[t][j][-3])
                                s = 0
                                for k in range(0,t1):
                                    s = s + ASG[k,t1]
                                if s == 1 and ASG[t2,t1] == 1:
                                    flag = 1
                                    for l in range(2,len(subpath_constraints[t][i])):
                                        subpath_constraints[t][j].append(subpath_constraints[t][i][l])
                            elif subpath_constraints[t][i][0] == subpath_constraints[t][j][-3] and subpath_constraints[t][i][1] == subpath_constraints[t][j][-2] and subpath_constraints[t][i][2] == subpath_constraints[t][j][-1]:   
                                if len(subpath_constraints[t][j]) == 3:
                                    flag = 1
                                elif len(subpath_constraints[t][j]) == 4:
                                    t1 = exon_bin[t].index(subpath_constraints[t][j][-3])
                                    t2 = exon_bin[t].index(subpath_constraints[t][j][-4])
                                    s = 0
                                    for k in range(0,t1):                     
                                        s = s + ASG[k,t1]
                                    if s == 1 and ASG[t2,t1] == 1:
                                        flag = 1
                                        for l in range(3,len(subpath_constraints[t][i])):
                                            subpath_constraints[t][j].append(subpath_constraints[t][i][l])
                        if flag == 1:               #
                            subpath_constraints[t][i] = []

                    for i in range(len(subpath_constraints[t]) - 1,0,-1):
                        if subpath_constraints[t][i] == []:
                            subpath_constraints[t].pop(i)
                    temp = copy.deepcopy(subpath_constraints[t])
                    subpath_constraints[t] = []
                    for i in temp:
                        if i not in subpath_constraints[t]:
                            subpath_constraints[t].append(i)
                    new_edge[t] = edge[t] + subpath_constraints[t]
                    ASG_0 = copy.deepcopy(ASG)
                    exon_bin_1 = copy.deepcopy(exon_bin[t])
                    for i in range(0,len(subpath_constraints[t])):
                        for j in range(0,len(subpath_constraints[t][i]) - 1):
                            t1 = exon_bin[t].index(subpath_constraints[t][i][j])
                            t2 = exon_bin[t].index(subpath_constraints[t][i][j + 1])
                            ASG[t1,t2] = 0
                        t1 = exon_bin[t].index(subpath_constraints[t][i][0])
                        t2 = exon_bin[t].index(subpath_constraints[t][i][-1])
                        ASG[t1,t2] = 1
                           

                    
                    isolate = []
                    for i in range(1,node_num_not_isolated1 - 1):
                        if sum(ASG[:,i]) == 0 and sum(ASG[i,:]) == 0:
                            isolate.append(i)
                    isolate.reverse()
                    
                    for i in isolate:
                        exon_bin[t].pop(i)      
                    ASG = np.delete(ASG,isolate,0)
                    ASG = np.delete(ASG,isolate,1)
                    
                    node_num_not_isolated2 = len(exon_bin[t])
                    edge_num =  int(np.sum(ASG))
                    new_node_num = node_num_not_isolated2 + edge_num
            
                    while np.sum(ASG) < 4 * edge_num:   
                        temp = 0
                        for i in range(0,len(ASG)):
                            for j in range(len(ASG) - 1,-1,-1):
                                if int(ASG[i,j]) == 1:
                                    temp = 1
                                    break 
                            if temp == 1:                       
                                break

                        ASG[i,j] = 0
                        ASG = np.insert(ASG,i + 1,values = 0,axis = 0)                  
                        ASG = np.insert(ASG,i + 1,values = 0,axis = 1)
                        ASG[i,i + 1] = 2
                        ASG[i + 1,j + 1] = 2
                        exon_bin[t].insert(i + 1,[])
                    for i in range(0,new_node_num):                 
                        for j in range(0,new_node_num):
                            if ASG[i,j] == 2:
                                ASG[i,j] = 1
                

                    ASG_tc = copy.deepcopy(ASG)      

                    for k in range(0,new_node_num):
                        for i in range(0,new_node_num):
                            for j in range(0,new_node_num):
                                if ASG_tc[i,k] + ASG_tc[k,j] == 2:
                                    ASG_tc[i,j] = 1

                    M = max_match(ASG_tc)                 
                    path = []         
                    x = []
                    for i in range(0,new_node_num):
                        x.append(1)             
                    while sum(x) > 0:          
                        a = x.index(1)     
                        x[a] = 0     
                        temp = [a]
                        while sum(M[a,:]) > 0:
                            for j in range(0,new_node_num):
                                if M[a,j] == 1:
                                    break
                            x[j] = 0
                            temp.append(j)
                            a = j
                        path.append(temp)

                    for i in range(0,len(path)):  
                        a = path[i][0]
                        b = path[i][-1]
                        if exon_bin[t][a] == []:
                            for j in range(0,new_node_num):
                                if ASG[j,a] != 0:
                                    break
                            path[i].insert(0,j) 
                        if exon_bin[t][b] == []:
                            for j in range(0,new_node_num):
                                if ASG[b,j] != 0:
                                    break
                            path[i].append(j)
                    for i in range(len(path) - 1,-1,-1):
                        for j in range(len(path[i]) - 1,0,-1):    
                            a = path[i][j] 
                            b = path[i][j - 1]
                            if exon_bin[t][a] == [] and exon_bin[t][b] == []:
                                path.pop(i)
                                break

                                

                    path_exon[t] = [None] * len(path)           
                    for i in range(0,len(path)):
                        temp = []
                        path_exon[t][i] = [[]]
                        for j in range(0,len(path[i])):
                            a = path[i][j]
                            temp.append(exon_bin[t][a])             
                        path_exon[t][i][0].append(temp[0])
                        for j in range(1,len(path[i])):
                            if temp[j] != []:
                                if temp[j - 1] == []:
                                    for k in path_exon[t][i]:
                                        k.append(temp[j])
                                else:
                                    a = exon_bin_1.index(temp[j - 1])
                                    b = exon_bin_1.index(temp[j])
                                    if ASG_0[a,b] == 1:
                                        for k in path_exon[t][i]:
                                            k.append(temp[j])
                                    else:
                                        ab_path = findPath(ASG_0,[a],b)
                                        path_len = []
                                        for k in ab_path:
                                            path_len.append(len(k))
                                        min_num = min(path_len)
                                        ab_min_path = ab_path[path_len.index(min_num)]
                                        for k in path_exon[t][i]:
                                            for m in range(1,len(ab_min_path)):
                                                k.append(exon_bin_1[ab_min_path[m]])
                                                
                            else:
                                a =  temp[j - 1]
                                b =  temp[j + 1]
                                l = 0
                                temp1 = []
                                for k in range(0,len(new_edge[t])):
                                    if new_edge[t][k][0] == a and new_edge[t][k][-1] == b:
                                        l = l + 1
                                        temp2 = copy.deepcopy(new_edge[t][k])
                                        temp2.pop(0)
                                        temp2.pop(-1)
                                        temp1.append(temp2)
                                p = len(path_exon[t][i])
                                path_exon[t][i] = path_exon[t][i] * l
                                for n in range(0,l):
                                    for m in range(0,p):
                                         q = n * p + m
                                         path_exon[t][i][q] = path_exon[t][i][q] + temp1[n]
                                    
                
                    temp = copy.deepcopy(path_exon[t])
                    path_exon[t] = []
                    for i in temp:
                        for j in i:
                            path_exon[t].append(j)
                    path_head = []
                    path_tail = []
                    for k in range(0,len(path_exon[t])):
                        if path_exon[t][k][0][0] != start and path_exon[t][k][0] not in path_head:
                            path_head.append(path_exon[t][k][0])
                        if path_exon[t][k][-1][-1] != end and path_exon[t][k][-1] not in path_tail:
                            path_tail.append(path_exon[t][k][-1])           
                    #head_to_start   
                    if path_head != []:
                        head_to_start =  [None] * len(path_head)
                        for k in range(0,len(path_head)):
                            b = exon_bin_1.index(path_head[k])
                            head_to_start[k] = findPath(ASG_0,[0],b)
                            path_len = []
                            for i in range(0,len(head_to_start[k])):
                                path_len.append(len(head_to_start[k][i]))
                            min_num = min(path_len)
                            for i in range(len(head_to_start[k]) - 1,-1,-1):
                                if len(head_to_start[k][i]) > min_num:
                                    head_to_start[k].pop(i)                 
                            for i in range(0,len(head_to_start[k])):
                                for j in range(0,len(head_to_start[k][i])):
                                    head_to_start[k][i][j] = exon_bin_1[head_to_start[k][i][j]]             
                    #tail_to_end    
                    if path_tail != []:
                        tail_to_end = [None] * len(path_tail)
                        for k in range(0,len(path_tail)):
                            a = exon_bin_1.index(path_tail[k])
                            tail_to_end[k] = findPath(ASG_0,[a],len(exon_bin_1) - 1)
                            path_len = []
                            for i in range(0,len(tail_to_end[k])):
                                path_len.append(len(tail_to_end[k][i]))
                            min_num = min(path_len)
                            for i in range(len(tail_to_end[k]) - 1,-1,-1):
                                if len(tail_to_end[k][i]) > min_num:
                                    tail_to_end[k].pop(i)                   
                            for i in range(0,len(tail_to_end[k])):
                                for j in range(0,len(tail_to_end[k][i])):
                                    tail_to_end[k][i][j] = exon_bin_1[tail_to_end[k][i][j]]

                    transcript_bin[t] = []
                    for k in range(0,len(path_exon[t])):
                        if path_exon[t][k][0][0] ==  start and path_exon[t][k][-1][-1] == end: 
                            transcript_bin[t].append(path_exon[t][k])
                        elif path_exon[t][k][0][0] !=  start and path_exon[t][k][-1][-1] == end: 
                            a = path_head.index(path_exon[t][k][0])
                            path_exon[t][k].pop(0)
                            for i in range(0,len(head_to_start[a])):                        
                                transcript_bin[t].append(head_to_start[a][i] + path_exon[t][k])
                        elif path_exon[t][k][0][0] ==  start and path_exon[t][k][-1][-1] != end: 
                            a = path_tail.index(path_exon[t][k][-1])
                            path_exon[t][k].pop(-1)
                            for i in range(0,len(tail_to_end[a])):                      
                                transcript_bin[t].append(path_exon[t][k] + tail_to_end[a][i])
                        elif path_exon[t][k][0][0] !=  start and path_exon[t][k][-1][-1] != end:
                            a = path_head.index(path_exon[t][k][0])
                            b = path_tail.index(path_exon[t][k][-1])
                            if len(path_exon[t][k]) == 1:
                                path_exon[t][k].pop(0)
                            else:
                                path_exon[t][k].pop(0)
                                path_exon[t][k].pop(-1)
                            temp = []
                            p = len(head_to_start[a])
                            l = len(tail_to_end[b])
                            for i in range(0,p):                        
                                temp.append(head_to_start[a][i] + path_exon[t][k])
                            temp = temp * l
                            for n in range(0,l):
                                for m in range(0,p):
                                    q = n * p + m
                                    transcript_bin[t].append(temp[q] + tail_to_end[b][n])
                    isoform[t] = []
                    for i in transcript_bin[t]:
                        if i not in isoform[t]:
                            isoform[t].append(i)
                    isoform_num[t] = len(isoform[t])    
    isoform_num_0 = copy.deepcopy(isoform_num)
    for t in range(circ_spot_num - 1,-1,-1):
        if isoform_num[t] == 0:
            isoform_num.pop(t)
            isoform.pop(t)
    isoform_num_sum = sum(isoform_num)
    isoform_len = []
    matrix = []
    print("我运行到定量前一步的前一步了，isoform_num_sum = " + str(isoform_num_sum))
    if isoform_num_sum > 0:
        for t in range(0,len(isoform_num)):     
            for i in range(0,isoform_num[t]):
                matrix.append(isoform[t][i])     
                length = 0
                for j in range(0,len(isoform[t][i])):        
                    length = length + int(isoform[t][i][j][1]) - int(isoform[t][i][j][0]) + 1
                isoform_len.append(length)
        print("我运行到定量前一步了，len(matrix) = " + str(len(matrix)))
        FPKM = quantification(matrix,isoform_len)
        if sum(FPKM) > 0:       
            matrix_0 = copy.deepcopy(matrix)
            matrix_new = []
            FPKM_new = []
            isoform_len_new = []
            for t in range(0,len(isoform_num)):    
                matrix_temp = []
                FPKM_temp = np.zeros(isoform_num[t])
                for i in range(0,isoform_num[t]):
                    matrix_temp.append(matrix[0])
                    FPKM_temp[i] = FPKM[0]
                    matrix.pop(0)
                    FPKM = np.delete(FPKM,0)
                FPKM_temp_max = np.max(FPKM_temp)
                for i in range(0,isoform_num[t]):
                    if FPKM_temp[i] > FPKM_temp_max * 0.001:
                        matrix_new.append(matrix_temp[i])
                        FPKM_new = np.concatenate((FPKM_new,[FPKM_temp[i]]))
                        temp = matrix_0.index(matrix_temp[i])
                        isoform_len_new.append(isoform_len[temp]) 
                        
            results = open("result_table.txt","a")
            for i in range(0,len(matrix_new)):
                temp =  [0] * 2
                temp[0] = matrix_new[i][0][0]
                temp[1] = matrix_new[i][-1][-1]
                junction_reads = circ_num[circ_spot.index(temp)]
                strd = strand[circ_spot.index(temp)]
                exon_len = []
                exon_offset =[]
                for k in range(0,len(matrix_new[i])):
                    exon_len.append(int(matrix_new[i][k][1]) - int(matrix_new[i][k][0]) + 1)
                    exon_offset.append(int(matrix_new[i][k][0]) - int(matrix_new[i][0][0]))
                match_frag = round(FPKM_new[i] * sum(exon_len) * match_num_total * (10 ** (-9)))
                match_frag = int(match_frag)
                if match_frag == 0:
                    match_frag = 1
                for k in range(0,len(exon_len)):
                    exon_len[k] = str(exon_len[k])
                for k in range(0,len(exon_offset)):
                    exon_offset[k] = str(exon_offset[k])                
                results.write(chr_name + "\t" + str(matrix_new[i][0][0]) + "\t" + str(matrix_new[i][-1][-1]) + "\t" + gene_id + "\t" + str(strd) + "\t" + str(len(matrix_new[i])) + "\t" + ','.join(exon_len) + "\t" + ','.join(exon_offset) + "\t" + str(junction_reads) + "\t" + str(FPKM_new[i]) + "\t" + str(match_frag) + "\n")
            results.close()


    os.remove(tempfile_path + gene_id + ".gtf")
    os.remove(tempfile_path + gene_id + ".sam")
    os.remove(tempfile_path + gene_id + "_A.sam")
    os.remove(tempfile_path + gene_id + "_B.sam")

    for t in range(0,circ_spot_num):
        os.remove(tempfile_path + gene_id + "_" + str(t) + ".sam")  
def findExonSequence(item,circ_chr,fasta):
    chrName = item[1]
    start = int(item[2])
    end = int(item[3])
    in_gtf = open(tempfile_path + chrName + ".gtf", "r")
    exon = {}
    exon_num = 1
    sequence = ""
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        if len(array[0]) < 6:
            if array[2] == "exon" and int(array[3]) >= start and int(array[4]) <= end:
                k = circ_chr.index(array[0])
                exon_left = int(array[3])
                exon_right = int(array[4])
                info = divmod(exon_left, int(len_fasta))
                info_right = divmod(exon_right, int(len_fasta))
                if(int(array[4]) - int(array[3]) > (60 -int(info[1]))):
                    for i in range(info_right[0] - info[0] + 1):
                        if i == 0:
                            exon_sequence = fasta[k][info[0] + i][info[1]:]
                            continue
                        if i == info_right[0] - info[0]:
                            exon_sequence += fasta[k][info[0] + i][:info_right[1] + 1]
                            continue
                        exon_sequence += fasta[k][info[0] + i][:]
                if (int(array[4]) - int(array[3]) <= (60 - int(info[1]))):
                        span = int(array[4]) - int(array[3])
                        exon_sequence = fasta[k][info[0]][info[1]:info[1] + span + 1]
                if exon_num == 1:
                    exon[exon_num] = str(array[3]) + '-' + str(array[4])
                    sequence += exon_sequence
                    exon_num += 1
                elif exon_num > 1:
                    if int(array[3]) > int(exon[int(exon_num) - 1].split('-')[1]):  ##如果有多个exon，第二个exon的起点需大于第一个exon的终点
                        exon[exon_num] = str(array[3]) + '-' + str(array[4])
                        sequence += exon_sequence
                        exon_num += 1
    return sequence
def calculate_mixed(gene,gene_name, chr_no, match_num_total, head_num, inputfrag_len, tempfile_path):  ############   inputfrag_len = 2*150 + 100
    gene_id = gene
    i = gene_name.index(gene_id)
    chr_name = chr_no[i]
    bound = 3
    def max_match(A):
        dim = A.shape
        m = dim[0]
        n = dim[1]
        M = np.zeros((m,n))
        for i in range(0,m):
            for j in range(0,n):
                if A[i,j] == 1:
                    M[i,j] = 1
                    break
            if np.sum(M) > 0:
                break

        while 1:
            x = []
            for i in range(0,m):
                x.append(0)
            y = []
            for j in range(0,n):
                y.append(0)

            for i in range(0,m):
                pd = 1
                for j in range(0,n):
                    if M[i,j] == 1:
                        pd = 0
                if pd == 1:
                    x[i] = - n 

            pd = 0
            while 1:
                xi = -1
                for i in range(0,m):
                    if x[i] < 0:
                        xi = i
                        break
                if xi == -1:
                    pd = 1
                    break
                x[xi] = x[xi] * (-1)
                yy = []
                for j in range(0,n):
                    if A[xi,j] == 1 and y[j] == 0:
                        y[j] = xi
                        yy.append(j)

                if len(yy) >= 1:
                    for j in range(0,len(yy)):
                        pdd = 1
                        for i in range(0,m):
                            if M[i,yy[j]] == 1:
                                x[i] = -yy[j]
                                pdd = 0
                                break
                        if pdd == 1:
                            break
                    if pdd == 1:
                        j = yy[j]
                        P1 = []
                        P2 = []
                        while 1:
                            P2.append(j)
                            P1.append(y[j])
                            j = abs(x[y[j]])
                            if j == n:
                                break
                        for i in range(0,len(P1)):
                            if M[P1[i],P2[i]] == 1:
                                M[P1[i],P2[i]] = 0
                            else:
                                M[P1[i],P2[i]] = 1
                        break
            if pd == 1:
                break
        return M
    def findPath(Graph,partialPath,destination):  
        pathLength = len(partialPath)
        lastNode = partialPath[-1]
        possiablePaths = []
        nextNodes = []
        for i in range(lastNode,destination + 1):
            if Graph[lastNode,i] == 1:
                nextNodes.append(i)
        for i in range(0,len(nextNodes)):
            if destination == nextNodes[i]:
                tmpPath = partialPath + [destination]
                possiablePaths.append(tmpPath)
                nextNodes[i] = 0
        for i in range(len(nextNodes) - 1,-1,-1):
            if nextNodes[i] == 0:
                nextNodes.pop(i)
        for i in range(0,len(nextNodes)):
            tmpPath = partialPath + [nextNodes[i]]
            tmpPsbPaths = findPath(Graph, tmpPath, destination)     
            possiablePaths = possiablePaths + tmpPsbPaths
        return possiablePaths
    def selection_mixed(matrix_mixed, isoform_len_mixed,circAST_spot_num,liner_spot_num):  #####################################matrix是[i][j][0],[i][j][1]为第i个isoform，第j个exon，0是exon起点，1是exon终点
        k = 0
        reads_name = []
        in_sam = open(tempfile_path + gene_id + "_A.sam", "r")
        while k < head_num:
            line1 = in_sam.readline()  ##########################用while要readline一下跳过前面的行数
            k += 1
        k = 0
        while k != sam_row_num_A:
            # print(k)
            line1 = in_sam.readline()
            line1 = line1.strip('\n')
            array1 = line1.split()
            data11 = re.findall("\d+", array1[5])
            data12 = re.findall("\D", array1[5])
            line2 = in_sam.readline()
            line2 = line2.strip('\n')
            array2 = line2.split()
            data21 = re.findall("\d+", array2[5])
            data22 = re.findall("\D", array2[5])
            count_array = np.zeros(len(matrix_mixed))
            frag_len_isoform = np.zeros(len(matrix_mixed))
            if data12.count("N") == 0 and data22.count("N") == 0:
                match1 = 0
                for i in range(0, len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match1 = match1 + int(data11[i])  #######match1记录左端match上的长度
                match2 = 0
                for i in range(0, len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match2 = match2 + int(data21[i])  #######match2记录右端match上的长度
                t11 = int(array1[3])
                t12 = t11 + match1 - 1
                t21 = int(array2[3])
                t22 = t21 + match2 - 1
                for i in range(0, len(matrix_mixed)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0, len(matrix_mixed[i])):
                        if t11 >= int(matrix_mixed[i][j][0]) - 2 and t12 <= int(matrix_mixed[i][j][1]) + 2:  #######左端完全落在某一个exon里
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix_mixed[i][j][1]) - t11 + 1  ########################没看懂为啥是exon终点-左端reads长度
                            break
                    if t > -1:
                        if t21 >= int(matrix_mixed[i][t][0]) - 2 and t22 <= int(matrix_mixed[i][t][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = t21 - t11 + match2  ###########################frag_len_isoform[i]应该是在第i个isoform上fragment的长度
                        else:
                            for j in range(t + 1, len(matrix_mixed[i])):
                                if t21 >= int(matrix_mixed[i][j][0]) - 2 and t22 <= int(matrix_mixed[i][j][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + t22 - int(matrix_mixed[i][j][0]) + 1  ###########明白了
                                    break
                                else:
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - int(matrix_mixed[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000  #######如果友谊端匹配不上，设一个极大值
            elif data12.count("N") == 0 and data22.count("N") == 1:
                match1 = 0
                for i in range(0, len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match1 = match1 + int(data11[i])
                t11 = int(array1[3])
                t12 = t11 + match1 - 1
                N_id2 = data22.index("N")  ##########如果有多个N出现的情况呢？？？？
                match21 = 0
                for i in range(0, N_id2):
                    if data22[i] == "M" or data22[i] == "D":
                        match21 = match21 + int(data21[i])
                match22 = 0
                for i in range(N_id2 + 1, len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match22 = match22 + int(data21[i])
                t21 = int(array2[3])
                t22 = t21 + match21 - 1
                t23 = t22 + int(data21[N_id2]) + 1
                t24 = t23 + match22 - 1
                for i in range(0, len(matrix_mixed)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0, len(matrix_mixed[i]) - 1):
                        if t11 >= int(matrix_mixed[i][j][0]) - 2 and t12 <= int(matrix_mixed[i][j][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix_mixed[i][j][1]) - t11 + 1
                            break
                    if t > -1:
                        if t21 >= int(matrix_mixed[i][t][0]) - 2 and abs(t22 - int(matrix_mixed[i][t][1])) <= 2 and abs(t23 - int(matrix_mixed[i][t + 1][0])) <= 2 and t24 <= int(matrix_mixed[i][t + 1][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = int(matrix_mixed[i][t][1]) - t11 + 1 + match22
                        else:
                            for j in range(t + 1, len(matrix_mixed[i]) - 1):
                                if t21 >= int(matrix_mixed[i][j][0]) - 2 and abs(t22 - int(matrix_mixed[i][j][1])) <= 2 and abs(t23 - int(matrix_mixed[i][j + 1][0])) <= 2 and t24 <= int(matrix_mixed[i][j + 1][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - int(matrix_mixed[i][j][0]) + 1 + match22
                                    break
                                else:
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - int(matrix_mixed[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000

            elif data12.count("N") == 1 and data22.count("N") == 0:
                N_id1 = data12.index("N")
                match11 = 0
                for i in range(0, N_id1):
                    if data12[i] == "M" or data12[i] == "D":
                        match11 = match11 + int(data11[i])
                match12 = 0
                for i in range(N_id1 + 1, len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match12 = match12 + int(data11[i])
                t11 = int(array1[3])
                t12 = t11 + match11 - 1
                t13 = t12 + int(data11[N_id1]) + 1
                t14 = t13 + match12 - 1
                match2 = 0
                for i in range(0, len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match2 = match2 + int(data21[i])
                t21 = int(array2[3])
                t22 = t21 + match2 - 1

                for i in range(0, len(matrix_mixed)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0, len(matrix_mixed[i]) - 1):
                        if t11 >= int(matrix_mixed[i][j][0]) - 2 and abs(t12 - int(matrix_mixed[i][j][1])) <= 2 and abs(t13 - int(matrix_mixed[i][j + 1][0])) <= 2 and t14 <= int(matrix_mixed[i][j + 1][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix_mixed[i][j][1]) - t11 + 1
                            break
                    if t > -1:
                        for j in range(t + 1, len(matrix_mixed[i])):
                            if t21 >= int(matrix_mixed[i][j][0]) - 2 and t22 <= int(matrix_mixed[i][j][1]) + 2:
                                match[1] = 1
                                frag_len_isoform[i] = frag_len_isoform[i] + t21 - int(matrix_mixed[i][j][0]) + match2
                                break
                            else:
                                frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - int(matrix_mixed[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000

            elif data12.count("N") == 1 and data22.count("N") == 1:
                N_id1 = data12.index("N")
                match11 = 0
                for i in range(0, N_id1):
                    if data12[i] == "M" or data12[i] == "D":
                        match11 = match11 + int(data11[i])
                match12 = 0
                for i in range(N_id1 + 1, len(data12)):
                    if data12[i] == "M" or data12[i] == "D":
                        match12 = match12 + int(data11[i])
                t11 = int(array1[3])
                t12 = t11 + match11 - 1
                t13 = t12 + int(data11[N_id1]) + 1
                t14 = t13 + match12 - 1
                N_id2 = data22.index("N")
                match21 = 0
                for i in range(0, N_id2):
                    if data22[i] == "M" or data22[i] == "D":
                        match21 = match21 + int(data21[i])
                match22 = 0
                for i in range(N_id2 + 1, len(data22)):
                    if data22[i] == "M" or data22[i] == "D":
                        match22 = match22 + int(data21[i])
                t21 = int(array2[3])
                t22 = t21 + match21 - 1
                t23 = t22 + int(data21[N_id2]) + 1
                t24 = t23 + match22 - 1

                for i in range(0, len(matrix_mixed)):
                    match = np.zeros(2)
                    t = -1
                    for j in range(0, len(matrix_mixed[i]) - 2):
                        if t11 >= int(matrix_mixed[i][j][0]) - 2 and abs(t12 - int(matrix_mixed[i][j][1])) <= 2 and abs(t13 - int(matrix_mixed[i][j + 1][0])) <= 2 and t14 <= int(matrix_mixed[i][j + 1][1]) + 2:
                            t = j
                            match[0] = 1
                            frag_len_isoform[i] = int(matrix_mixed[i][j][1]) - t11 + 1 + int(matrix_mixed[i][j + 1][1]) - int(matrix_mixed[i][j + 1][0]) + 1
                            break
                    if t > -1 and t + 1 < len(matrix_mixed[i]):
                        t = t + 1
                        if t21 >= int(matrix_mixed[i][t - 1][0]) - 2 and abs(t22 - int(matrix_mixed[i][t - 1][1])) <= 2 and abs(t23 - int(matrix_mixed[i][t][0])) <= 2 and t24 <= int(matrix_mixed[i][t][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = int(matrix_mixed[i][t - 1][1]) - t11 + 1 + match22
                        elif t21 >= int(matrix_mixed[i][t][0]) - 2 and abs(t22 - int(matrix_mixed[i][t][1])) <= 2 and abs(t23 - int(matrix_mixed[i][t + 1][0])) <= 2 and t24 <= int(matrix_mixed[i][t + 1][1]) + 2:
                            match[1] = 1
                            frag_len_isoform[i] = frag_len_isoform[i] + match22
                        else:
                            for j in range(t + 1, len(matrix_mixed[i]) - 1):
                                if t21 >= int(matrix_mixed[i][j][0]) and abs(t22 - int(matrix_mixed[i][j][1])) <= 2 and abs(t23 - int(matrix_mixed[i][j + 1][0])) <= 2 and t24 <= int(matrix_mixed[i][j + 1][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - int(matrix_mixed[i][j][0]) + match22
                                    break
                                else:
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - int(matrix_mixed[i][j][0]) + 1
                    count_array[i] = match[0] * match[1]
                    if count_array[i] == 0:
                        frag_len_isoform[i] = 1000000
            else:
                for i in range(0, len(matrix_mixed)):  ############为何其他情况长度都记为1000000？？中间跨越多个N的情况呢？
                    frag_len_isoform[i] = 1000000
            if sum(count_array) > 0:
                match_count = open(tempfile_path + gene_id + "_match.txt", "a")  # a属性，每次打开后接着写，所以用完应该删除
                for i in range(0, len(matrix_mixed)):
                    match_count.write(str(count_array[i]) + "\t")
                match_count.write("\n")
                match_count.close()
                fraglen = open(tempfile_path + gene_id + "_fraglen.txt", "a")
                for i in range(0, len(matrix_mixed)):
                    fraglen.write(str(frag_len_isoform[i]) + "\t")
                fraglen.write("\n")
                fraglen.close()
                reads_name.append(array1[0])
            k = k + 2
        in_sam.close()
        ###################################################################################
        if sam_row_num_B > 0:
            k = 0
            in_sam = open(tempfile_path + gene_id + "_B.sam", "r")
            while k < head_num:
                line1 = in_sam.readline()
                k += 1
            k = 0
            while k != sam_row_num_B:
                line1 = in_sam.readline()
                line1 = line1.strip('\n')
                array1 = line1.split()
                data11 = re.findall("\d+", array1[5])
                data12 = re.findall("\D", array1[5])
                line2 = in_sam.readline()
                line2 = line2.strip('\n')
                array2 = line2.split()
                data21 = re.findall("\d+", array2[5])
                data22 = re.findall("\D", array2[5])
                count_array = np.zeros(len(matrix_mixed))
                frag_len_isoform = np.zeros(len(matrix_mixed))
                if data12.count("N") == 0 and data22.count("N") == 0:
                    match1 = 0
                    for i in range(0, len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match1 = match1 + int(data11[i])
                    match2 = 0
                    for i in range(0, len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match2 = match2 + int(data21[i])
                    t11 = int(array1[3])
                    t12 = t11 + match1 - 1
                    t21 = int(array2[3])
                    t22 = t21 + match2 - 1
                    for i in range(0, circAST_spot_num):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0, len(matrix_mixed[i])):
                            if t11 >= int(matrix_mixed[i][j][0]) - 2 and t12 <= int(matrix_mixed[i][j][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix_mixed[i][j][0])
                                for m in range(0, t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                break
                        if t > -1:
                            for j in range(t, len(matrix_mixed[i])):
                                if t21 >= int(matrix_mixed[i][j][0]) - 2 and t22 <= int(matrix_mixed[i][j][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - t22
                                    for m in range(j + 1, len(matrix_mixed[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match1 + match2
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000
                    for p in range(circAST_spot_num, len(matrix_mixed)):
                        frag_len_isoform[p] = 1000000

                elif data12.count("N") == 0 and data22.count("N") == 1:
                    match1 = 0
                    for i in range(0, len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match1 = match1 + int(data11[i])
                    t11 = int(array1[3])
                    t12 = t11 + match1 - 1
                    N_id2 = data22.index("N")
                    match21 = 0
                    for i in range(0, N_id2):
                        if data22[i] == "M" or data22[i] == "D":
                            match21 = match21 + int(data21[i])
                    match22 = 0
                    for i in range(N_id2 + 1, len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match22 = match22 + int(data21[i])
                    t21 = int(array2[3])
                    t22 = t21 + match21 - 1
                    t23 = t22 + int(data21[N_id2]) + 1
                    t24 = t23 + match22 - 1
                    for i in range(0, circAST_spot_num):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0, len(matrix_mixed[i]) - 1):
                            if t11 >= int(matrix_mixed[i][j][0]) - 2 and t12 <= int(matrix_mixed[i][j][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix_mixed[i][j][0])
                                for m in range(0, t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                break
                        if t > -1:
                            for j in range(t, len(matrix_mixed[i]) - 1):
                                if t21 >= int(matrix_mixed[i][j][0]) - 2 and abs(t22 - int(matrix_mixed[i][j][1])) <= 2 and abs(t23 - int(matrix_mixed[i][j + 1][0])) <= 2 and t24 <= int(matrix_mixed[i][j + 1][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j + 1][1]) - t24
                                    for m in range(j + 2, len(matrix_mixed[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match1 + match21 + match22
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000
                    for p in range(circAST_spot_num, len(matrix_mixed)):
                        frag_len_isoform[p] = 1000000

                elif data12.count("N") == 1 and data22.count("N") == 0:
                    N_id1 = data12.index("N")
                    match11 = 0
                    for i in range(0, N_id1):
                        if data12[i] == "M" or data12[i] == "D":
                            match11 = match11 + int(data11[i])
                    match12 = 0
                    for i in range(N_id1 + 1, len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match12 = match12 + int(data11[i])
                    t11 = int(array1[3])
                    t12 = t11 + match11 - 1
                    t13 = t12 + int(data11[N_id1]) + 1
                    t14 = t13 + match12 - 1
                    match2 = 0
                    for i in range(0, len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match2 = match2 + int(data21[i])
                    t21 = int(array2[3])
                    t22 = t21 + match2 - 1

                    for i in range(0, circAST_spot_num):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0, len(matrix_mixed[i]) - 1):
                            if t11 >= int(matrix_mixed[i][j][0]) - 2 and abs(t12 - int(matrix_mixed[i][j][1])) <= 2 and abs(t13 - int(matrix_mixed[i][j + 1][0])) <= 2 and t14 <= int(matrix_mixed[i][j + 1][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix_mixed[i][j][0])
                                for m in range(0, t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                break
                        if t > -1:
                            for j in range(t + 1, len(matrix_mixed[i])):
                                if t21 >= int(matrix_mixed[i][j][0]) - 2 and t22 <= int(matrix_mixed[i][j][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j][1]) - t22
                                    for m in range(j + 1, len(matrix_mixed[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match11 + match12 + match2
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000
                    for p in range(circAST_spot_num, len(matrix_mixed)):
                        frag_len_isoform[p] = 1000000

                elif data12.count("N") == 1 and data22.count("N") == 1:
                    N_id1 = data12.index("N")
                    match11 = 0
                    for i in range(0, N_id1):
                        if data12[i] == "M" or data12[i] == "D":
                            match11 = match11 + int(data11[i])
                    match12 = 0
                    for i in range(N_id1 + 1, len(data12)):
                        if data12[i] == "M" or data12[i] == "D":
                            match12 = match12 + int(data11[i])
                    t11 = int(array1[3])
                    t12 = t11 + match11 - 1
                    t13 = t12 + int(data11[N_id1]) + 1
                    t14 = t13 + match12 - 1
                    N_id2 = data22.index("N")
                    match21 = 0
                    for i in range(0, N_id2):
                        if data22[i] == "M" or data22[i] == "D":
                            match21 = match21 + int(data21[i])
                    match22 = 0
                    for i in range(N_id2 + 1, len(data22)):
                        if data22[i] == "M" or data22[i] == "D":
                            match22 = match22 + int(data21[i])
                    t21 = int(array2[3])
                    t22 = t21 + match21 - 1
                    t23 = t22 + int(data21[N_id2]) + 1
                    t24 = t23 + match22 - 1

                    for i in range(0, circAST_spot_num):
                        match = np.zeros(2)
                        t = -1
                        for j in range(0, len(matrix_mixed[i]) - 2):
                            if t11 >= int(matrix_mixed[i][j][0]) - 2 and abs(t12 - int(matrix_mixed[i][j][1])) <= 2 and abs(t13 - int(matrix_mixed[i][j + 1][0])) <= 2 and t14 <= int(matrix_mixed[i][j + 1][1]) + 2:
                                t = j
                                match[0] = 1
                                frag_len_isoform[i] = t11 - int(matrix_mixed[i][j][0])
                                for m in range(0, t):
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                break
                        if t > -1 and t + 1 < len(matrix_mixed[i]):
                            for j in range(t, len(matrix_mixed[i]) - 1):
                                if t21 >= int(matrix_mixed[i][j][0]) and abs(t22 - int(matrix_mixed[i][j][1])) <= 2 and abs(t23 - int(matrix_mixed[i][j + 1][0])) <= 2 and t24 <= int(matrix_mixed[i][j + 1][1]) + 2:
                                    match[1] = 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][j + 1][1]) - t24
                                    for m in range(j + 2, len(matrix_mixed[i])):
                                        frag_len_isoform[i] = frag_len_isoform[i] + int(matrix_mixed[i][m][1]) - int(matrix_mixed[i][m][0]) + 1
                                    frag_len_isoform[i] = frag_len_isoform[i] + match11 + match12 + match21 + match22
                                    break
                        count_array[i] = match[0] * match[1]
                        if count_array[i] == 0:
                            frag_len_isoform[i] = 1000000
                    for p in range(circAST_spot_num, len(matrix_mixed)):
                        frag_len_isoform[p] = 1000000
                else:
                    for i in range(0, len(matrix_mixed)):
                        frag_len_isoform[i] = 1000000
                if sum(count_array) > 0:
                    match_count = open(tempfile_path + gene_id + "_match.txt", "a")
                    for i in range(0, len(matrix_mixed)):
                        match_count.write(str(count_array[i]) + "\t")
                    match_count.write("\n")
                    match_count.close()
                    fraglen = open(tempfile_path + gene_id + "_fraglen.txt", "a")
                    for i in range(0, len(matrix_mixed)):
                        fraglen.write(str(frag_len_isoform[i]) + "\t")
                    fraglen.write("\n")
                    fraglen.close()
                    reads_name.append(array1[0])
                k = k + 2
            in_sam.close()

        FPKM = np.zeros(len(matrix_mixed))
        # selected_isoform_circ = []    ####修改得到筛选isoform的形式
        # selected_isoform_liner = []
        circ_reads = []
        liner_reads = []
        selected_isoform_circ = []
        selected_isoform_liner = []
        if os.path.exists(tempfile_path + gene_id + "_fraglen.txt"):
            y = np.loadtxt(tempfile_path + gene_id + "_fraglen.txt")  #
            dim = np.shape(y)
            if dim == ():  ##只有一个fraglen值
                os.remove(tempfile_path + gene_id + "_match.txt")
                os.remove(tempfile_path + gene_id + "_fraglen.txt")
                return 0,0,0,0, circAST_spot_num, liner_spot_num
            if dim != ():
                # label = np.zeros((dim[0], 2))
                if len(dim) == 1:  #############################numpy数组的形状是一维的
                    if len(matrix_mixed) == 1:  #######################matrix中只有一个isoform
                        if circAST_spot_num == 0:
                            for j in range(0, dim[0]):
                                liner_reads.append(reads_name[j])
                            if(isoform_len_mixed[0] > inputfrag_len):
                                FPKM[0] = 10 ** 9 * (dim[0]) / ((isoform_len_mixed[0] - inputfrag_len) * match_num_total * 1.0)
                                FPKM[0] = float('%.4f' % FPKM[0])
                            else:
                                FPKM[0] = 10 ** 9 * (dim[0]) / ((isoform_len_mixed[0]) * match_num_total * 1.0)
                                FPKM[0] = float('%.4f' % FPKM[0])
                        else:
                            for j in range(0, dim[0]):
                                circ_reads.append(reads_name[j])
                            for k in range(0, ciri_spot_num):
                                if circ_spot[k][0] == matrix_mixed[0][0][0] and circ_spot[k][1] == matrix_mixed[0][-1][-1]:
                                    FPKM[0] = 10 ** 9 * (circ_num[k] + dim[0]) / ((isoform_len_mixed[0] - bound * 2 + 2.0) * match_num_total * 1.0)
                                    FPKM[0] = float('%.4f' % FPKM[0])
                    if len(matrix_mixed) != 1:  #######################matrix中只有一个reads
                        os.remove(tempfile_path + gene_id + "_match.txt")
                        os.remove(tempfile_path + gene_id + "_fraglen.txt")
                        return 0,0,0,0, circAST_spot_num, liner_spot_num

                elif len(dim) == 2:
                    match_num_y = len(y)
                    gaussian_dist = np.random.normal(loc=inputfrag_len,scale=100,size=1000000)
                    gaussian_dist = np.around(gaussian_dist, 0)
                    for i in range(0, ciri_spot_num):
                        count_array = np.zeros(len(matrix_mixed))
                        frag_len_isoform = np.zeros(len(matrix_mixed))
                        for k in range(0, circAST_spot_num):
                            if circ_spot[i][0] == matrix_mixed[k][0][0] and circ_spot[i][1] == matrix_mixed[k][-1][-1]:
                                count_array[k] = 1
                                frag_len_isoform[k] = inputfrag_len
                            else:
                                frag_len_isoform[k] = 1000000
                        for p in range(circAST_spot_num, len(matrix_mixed)):
                            frag_len_isoform[p] = 1000000
                        if sum(count_array) > 0:
                            match_count = open(tempfile_path + gene_id + "_match.txt", "a")
                            for j in range(0, circ_num[i]):
                                for t in range(0, len(matrix_mixed)):
                                    match_count.write(str(count_array[t]) + "\t")
                                match_count.write("\n")
                            match_count.close()
                            fraglen = open(tempfile_path + gene_id + "_fraglen.txt", "a")
                            for j in range(0, circ_num[i]):
                                for t in range(0, len(matrix_mixed)):
                                    fraglen.write(str(frag_len_isoform[t]) + "\t")
                                fraglen.write("\n")
                            fraglen.close()
                    new_frag = np.loadtxt(tempfile_path + gene_id + "_fraglen.txt")
                    match_num = len(new_frag)
                    dim_new = np.shape(new_frag)
                    label = np.zeros((dim_new[0], 2))

                    lt = np.array(isoform_len_mixed)
                    lu = np.zeros(len(matrix_mixed))
                    lu_iter = np.zeros(len(matrix_mixed))
                    error = np.zeros(len(matrix_mixed))
                    F = np.zeros((match_num, len(matrix_mixed)))  #########################F与new_frag文件维度相同，存放高斯概率
                    for j in range(0, len(matrix_mixed)):
                        lu[j] = 1.0 / len(matrix_mixed)  ################lu值不变  原来：lu[j] = 1.0 / len(matrix_mixed) 修改：
                    i = 0
                    while (i < len(F)):
                        F_sum = 0
                        for j in range(0, len(matrix_mixed)):
                            if new_frag[i, j] == 1000000:
                                F[i, j] = 0
                            else:
                                F[i, j] = len(gaussian_dist[np.where(gaussian_dist == new_frag[i, j])]) / 1000000  #####这步是在计算泊松分布概率,可以修改为（左-1 + 中 + 右+1）/3000000
                            F_sum += F[i, j]
                        if F_sum == 0:
                            F = np.delete(F, i, axis=0)
                            new_frag = np.delete(new_frag, i, axis=0)
                            match_num  = match_num - 1
                            i -= 1
                        i += 1
                    dim_new = np.shape(new_frag)
                    label = np.zeros((dim_new[0], 2))
                    flag = 1
                    flag1 = 0
                    z = np.zeros((match_num, len(matrix_mixed)))
                    iterate = 0
                    lul = np.zeros(len(matrix_mixed))
                    alpha = np.zeros(len(matrix_mixed))
                    alpha_iter = np.zeros(len(matrix_mixed))
                    al_lt = np.zeros(len(matrix_mixed))
                    z_fenzi = np.zeros((match_num, len(matrix_mixed)))
                    while flag == 1:
                        t = 0
                        for j in range(0, len(matrix_mixed)):
                            if j < circAST_spot_num:
                                lul[j] = lu[j] * (lt[j] - bound * 2 + 2.0)  ####lt是每一个isoform的长度  * lu[j] =  1.0 / len(matrix)   均一化
                            if j >= circAST_spot_num:
                                if(lt[j] > inputfrag_len):
                                    lul[j] = lu[j] * (lt[j] - inputfrag_len)
                                else:
                                    lul[j] = lu[j] * (lt[j])
                        t = sum(lul)  #################   t是isoform的均值
                        for j in range(0, len(matrix_mixed)):
                            alpha[j] = lul[j] / (t * 1.0)  ##############每一个lul占均值的比值，每一个isoform长度占所有isoform长度总和的比例
                        for i in range(0, match_num):
                            sum1 = 0
                            distribution = []
                            for j in range(0, len(matrix_mixed)):
                                if j < circAST_spot_num:
                                    z_fenzi[i, j] = alpha[j] * F[i, j] / (lt[j] - bound * 2 + 2.0)  ############
                                if j >= circAST_spot_num:
                                    if(lt[j] > new_frag[i, j]):
                                        z_fenzi[i, j] = alpha[j] * F[i, j] / (lt[j] - new_frag[i, j])  ############bound是什么？bound = 3？
                                    else:
                                        z_fenzi[i, j] = alpha[j] * F[i, j] / (lt[j])
                                sum1 = sum1 + z_fenzi[i, j]
                            if sum1 != 0:
                                for j in range(0, len(matrix_mixed)):  ####################################对一行进行操作，即对每一个fragment，计算他们在不同isoform的概率
                                    z[i, j] = z_fenzi[i, j] * 1.0 / sum1  ################这里的z[i,j]是后验概率
                            else:
                                flag = 0  #################什么时候才会跳出while循环？sum1怎么为0       如果这一个fragment没有比对到任何isoform
                                break
                        if flag == 1:
                            n = sum(z)
                            for j in range(0, len(matrix_mixed)):
                                # alpha_iter[j] =  post_count[j] * 1.0 / match_num  ###平均z[,j]
                                alpha_iter[j] = n[j] * 1.0 / match_num  ##这样来获得概率值误差更小？
                            for j in range(0, len(matrix_mixed)):
                                if j < circAST_spot_num:
                                    al_lt[j] = alpha_iter[j] * 1.0 / (lt[j] - bound * 2 + 2.0)  ####平均z[,j]对第j个isoform的长度取平均
                                if j >= circAST_spot_num:
                                    if(lt[j] > inputfrag_len):
                                        al_lt[j] = alpha_iter[j] * 1.0 / (lt[j] - inputfrag_len)  ####平均z[,j]对第j个isoform的长度取平均
                                    else:
                                        al_lt[j] = alpha_iter[j] * 1.0 / (lt[j])  ####平均z[,j]对第j个isoform的长度取平均
                            k = sum(al_lt)  ########################### k是
                            for j in range(0, len(matrix_mixed)):
                                lu_iter[j] = al_lt[j] * 1.0 / k  #########################计算后每个isoform的权重
                                error[j] = abs(lu[j] - lu_iter[j])
                            iterate = iterate + 1
                            # print(error)
                            if sum(error) < 0.001:
                                flag = 0
                                flag1 = 1
                            else:
                                for j in range(0, len(matrix_mixed)):
                                    lu[j] = lu_iter[j]
                    if flag1 == 1:
                        total_count = np.zeros(len(matrix_mixed))
                        # for i in range(match_num):
                            # if label[i][0] == 0:
                                # circ_reads.append(reads_name[i])
                            # if label[i][0] == 1:
                                # liner_reads.append(reads_name[i])
                        for j in range(0, len(matrix_mixed)):
                            total_count[j] = alpha_iter[j] * match_num
                            if j < circAST_spot_num:
                                selected_isoform_circ.append(matrix_mixed[j])
                                FPKM[j] = 10 ** 9 * alpha_iter[j] * match_num / (int(lt[j] - bound * 2 + 2.0) * match_num_total * 1.0)
                                FPKM[j] = float('%.4f' % FPKM[j])
                            if j >= circAST_spot_num:
                                selected_isoform_liner.append(matrix_mixed[j])
                                if(lt[j] > inputfrag_len):
                                    FPKM[j] = 10 ** 9 * alpha_iter[j] * match_num / (int(lt[j] - inputfrag_len) * match_num_total * 1.0)
                                else:
                                    FPKM[j] = 10 ** 9 * alpha_iter[j] * match_num / (int(lt[j]) * match_num_total * 1.0)
                                FPKM[j] = float('%.4f' % FPKM[j])
                        post_count_ratio = np.zeros(len(matrix_mixed))
                        post_count_sum_circ = 0
                        post_count_sum_liner = 0
                        recount = 0
                        for j in range(0, len(matrix_mixed)):  ##希望能够筛掉一些假阳性reads
                            if j < circAST_spot_num:
                                post_count_sum_circ += total_count[j]
                            if j >= circAST_spot_num:
                                post_count_sum_liner += total_count[j]
                        for j in range(0, len(matrix_mixed)):  ##希望能够筛掉一些假阳性reads
                            if j < circAST_spot_num:
                                post_count_ratio[j] = total_count[j] / post_count_sum_circ
                            if j >= circAST_spot_num:
                                post_count_ratio[j] = total_count[j] / post_count_sum_liner
                        p = 0
                        post_count_ratio = list(post_count_ratio)
                        # print("post_count_ratio_before = " + str(post_count_ratio) + "matrix_mixed_before = " + str(matrix_mixed))
                        while (p < len(post_count_ratio)):
                            if float(post_count_ratio[p]) < 0.04:
                                post_count_ratio.pop(p)
                                matrix_mixed.pop(p)
                                isoform_len_mixed.pop(p)
                                if p < circAST_spot_num:
                                    circAST_spot_num -= 1
                                if p >= circAST_spot_num:
                                    liner_spot_num -= 1
                                recount = 1
                                p -= 1
                            p += 1
                        # print("post_count_ratio_after = " + str(post_count_ratio) + "matrix_mixed_after = " + str(matrix_mixed))
                        if recount == 1:  #########################################想要过滤掉假阳性的环形isoform
                            os.remove(tempfile_path + gene_id + "_match.txt")
                            os.remove(tempfile_path + gene_id + "_fraglen.txt")
                            return circ_reads, FPKM, matrix_mixed, isoform_len_mixed, circAST_spot_num, liner_spot_num
            print("gene_id = " + str(gene_id) + "\t" + "circAST_spot_num = " + str(circAST_spot_num) + "\t" + "liner_spot_num = " + str(liner_spot_num) + "\t" + "matrix_mixed = " + str(matrix_mixed) + "\n")
            os.remove(tempfile_path + gene_id + "_match.txt")
            os.remove(tempfile_path + gene_id + "_fraglen.txt")
            # out_sam = open(tempfile_path + "split_record", "a")
            # out_sam.write("gene_id = " + str(gene_id) + "\t" + "circ_reads_len = " + str(len(circ_reads)) + "\t" + "liner_reads_len = " + str(len(liner_reads)) + "\t" + "total_reads_len = " + str(sam_row_num_A) + "\t" + "reads_name = " + str(len(reads_name)) + "\t" + "isoform_len = " + str(isoform_len_mixed) + "\t" + "circAST_spot_num = " + str(circAST_spot_num) + "\t" + "liner_spot_num = " + str(len(liner_spot)) + "\t" + "post_count" + str(post_count) + "\t" + "lu" + str(lu) + "\n")
            # out_sam.close()
            # print("selected_isoform_circ =  ",selected_isoform_circ,"selected_isoform_liner =  ",selected_isoform_liner,"post_count = ",post_count)
        return circ_reads, FPKM, matrix_mixed, isoform_len_mixed, circAST_spot_num, liner_spot_num

    #######在gtf文件中找属于该染色体，该基因位点的exon条目，写入gene_id.gtf，记录基因的起始和中止坐标
    gene_coordinate_min = 0
    gene_coordinate_max = 0
    row = 0
    in_gtf = open(tempfile_path + chr_name + ".gtf")
    out_gtf = open(tempfile_path + gene_id + ".gtf", "w")
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        t = array.index("gene_id")
        temp_name = array[t + 1][1:-2]
        if array[0] == chr_name and array[2] == "exon" and temp_name == gene_id:
            out_gtf.write(line + "\n")
            row += 1
            if row == 1:
                gene_coordinate_min = int(array[3])
            if int(array[4]) > gene_coordinate_max:
                gene_coordinate_max = int(array[4])
    in_gtf.close()
    out_gtf.close()
    ##处理selected_circAST，先组装
    sam_row_num_A = 0
    i = 0
    in_sam = open(tempfile_path + "selected_" + chr_name + "_A.sam","r")
    out_sam = open(tempfile_path + gene_id + "_A.sam","a")
    for line in in_sam:
        line = line.strip('\n')
        if i < head_num:
            out_sam.write(line + "\n")
        else:
            k = i - head_num
            line = line.strip('\n')
            array = line.split()
            if k % 2 == 1:
                t1 = int(array[7])
                data = re.findall("\d+\D",array[5])
                s = 0
                for j in data:
                    if j[-1] in CIGAR_set1:
                        s = s + int(j[0:-1])
                t2 = int(array[3]) + s - 1
                if t1 >= gene_coordinate_min - 2 and t2 <= gene_coordinate_max + 2:
                    out_sam.write(line1 + "\n")
                    out_sam.write(line + "\n")
                    sam_row_num_A += 2
            line1 = line
        i = i + 1
    in_sam.close()
    out_sam.close()

    shutil.copyfile(tempfile_path + gene_id + "_A.sam", tempfile_path + gene_id + ".sam")               
    sam_row_num_B = 0
    i = 0
    in_sam = open(tempfile_path + "selected_" + chr_name + "_B.sam","r") 
    out_sam1 = open(tempfile_path + gene_id + "_B.sam","w")
    out_sam2 = open(tempfile_path + gene_id + ".sam","a")
    for line in in_sam:
        line = line.strip('\n')
        if i < head_num:
            out_sam1.write(line + "\n")
        else:
            k = i - head_num
            line = line.strip('\n')
            array = line.split()
            if k % 2 == 1:
                t1 = int(array[7])
                data = re.findall("\d+\D",array[5])
                s = 0
                for j in data:
                    if j[-1] in CIGAR_set1:
                        s = s + int(j[0:-1])
                t2 = int(array[3]) + s - 1
                if t1 >= gene_coordinate_min - 2 and t2 <= gene_coordinate_max + 2:
                    out_sam1.write(line1 + "\n")
                    out_sam1.write(line + "\n")
                    sam_row_num_B += 2
                    out_sam2.write(line1 + "\n")
                    out_sam2.write(line + "\n")
            line1 = line
        i = i + 1
    in_sam.close()
    out_sam1.close()
    out_sam2.close()
    
    ##组装
    circ_spot = []
    circ_num = []
    strand = []
    in_txt = open(tempfile_path + "circRNA_list_selected.txt","r")
    for line in in_txt:
        line = line.strip('\n')
        array = line.split()
        if array[0] == chr_name and array[4] == gene_id:
            temp =  [0] * 2
            temp[0] = int(array[1]) + 1 
            temp[1] = int(array[2])
            circ_spot.append(temp)
            circ_num.append(int(float(array[6])))
            strand.append(array[3])
    in_txt.close()
    
    circ_spot_num = len(circ_spot)
    sam_num = [0] * circ_spot_num     
    for t in range(0,circ_spot_num):
        i = 0
        sam_num[t] = 0
        in_sam = open(tempfile_path + gene_id + ".sam","r")
        out_sam = open(tempfile_path + gene_id + "_" + str(t) + ".sam","w")   
        for line in in_sam:         
            if i > head_num - 1:
                k = i - head_num
                line = line.strip('\n')
                array = line.split()
                if k % 2 == 1:      
                    t1 = int(array[7])
                    data = re.findall("\d+\D",array[5])
                    s = 0               
                    for j in data:
                        if j[-1] in CIGAR_set1:
                            s = s + int(j[0:-1])
                    t2 = int(array[3]) + s - 1
                    if t1 >= circ_spot[t][0] - 2 and t2 <= circ_spot[t][1] + 2:
                        out_sam.write(line1 + "\n")
                        out_sam.write(line + "\n")
                        sam_num[t] += 2
                line1 = line
            i = i + 1
        in_sam.close()
        out_sam.close()

    #construct DAGs
    exon_bin_0 = [None] * circ_spot_num
    exon_bin = [None] * circ_spot_num
    transcript_bin_idx = [None] * circ_spot_num
    transcript_bin = [None] * circ_spot_num
    subpath_constraints = [None] * circ_spot_num   
    edge = [None] * circ_spot_num
    new_edge = [None] * circ_spot_num
    isoform = [None] * circ_spot_num
    path_exon = [None] * circ_spot_num
    isoform = [None] * circ_spot_num
    isoform_num = [0] * circ_spot_num
    for t in range(0,circ_spot_num):
        exon_bin_0[t] = []
        transcript_bin[t] = []
        isoform[t] = []
        start = circ_spot[t][0]   
        end = circ_spot[t][1]    
        in_gtf = open(tempfile_path + gene_id + ".gtf","r")     
        for line in in_gtf:         
            line = line.strip('\n')
            array = line.split()
            if int(array[3]) >= start and int(array[4]) <= end:
                temp =  [0] * 2
                temp[0] =  int(array[3])
                temp[1] =  int(array[4])
                if temp not in exon_bin_0[t]:
                    exon_bin_0[t].append(temp)
        in_gtf.close()
        exon_bin_0[t].sort(key=lambda x:(x[0],x[1]))     
        exon_lb = []
        exon_rb = []
        for i in range(0,len(exon_bin_0[t])):
            exon_lb.append(exon_bin_0[t][i][0])
            exon_rb.append(exon_bin_0[t][i][1])
        if start in exon_lb and end in exon_rb: 
            while 1:
                if exon_bin_0[t][-1][-1] != end:
                    exon_bin_0[t].pop(-1)
                else:
                    break

            exon_bin[t] = copy.deepcopy(exon_bin_0[t])
                        

            if len(exon_bin[t]) == 1: 
                isoform[t] = [exon_bin[t]]
                isoform_num[t] = 1
            if len(exon_bin[t]) == 2:
                if exon_bin[t][0][0] == start and exon_bin[t][0][1] < exon_bin[t][1][0] and exon_bin[t][1][1] == end:
                    isoform[t] = [exon_bin[t]]
                    isoform_num[t] = 1
                elif exon_bin[t][0][0] == start and exon_bin[t][0][1] == end:
                    exon_bin[t].pop(1)
                    isoform[t] = [exon_bin[t]]
                    isoform_num[t] = 1
                elif exon_bin[t][1][0] == start and exon_bin[t][1][1] == end:
                    exon_bin[t].pop(0)
                    isoform[t] = [exon_bin[t]]
                    isoform_num[t] = 1

            if len(exon_bin[t]) >= 3 and exon_bin[t][0][0] < exon_bin[t][1][0] and exon_bin[t][-2][1] < exon_bin[t][-1][1]:
                node_num = len(exon_bin[t])
                ASG = np.zeros((node_num,node_num))
                subpath_constraints[t] = []
                in_sam = open(tempfile_path + gene_id + "_" + str(t) + ".sam","r")
                for line in in_sam:
                    line = line.strip('\n')
                    array = line.split()
                    data1 = re.findall("\d+",array[5])
                    data2 = re.findall("\D",array[5])
                    if data2.count("N") == 1:                       
                        N_id = data2.index("N")
                        match1 = 0
                        for i in range(0,N_id):
                            if data2[i] == "M" or data2[i] == "D":
                                match1 = match1 + int(data1[i])                     
                        match2 = 0
                        for i in range(N_id + 1,len(data1)):
                            if data2[i] == "M" or data2[i] == "D":
                                match2 = match2 + int(data1[i])
                        if match1 >= bound and match2 >= bound:    
                            t1 = int(array[3])
                            t2 = t1 + match1 - 1
                            t3 = t2 + int(data1[N_id]) + 1
                            t4 = t3 + match2 - 1
                            for i in range(0,node_num - 1):
                                if t1 >= exon_bin[t][i][0] - 2 and abs(t2 - exon_bin[t][i][1]) <= 2 :
                                    for j in range(i + 1,node_num):
                                        if abs(t3 - exon_bin[t][j][0]) <=2 and t4 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1                       
                    elif data2.count("N") == 2:     
                        N_id1 = data2.index("N")
                        data2[N_id1] = "X"
                        N_id2 = data2.index("N")
                        data2[N_id1] = "N"
                        match1 = 0
                        for i in range(0,N_id1):
                            if data2[i] == "M" or data2[i] == "D":
                                match1 = match1 + int(data1[i]) 
                        match2 = 0
                        for i in range(N_id1 + 1,N_id2):
                            if data2[i] == "M" or data2[i] == "D":
                                match2 = match2 + int(data1[i])
                        match3 = 0
                        for i in range(N_id2 + 1,len(data1)):
                            if data2[i] == "M" or data2[i] == "D":
                                match3 = match3 + int(data1[i])
                        t1 = int(array[3])
                        t2 = t1 + match1 - 1
                        t3 = t2 + int(data1[N_id1]) + 1
                        t4 = t3 + match2 - 1
                        t5 = t4 + int(data1[N_id2]) + 1
                        t6 = t5 + match3 - 1
                        if match1 < bound and match3 >= bound:      
                            for i in range(1,node_num - 1):
                                if abs(t3 - exon_bin[t][i][0]) <= 2 and abs(t4 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num):
                                        if abs(t5 - exon_bin[t][j][0]) <= 2 and t6 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1                           
                        elif match1 >= bound and match3 < bound:  
                            for i in range(0,node_num - 2):
                                if t1 >= exon_bin[t][i][0] - 2 and abs(t2 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t3 - exon_bin[t][j][0]) <= 2 and abs(t4 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                        elif match1 >= bound and match3 >= bound:  
                            for i in range(0,node_num - 2):
                                if t1 >= exon_bin[t][i][0] - 2 and abs(t2 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t3 - exon_bin[t][j][0]) <= 2 and abs(t4 - exon_bin[t][j][1]) <= 2:
                                            for k in range(j + 1,node_num):
                                                if abs(t5 - exon_bin[t][k][0]) <= 2 and t6 <= exon_bin[t][k][1] + 2:
                                                    ASG[i,j] += 1
                                                    ASG[j,k] += 1
                                                    temp =  [None] * 3
                                                    temp[0] =  exon_bin[t][i]
                                                    temp[1] =  exon_bin[t][j]
                                                    temp[2] =  exon_bin[t][k]
                                                    if temp not in subpath_constraints[t]:
                                                        subpath_constraints[t].append(temp)
                in_sam.close()
                m = 0
                in_sam = open(tempfile_path + gene_id + "_" + str(t) + ".sam","r")
                while m != sam_num[t]:
                    line1 = in_sam.readline()
                    line1 = line1.strip('\n')
                    array1 = line1.split()
                    data11 = re.findall("\d+",array1[5])
                    data12 = re.findall("\D",array1[5])
                    line2 = in_sam.readline()
                    line2 = line2.strip('\n')
                    array2 = line2.split()
                    data21 = re.findall("\d+",array2[5])
                    data22 = re.findall("\D",array2[5])
                    if data12.count("N") == 0 and data22.count("N") == 0:    
                        match1 = 0
                        for i in range(0,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match1 = match1 + int(data11[i])                        
                        match2 = 0
                        for i in range(0,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match2 = match2 + int(data21[i])                        
                        t11 = int(array1[3])
                        t12 = t11 + match1 - 1
                        t21 = int(array2[3])
                        t22 = t21 + match2 - 1
                        for i in range(0,node_num - 1):
                            if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2:
                                flag = 0
                                for j in range(i + 1,node_num):
                                    if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                        flag = 1
                                        break
                                if flag == 1:
                                    if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                    if j < node_num - 1:
                                        if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                            j = j + 1
                                            if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                                ASG[i,j] += 1
                                    
                    elif data12.count("N") == 1 and data22.count("N") == 1:         
                        N_id1 = data12.index("N")
                        match11 = 0
                        for i in range(0,N_id1):
                            if data12[i] == "M" or data12[i] == "D":
                                match11 = match11 + int(data11[i])                      
                        match12 = 0
                        for i in range(N_id1 + 1,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match12 = match12 + int(data11[i])
                        t11 = int(array1[3])
                        t12 = t11 + match11 - 1
                        t13 = t12 + int(data11[N_id1]) + 1
                        t14 = t13 + match12 - 1
                        N_id2 = data22.index("N")
                        match21 = 0
                        for i in range(0,N_id2):
                            if data22[i] == "M" or data22[i] == "D":
                                match21 = match21 + int(data21[i])                      
                        match22 = 0
                        for i in range(N_id2 + 1,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match22 = match22 + int(data21[i])
                        t21 = int(array2[3])
                        t22 = t21 + match21 - 1
                        t23 = t22 + int(data21[N_id2]) + 1
                        t24 = t23 + match22 - 1                     

                        if match11 >= bound and match22 >= bound:   
                            for i in range(0,node_num - 2):     
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            for k in range(j + 1,node_num):
                                                if abs(t23 - exon_bin[t][k][0]) <= 2 and t24 <= exon_bin[t][k][1] + 2:
                                                    temp =  [None] * 3
                                                    temp[0] =  exon_bin[t][i]
                                                    temp[1] =  exon_bin[t][j]
                                                    temp[2] =  exon_bin[t][k]
                                                    if temp not in subpath_constraints[t]:
                                                        subpath_constraints[t].append(temp)
                                                    if match12 < bound:
                                                        ASG[i,j] += 1
                                                    if match21 < bound:
                                                        ASG[j,k] += 1
                            if match12 >= bound and match21 >= bound: 
                                for i in range(0,node_num - 3):
                                    if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                        for j in range(i + 1,node_num - 2):
                                            if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2:
                                                flag = 0
                                                for k in range(j + 1,node_num - 1):
                                                    if exon_bin[t][k][0] >= exon_bin[t][j][1]:
                                                        flag = 1
                                                        break
                                                if flag == 1:                                           
                                                    if t21 >= exon_bin[t][k][0] - 2 and abs(t22 - exon_bin[t][k][1]) <= 2:
                                                        for l in range(k + 1,node_num):
                                                            if abs(t23 - exon_bin[t][l][0]) <= 2 and t24 <= exon_bin[t][l][1] + 2:
                                                                ASG[j,k] += 1
                                                                temp =  [None] * 4
                                                                temp[0] =  exon_bin[t][i]
                                                                temp[1] =  exon_bin[t][j]
                                                                temp[2] =  exon_bin[t][k]
                                                                temp[3] =  exon_bin[t][l]
                                                                if temp not in subpath_constraints[t]:
                                                                    subpath_constraints[t].append(temp)
                                                    if k < node_num - 2:
                                                        if exon_bin[t][k][0] == exon_bin[t][k + 1][0]:
                                                            k = k + 1
                                                            if t21 >= exon_bin[t][k][0] - 2 and abs(t22 - exon_bin[t][k][1]) <= 2:
                                                                for l in range(k + 1,node_num):
                                                                    if abs(t23 - exon_bin[t][l][0]) <= 2 and t24 <= exon_bin[t][l][1] + 2:
                                                                        ASG[j,k] += 1
                                                                        temp =  [None] * 4
                                                                        temp[0] =  exon_bin[t][i]
                                                                        temp[1] =  exon_bin[t][j]
                                                                        temp[2] =  exon_bin[t][k]
                                                                        temp[3] =  exon_bin[t][l]
                                                                        if temp not in subpath_constraints[t]:
                                                                            subpath_constraints[t].append(temp)         
                                                                            
                        if match12 < bound and match22 < bound:   
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    ASG[i,j] += 1                                           
                                    
                        if match11 < bound and match21 < bound:   
                            for i in range(1,node_num - 1):
                                if abs(t13 - exon_bin[t][i][0]) <= 2 and t14 <= exon_bin[t][i][1] + 2 and t21 >= exon_bin[t][i][0] - 2 and abs(t22 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1

                        if match11 < bound and match22 < bound:   
                            for i in range(1,node_num - 2):
                                if abs(t13 - exon_bin[t][i][0]) <= 2 and t14 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 <= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 <= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    ASG[i,j] += 1
                                        

                    elif data12.count("N") == 1 and data22.count("N") == 0:   
                        N_id1 = data12.index("N")
                        match11 = 0
                        for i in range(0,N_id1):
                            if data12[i] == "M" or data12[i] == "D":
                                match11 = match11 + int(data11[i])                      
                        match12 = 0
                        for i in range(N_id1 + 1,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match12 = match12 + int(data11[i])
                        t11 = int(array1[3])
                        t12 = t11 + match11 - 1
                        t13 = t12 + int(data11[N_id1]) + 1
                        t14 = t13 + match12 - 1
                        match2 = 0
                        for i in range(0,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match2 = match2 + int(data21[i])                        
                        t21 = int(array2[3])
                        t22 = t21 + match2 - 1
                        if match11 >= bound and match12 >= bound:
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    for j in range(i + 1,node_num - 1):
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2:
                                            flag = 0
                                            for k in range(j + 1,node_num):
                                                if exon_bin[t][k][0] >= exon_bin[t][j][1]:
                                                    flag = 1
                                                    break
                                            if flag == 1:
                                                if t21 >= exon_bin[t][k][0] - 2 and t22 <= exon_bin[t][k][1] + 2:
                                                        ASG[j,k] += 1   
                                                        temp =  [None] * 3
                                                        temp[0] =  exon_bin[t][i]
                                                        temp[1] =  exon_bin[t][j]
                                                        temp[2] =  exon_bin[t][k]
                                                        if temp not in subpath_constraints[t]:
                                                            subpath_constraints[t].append(temp)
                                                if k < node_num - 1:
                                                    if exon_bin[t][k][0] == exon_bin[t][k + 1][0]:
                                                        k = k + 1
                                                        if t21 >= exon_bin[t][k][0] - 2 and t22 <= exon_bin[t][k][1] + 2:
                                                            ASG[j,k] += 1   
                                                            temp =  [None] * 3
                                                            temp[0] =  exon_bin[t][i]
                                                            temp[1] =  exon_bin[t][j]
                                                            temp[2] =  exon_bin[t][k]
                                                            if temp not in subpath_constraints[t]:
                                                                subpath_constraints[t].append(temp)
                                                                                                
                        elif match11 >= bound and match12 < bound:
                            for i in range(0,node_num - 1):
                                if t11 >= exon_bin[t][i][0] - 2 and abs(t12 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and  t22 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t13 - exon_bin[t][j][0]) <= 2 and t14 <= exon_bin[t][j][1] + 2 and t21 >= exon_bin[t][j][0] - 2 and  t22 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1
                                                                                
                        elif match11 < bound and match12 >= bound:
                            for i in range(1,node_num - 1):
                                if abs(t13 - exon_bin[t][i][0]) <= 2 and t14 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 >= exon_bin[t][j][0] - 2 and t22 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1   
                                                    
                    elif data12.count("N") == 0 and data22.count("N") == 1:     
                        match1 = 0
                        for i in range(0,len(data12)):
                            if data12[i] == "M" or data12[i] == "D":
                                match1 = match1 + int(data11[i])                        
                        t11 = int(array1[3])
                        t12 = t11 + match1 - 1                      
                        N_id2 = data22.index("N")
                        match21 = 0
                        for i in range(0,N_id2):
                            if data22[i] == "M" or data22[i] == "D":
                                match21 = match21 + int(data21[i])                      
                        match22 = 0
                        for i in range(N_id2 + 1,len(data22)):
                            if data22[i] == "M" or data22[i] == "D":
                                match22 = match22 + int(data21[i])
                        t21 = int(array2[3])
                        t22 = t21 + match21 - 1
                        t23 = t22 + int(data21[N_id2]) + 1
                        t24 = t23 + match22 - 1     

                        if match21 >= bound and match22 >= bound:
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            for k in range(j + 1,node_num):
                                                if abs(t23 - exon_bin[t][k][0]) <= 2 and t24 <= exon_bin[t][k][1] + 2:
                                                    ASG[i,j] += 1
                                                    temp =  [None] * 3
                                                    temp[0] =  exon_bin[t][i]
                                                    temp[1] =  exon_bin[t][j]
                                                    temp[2] =  exon_bin[t][k]
                                                    if temp not in subpath_constraints[t]:
                                                        subpath_constraints[t].append(temp)
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    for k in range(j + 1,node_num):
                                                        if abs(t23 - exon_bin[t][k][0]) <= 2 and t24 <= exon_bin[t][k][1] + 2:
                                                            ASG[i,j] += 1
                                                            temp =  [None] * 3
                                                            temp[0] =  exon_bin[t][i]
                                                            temp[1] =  exon_bin[t][j]
                                                            temp[2] =  exon_bin[t][k]
                                                            if temp not in subpath_constraints[t]:
                                                                subpath_constraints[t].append(temp)                                         
                                                    
                        elif match21 >= bound and match22 < bound:
                            for i in range(0,node_num - 2):
                                if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1
                                            break
                                    if flag == 1:
                                        if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 2:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if t21 >= exon_bin[t][j][0] - 2 and abs(t22 - exon_bin[t][j][1]) <= 2:
                                                    ASG[i,j] += 1
                                                    
                        elif match21 < bound and match22 >= bound:
                            for i in range(0,node_num - 1):
                                if t11 >= exon_bin[t][i][0] - 2 and t12 <= exon_bin[t][i][1] + 2 and t21 >= exon_bin[t][i][0] - 2 and abs(t22 - exon_bin[t][i][1]) <= 2:
                                    flag = 0
                                    for j in range(i + 1,node_num - 1):
                                        if exon_bin[t][j][0] >= exon_bin[t][i][1]:
                                            flag = 1                                    
                                            break
                                    if flag == 1:
                                        if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                            ASG[i,j] += 1
                                        if j < node_num - 1:
                                            if exon_bin[t][j][0] == exon_bin[t][j + 1][0]:
                                                j = j + 1
                                                if abs(t23 - exon_bin[t][j][0]) <= 2 and t24 <= exon_bin[t][j][1] + 2:
                                                    ASG[i,j] += 1
                                    
                                    
                    m += 2
                in_sam.close()
                for i in range(0,node_num):
                    for j in range(0,node_num):
                        if ASG[i,j] < 3:
                            ASG[i,j] = 0
                        else:
                            ASG[i,j] = 1          
                
                
                while 1:
                    isolate = []
                    for i in range(1,node_num - 1):
                        if sum(ASG[:,i]) == 0 or sum(ASG[i,:]) == 0:
                            isolate.append(i)
                    if isolate != []:
                        isolate.reverse()
                        for i in isolate:
                            exon_bin[t].pop(i)      
                        ASG = np.delete(ASG,isolate,0)
                        ASG = np.delete(ASG,isolate,1)
                        node_num = len(exon_bin[t])
                    else:
                        break
                
                if sum(ASG[0,:]) > 0 and sum(ASG[:,-1]) > 0:                        
                    for i in range(len(subpath_constraints[t]) - 1,-1,-1):
                        for j in range(0,len(subpath_constraints[t][i])):
                            if subpath_constraints[t][i][j] not in exon_bin[t]:
                                subpath_constraints[t].pop(i)
                                break
                    node_num_not_isolated1 = len(exon_bin[t])  
                
                    edge[t] = []
                    for i in range(0,node_num_not_isolated1):
                        for j in range(0,node_num_not_isolated1):
                            if ASG[i,j] == 1:
                                temp = [0] * 2
                                temp[0] = exon_bin[t][i]
                                temp[1] = exon_bin[t][j]
                                edge[t].append(temp)


                    subpath_constraints[t].sort(key=lambda x:(x[0],x[1],x[2]))  
                    for i in range(len(subpath_constraints[t]) - 1,0,-1):
                        flag = 0
                        for j in range(i - 1,-1,-1):
                            if subpath_constraints[t][i][0] == subpath_constraints[t][j][-1]:       
                                t1 = exon_bin[t].index(subpath_constraints[t][j][-1])
                                t2 = exon_bin[t].index(subpath_constraints[t][j][-2])
                                s = 0
                                for k in range(0,t1):                     
                                    s = s + ASG[k,t1]
                                if s == 1 and ASG[t2,t1] == 1:
                                    flag = 1
                                    for l in range(1,len(subpath_constraints[t][i])):
                                        subpath_constraints[t][j].append(subpath_constraints[t][i][l])
                            elif subpath_constraints[t][i][0] == subpath_constraints[t][j][-2] and subpath_constraints[t][i][1] == subpath_constraints[t][j][-1]: 
                                t1 = exon_bin[t].index(subpath_constraints[t][j][-2])
                                t2 = exon_bin[t].index(subpath_constraints[t][j][-3])
                                s = 0
                                for k in range(0,t1):
                                    s = s + ASG[k,t1]
                                if s == 1 and ASG[t2,t1] == 1:
                                    flag = 1
                                    for l in range(2,len(subpath_constraints[t][i])):
                                        subpath_constraints[t][j].append(subpath_constraints[t][i][l])
                            elif subpath_constraints[t][i][0] == subpath_constraints[t][j][-3] and subpath_constraints[t][i][1] == subpath_constraints[t][j][-2] and subpath_constraints[t][i][2] == subpath_constraints[t][j][-1]:   
                                if len(subpath_constraints[t][j]) == 3:
                                    flag = 1
                                elif len(subpath_constraints[t][j]) == 4:
                                    t1 = exon_bin[t].index(subpath_constraints[t][j][-3])
                                    t2 = exon_bin[t].index(subpath_constraints[t][j][-4])
                                    s = 0
                                    for k in range(0,t1):                     
                                        s = s + ASG[k,t1]
                                    if s == 1 and ASG[t2,t1] == 1:
                                        flag = 1
                                        for l in range(3,len(subpath_constraints[t][i])):
                                            subpath_constraints[t][j].append(subpath_constraints[t][i][l])
                        if flag == 1:               #
                            subpath_constraints[t][i] = []

                    for i in range(len(subpath_constraints[t]) - 1,0,-1):
                        if subpath_constraints[t][i] == []:
                            subpath_constraints[t].pop(i)
                    temp = copy.deepcopy(subpath_constraints[t])
                    subpath_constraints[t] = []
                    for i in temp:
                        if i not in subpath_constraints[t]:
                            subpath_constraints[t].append(i)
                    new_edge[t] = edge[t] + subpath_constraints[t]
                    ASG_0 = copy.deepcopy(ASG)
                    exon_bin_1 = copy.deepcopy(exon_bin[t])
                    for i in range(0,len(subpath_constraints[t])):
                        for j in range(0,len(subpath_constraints[t][i]) - 1):
                            t1 = exon_bin[t].index(subpath_constraints[t][i][j])
                            t2 = exon_bin[t].index(subpath_constraints[t][i][j + 1])
                            ASG[t1,t2] = 0
                        t1 = exon_bin[t].index(subpath_constraints[t][i][0])
                        t2 = exon_bin[t].index(subpath_constraints[t][i][-1])
                        ASG[t1,t2] = 1
                           

                    
                    isolate = []
                    for i in range(1,node_num_not_isolated1 - 1):
                        if sum(ASG[:,i]) == 0 and sum(ASG[i,:]) == 0:
                            isolate.append(i)
                    isolate.reverse()
                    
                    for i in isolate:
                        exon_bin[t].pop(i)      
                    ASG = np.delete(ASG,isolate,0)
                    ASG = np.delete(ASG,isolate,1)
                    
                    node_num_not_isolated2 = len(exon_bin[t])
                    edge_num =  int(np.sum(ASG))
                    new_node_num = node_num_not_isolated2 + edge_num
            
                    while np.sum(ASG) < 4 * edge_num:   
                        temp = 0
                        for i in range(0,len(ASG)):
                            for j in range(len(ASG) - 1,-1,-1):
                                if int(ASG[i,j]) == 1:
                                    temp = 1
                                    break 
                            if temp == 1:                       
                                break

                        ASG[i,j] = 0
                        ASG = np.insert(ASG,i + 1,values = 0,axis = 0)                  
                        ASG = np.insert(ASG,i + 1,values = 0,axis = 1)
                        ASG[i,i + 1] = 2
                        ASG[i + 1,j + 1] = 2
                        exon_bin[t].insert(i + 1,[])
                    for i in range(0,new_node_num):                 
                        for j in range(0,new_node_num):
                            if ASG[i,j] == 2:
                                ASG[i,j] = 1
                

                    ASG_tc = copy.deepcopy(ASG)      

                    for k in range(0,new_node_num):
                        for i in range(0,new_node_num):
                            for j in range(0,new_node_num):
                                if ASG_tc[i,k] + ASG_tc[k,j] == 2:
                                    ASG_tc[i,j] = 1

                    M = max_match(ASG_tc)                 
                    path = []         
                    x = []
                    for i in range(0,new_node_num):
                        x.append(1)             
                    while sum(x) > 0:          
                        a = x.index(1)     
                        x[a] = 0     
                        temp = [a]
                        while sum(M[a,:]) > 0:
                            for j in range(0,new_node_num):
                                if M[a,j] == 1:
                                    break
                            x[j] = 0
                            temp.append(j)
                            a = j
                        path.append(temp)

                    for i in range(0,len(path)):  
                        a = path[i][0]
                        b = path[i][-1]
                        if exon_bin[t][a] == []:
                            for j in range(0,new_node_num):
                                if ASG[j,a] != 0:
                                    break
                            path[i].insert(0,j) 
                        if exon_bin[t][b] == []:
                            for j in range(0,new_node_num):
                                if ASG[b,j] != 0:
                                    break
                            path[i].append(j)
                    for i in range(len(path) - 1,-1,-1):
                        for j in range(len(path[i]) - 1,0,-1):    
                            a = path[i][j] 
                            b = path[i][j - 1]
                            if exon_bin[t][a] == [] and exon_bin[t][b] == []:
                                path.pop(i)
                                break

                                

                    path_exon[t] = [None] * len(path)           
                    for i in range(0,len(path)):
                        temp = []
                        path_exon[t][i] = [[]]
                        for j in range(0,len(path[i])):
                            a = path[i][j]
                            temp.append(exon_bin[t][a])             
                        path_exon[t][i][0].append(temp[0])
                        for j in range(1,len(path[i])):
                            if temp[j] != []:
                                if temp[j - 1] == []:
                                    for k in path_exon[t][i]:
                                        k.append(temp[j])
                                else:
                                    a = exon_bin_1.index(temp[j - 1])
                                    b = exon_bin_1.index(temp[j])
                                    if ASG_0[a,b] == 1:
                                        for k in path_exon[t][i]:
                                            k.append(temp[j])
                                    else:
                                        ab_path = findPath(ASG_0,[a],b)
                                        path_len = []
                                        for k in ab_path:
                                            path_len.append(len(k))
                                        min_num = min(path_len)
                                        ab_min_path = ab_path[path_len.index(min_num)]
                                        for k in path_exon[t][i]:
                                            for m in range(1,len(ab_min_path)):
                                                k.append(exon_bin_1[ab_min_path[m]])
                                                
                            else:
                                a =  temp[j - 1]
                                b =  temp[j + 1]
                                l = 0
                                temp1 = []
                                for k in range(0,len(new_edge[t])):
                                    if new_edge[t][k][0] == a and new_edge[t][k][-1] == b:
                                        l = l + 1
                                        temp2 = copy.deepcopy(new_edge[t][k])
                                        temp2.pop(0)
                                        temp2.pop(-1)
                                        temp1.append(temp2)
                                p = len(path_exon[t][i])
                                path_exon[t][i] = path_exon[t][i] * l
                                for n in range(0,l):
                                    for m in range(0,p):
                                         q = n * p + m
                                         path_exon[t][i][q] = path_exon[t][i][q] + temp1[n]
                                    
                
                    temp = copy.deepcopy(path_exon[t])
                    path_exon[t] = []
                    for i in temp:
                        for j in i:
                            path_exon[t].append(j)
                    path_head = []
                    path_tail = []
                    for k in range(0,len(path_exon[t])):
                        if path_exon[t][k][0][0] != start and path_exon[t][k][0] not in path_head:
                            path_head.append(path_exon[t][k][0])
                        if path_exon[t][k][-1][-1] != end and path_exon[t][k][-1] not in path_tail:
                            path_tail.append(path_exon[t][k][-1])           
                    #head_to_start   
                    if path_head != []:
                        head_to_start =  [None] * len(path_head)
                        for k in range(0,len(path_head)):
                            b = exon_bin_1.index(path_head[k])
                            head_to_start[k] = findPath(ASG_0,[0],b)
                            path_len = []
                            for i in range(0,len(head_to_start[k])):
                                path_len.append(len(head_to_start[k][i]))
                            min_num = min(path_len)
                            for i in range(len(head_to_start[k]) - 1,-1,-1):
                                if len(head_to_start[k][i]) > min_num:
                                    head_to_start[k].pop(i)                 
                            for i in range(0,len(head_to_start[k])):
                                for j in range(0,len(head_to_start[k][i])):
                                    head_to_start[k][i][j] = exon_bin_1[head_to_start[k][i][j]]             
                    #tail_to_end    
                    if path_tail != []:
                        tail_to_end = [None] * len(path_tail)
                        for k in range(0,len(path_tail)):
                            a = exon_bin_1.index(path_tail[k])
                            tail_to_end[k] = findPath(ASG_0,[a],len(exon_bin_1) - 1)
                            path_len = []
                            for i in range(0,len(tail_to_end[k])):
                                path_len.append(len(tail_to_end[k][i]))
                            min_num = min(path_len)
                            for i in range(len(tail_to_end[k]) - 1,-1,-1):
                                if len(tail_to_end[k][i]) > min_num:
                                    tail_to_end[k].pop(i)                   
                            for i in range(0,len(tail_to_end[k])):
                                for j in range(0,len(tail_to_end[k][i])):
                                    tail_to_end[k][i][j] = exon_bin_1[tail_to_end[k][i][j]]

                    transcript_bin[t] = []
                    for k in range(0,len(path_exon[t])):
                        if path_exon[t][k][0][0] ==  start and path_exon[t][k][-1][-1] == end: 
                            transcript_bin[t].append(path_exon[t][k])
                        elif path_exon[t][k][0][0] !=  start and path_exon[t][k][-1][-1] == end: 
                            a = path_head.index(path_exon[t][k][0])
                            path_exon[t][k].pop(0)
                            for i in range(0,len(head_to_start[a])):                        
                                transcript_bin[t].append(head_to_start[a][i] + path_exon[t][k])
                        elif path_exon[t][k][0][0] ==  start and path_exon[t][k][-1][-1] != end: 
                            a = path_tail.index(path_exon[t][k][-1])
                            path_exon[t][k].pop(-1)
                            for i in range(0,len(tail_to_end[a])):                      
                                transcript_bin[t].append(path_exon[t][k] + tail_to_end[a][i])
                        elif path_exon[t][k][0][0] !=  start and path_exon[t][k][-1][-1] != end:
                            a = path_head.index(path_exon[t][k][0])
                            b = path_tail.index(path_exon[t][k][-1])
                            if len(path_exon[t][k]) == 1:
                                path_exon[t][k].pop(0)
                            else:
                                path_exon[t][k].pop(0)
                                path_exon[t][k].pop(-1)
                            temp = []
                            p = len(head_to_start[a])
                            l = len(tail_to_end[b])
                            for i in range(0,p):                        
                                temp.append(head_to_start[a][i] + path_exon[t][k])
                            temp = temp * l
                            for n in range(0,l):
                                for m in range(0,p):
                                    q = n * p + m
                                    transcript_bin[t].append(temp[q] + tail_to_end[b][n])
                    isoform[t] = []
                    for i in transcript_bin[t]:
                        if i not in isoform[t]:
                            isoform[t].append(i)
                    isoform_num[t] = len(isoform[t])    
    isoform_num_0 = copy.deepcopy(isoform_num)
    for t in range(circ_spot_num - 1,-1,-1):
        if isoform_num[t] == 0:
            isoform_num.pop(t)
            isoform.pop(t)
    isoform_num_sum = sum(isoform_num)
    isoform_len = []
    matrix = []
    if isoform_num_sum > 0:
        for t in range(0,len(isoform_num)):     
            for i in range(0,isoform_num[t]):
                matrix.append(isoform[t][i])     
                length = 0
                for j in range(0,len(isoform[t][i])):        
                    length = length + int(isoform[t][i][j][1]) - int(isoform[t][i][j][0]) + 1
                isoform_len.append(length)
    
    os.remove(tempfile_path + gene_id + "_A.sam")
    os.remove(tempfile_path + gene_id + "_B.sam")
    os.remove(tempfile_path + gene_id + ".sam")
    
    
    ##############################定量
    sam_row_num_A = 0
    i = 0
    in_sam = open(tempfile_path + "selected_" + chr_name + "_A.sam", "r")
    out_sam = open(tempfile_path + gene_id + "_A.sam", "a")
    for line in in_sam:
        line = line.strip('\n')
        if i < head_num:
            out_sam.write(line + "\n")
        else:
            k = i - head_num
            line = line.strip('\n')
            array = line.split()
            # print('array0 = ',array[0],'k = ',k,'i = ', i)
            if k % 2 == 1:  ###################################取右端
                t1 = int(array[7])
                data = re.findall("\d+\D", array[5])
                s = 0
                for j in data:
                    if j[-1] in CIGAR_set1:
                        s = s + int(j[0:-1])
                t2 = int(array[3]) + s - 1
                if t1 >= gene_coordinate_min - 2 and t2 <= gene_coordinate_max + 2:  #############如果双端reads起始终止在该基因区间内，写出
                    out_sam.write(line1 + "\n")
                    out_sam.write(line + "\n")
                    sam_row_num_A += 2
            line1 = line  #########################################记录左端的line
        i = i + 1
    in_sam.close()
    out_sam.close()
    shutil.copyfile(tempfile_path + gene_id + "_A.sam", tempfile_path + gene_id + ".sam")
    #####################################################################
    sam_row_num_B = 0
    i = 0
    in_sam = open(tempfile_path + "selected_" + chr_name + "_B.sam", "r")
    out_sam1 = open(tempfile_path + gene_id + "_B.sam", "w")
    out_sam2 = open(tempfile_path + gene_id + ".sam", "a")
    for line in in_sam:
        line = line.strip('\n')
        if i < head_num:
            out_sam1.write(line + "\n")
        else:
            k = i - head_num
            line = line.strip('\n')
            array = line.split()
            if k % 2 == 1:
                t1 = int(array[7])
                data = re.findall("\d+\D", array[5])
                s = 0
                for j in data:
                    if j[-1] in CIGAR_set1:
                        s = s + int(j[0:-1])
                t2 = int(array[3]) + s - 1
                if t1 >= gene_coordinate_min - 2 and t2 <= gene_coordinate_max + 2:
                    out_sam1.write(line1 + "\n")
                    out_sam1.write(line + "\n")
                    sam_row_num_B += 2
                    out_sam2.write(line1 + "\n")
                    out_sam2.write(line + "\n")
            line1 = line
        i = i + 1
    in_sam.close()
    out_sam1.close()
    out_sam2.close()
    #################################################################把circRNA位点信息按起始和终止分为不同的spot，
    circ_spot = []
    circ_num = []
    strand = []
    in_txt = open(tempfile_path + "circRNA_list_selected.txt", "r")
    for line in in_txt:
        line = line.strip('\n')
        array = line.split()
        if array[0] == chr_name and array[4] == gene_id:  #############################这里确定取出来的spot位点是在本chrname和geneid上
            temp = [0] * 2
            temp[0] = int(array[1]) + 1
            temp[1] = int(array[2])
            circ_spot.append(temp)
            circ_num.append(int(float(array[6])))
            strand.append(array[3])
    in_txt.close()
    ciri_spot_num = len(circ_spot)
    #####################################################################################################按线性gtf结果得到matrix_liner,isoform_len_liner
    liner_spot = []
    strand_liner = []
    in_gtf = open(tempfile_path + chr_name + "_stringtie.gtf", "r")
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        if array[2] == "transcript" and int(array[3]) >= gene_coordinate_min - 2 and int(array[4]) <= gene_coordinate_max + 2:  #############################这里确定取出来的spot位点是在本chrname和geneid上
            temp = [0] * 2
            temp[0] = int(array[3])
            temp[1] = int(array[4])
            liner_spot.append(temp)
            strand_liner.append(array[6])
    in_gtf.close()
    liner_spot_num = len(liner_spot)

    matrix_liner = [None] * len(liner_spot)
    isoform_len_liner = []
    in_gtf = open(tempfile_path + chr_name + "_stringtie.gtf", "r")
    t = 0
    mark = 0
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        if mark == 1 and array[2] == "exon":
            temp = [0] * 2
            temp[0] = int(array[3])
            temp[1] = int(array[4])
            matrix_liner[t].append(temp)
            lenth = lenth + int(array[4]) - int(array[3]) + 1
        if mark == 1 and array[2] == "transcript":
            mark = 0
            isoform_len_liner.append(lenth)
            t = t + 1
        if array[2] == "transcript" and int(array[3]) >= gene_coordinate_min - 2 and int(array[4]) <= gene_coordinate_max + 2:  #############################这里确定取出来的liner_spot位点是在本chrname和geneid上
            mark = 1
            lenth = 0
            matrix_liner[t] = []
    in_gtf.close()
    if len(isoform_len_liner) != len(matrix_liner):
        isoform_len_liner.append(lenth)
    ###################################################################################################################按circAST结果得到matrix,isoform_len
    circAST_spot_num = len(matrix)
    
    
    matrix_mixed = matrix + matrix_liner
    isoform_len_mixed = isoform_len + isoform_len_liner
    length_first_calculate = len(isoform_len) + len(isoform_len_liner)
    # print("circAST_spot_num    =    ", circAST_spot_num,"   liner_spot_num   =    ",len(liner_spot))
    circ_reads, FPKM, matrix_mixed_, isoform_len_mixed_, circAST_spot_num, liner_spot_num= selection_mixed(matrix_mixed, isoform_len_mixed,circAST_spot_num,liner_spot_num)
    if(matrix_mixed_ == 0):
        os.remove(tempfile_path + gene_id + ".gtf")
        os.remove(tempfile_path + gene_id + "_A.sam")
        os.remove(tempfile_path + gene_id + "_B.sam")
        os.remove(tempfile_path + gene_id + ".sam")
        return 0
    # split = circAST_spot_num ##怕第二次定量circAST_spot_num会变
    print("length_first_calculate = " + str(length_first_calculate) + ",len(matrix_mixed_) = " + str(len(matrix_mixed_)))
    if length_first_calculate != len(matrix_mixed_):  ####这样来实现筛选假阳性
        circ_reads, FPKM, matrix_mixed_, isoform_len_mixed_, circAST_spot_num, liner_spot_num= selection_mixed(matrix_mixed_, isoform_len_mixed_,circAST_spot_num,liner_spot_num)
    if sum(FPKM) > 0:
        # exon_num = [0] * len(selected_isoform_circ)
        # for t in range(0, len(selected_isoform_circ)):
        #     exon_num[t] = len(selected_isoform_circ[t])
        #
        # matrix_calculate = []
        # for t in range(0, len(selected_isoform_circ)):
        #     for i in range(0, exon_num[t]):
        #         matrix_calculate.append(selected_isoform_circ[t][i])
        # print(len(matrix_calculate))
        # matrix_0 = copy.deepcopy(matrix_calculate)
        # matrix_new = []
        # FPKM_new = []
        # isoform_len_new = []
        #
        # for t in range(0, len(exon_num)):
        #     matrix_temp = []
        #     FPKM_temp = np.zeros(exon_num[t])
        #     for i in range(0, exon_num[t]):
        #         matrix_temp.append(matrix_calculate[0])
        #         FPKM_temp[i] = FPKM[0]
        #         matrix_calculate.pop(0)
        #         FPKM = np.delete(FPKM, 0)
        #     FPKM_temp_max = np.max(FPKM_temp)
        #     for i in range(0, exon_num[t]):
        #         if FPKM_temp[i] > FPKM_temp_max * 0.001:
        #             matrix_new.append(matrix_temp[i])
        #             FPKM_new = np.concatenate((FPKM_new, [FPKM_temp[i]]))
        #             temp = matrix_0.index(matrix_temp[i])
        #             isoform_len_new.append(isoform_len_mixed[temp])

        results = open("result_table.txt", "a")
        results_liner = open("result_liner.txt","a")
        for i in range(0, len(matrix_mixed_)):
            if i < circAST_spot_num:
                temp = [0] * 2
                temp[0] = matrix_mixed_[i][0][0]
                temp[1] = matrix_mixed_[i][-1][-1]
                junction_reads = circ_num[circ_spot.index(temp)]
                strd = strand[circ_spot.index(temp)]
                exon_len = []
                exon_offset = []
                for k in range(0, len(matrix_mixed_[i])):
                    exon_len.append(int(matrix_mixed_[i][k][1]) - int(matrix_mixed_[i][k][0]) + 1)
                    exon_offset.append(int(matrix_mixed_[i][k][0]) - int(matrix_mixed_[i][0][0]))
                match_frag = round(FPKM[i] * sum(exon_len) * match_num_total * (10 ** (-9)))
                match_frag = int(match_frag)
                if match_frag == 0:
                    match_frag = 1
                for k in range(0, len(exon_len)):
                    exon_len[k] = str(exon_len[k])
                for k in range(0, len(exon_offset)):
                    exon_offset[k] = str(exon_offset[k])
                results.write(chr_name + "\t" + str(matrix_mixed_[i][0][0]) + "\t" + str(matrix_mixed_[i][-1][-1]) + "\t" + gene_id + "\t" + str(strd) + "\t" + str(len(matrix_mixed_[i])) + "\t" + ','.join(exon_len) + "\t" + ','.join(exon_offset) + "\t" + str(junction_reads) + "\t" + str(FPKM[i]) + "\t" + str(match_frag) + "\t" + "circ" + "\n")
            if i >= circAST_spot_num:
                temp = [0] * 2
                temp[0] = matrix_mixed_[i][0][0]
                temp[1] = matrix_mixed_[i][-1][-1]
                # junction_reads = circ_num[circ_spot.index(temp)]
                strd_liner = strand_liner[liner_spot.index(temp)]
                exon_len = []
                exon_offset = []
                for k in range(0, len(matrix_mixed_[i])):
                    exon_len.append(int(matrix_mixed_[i][k][1]) - int(matrix_mixed_[i][k][0]) + 1)
                    exon_offset.append(int(matrix_mixed_[i][k][0]) - int(matrix_mixed_[i][0][0]))
                match_frag = round(FPKM[i] * sum(exon_len) * match_num_total * (10 ** (-9)))
                match_frag = int(match_frag)
                if match_frag == 0:
                    match_frag = 1
                for k in range(0, len(exon_len)):
                    exon_len[k] = str(exon_len[k])
                for k in range(0, len(exon_offset)):
                    exon_offset[k] = str(exon_offset[k])
                results_liner.write(chr_name + "\t" + str(matrix_mixed_[i][0][0]) + "\t" + str(matrix_mixed_[i][-1][-1]) + "\t" + gene_id + "\t" + str(strd_liner) + "\t" + str(len(matrix_mixed_[i])) + "\t" + ','.join(exon_len) + "\t" + ','.join(exon_offset) + "\t" + "NA" + "\t" + str(FPKM[i]) + "\t" + str(match_frag) + "\n")
        results.close()
        results_liner.close()
    os.remove(tempfile_path + gene_id + ".gtf")
    os.remove(tempfile_path + gene_id + "_A.sam")
    os.remove(tempfile_path + gene_id + "_B.sam")
    os.remove(tempfile_path + gene_id + ".sam")
    for t in range(0,circ_spot_num):
        os.remove(tempfile_path + gene_id + "_" + str(t) + ".sam")  
    return circ_reads
def proceccSamFile(item,tempfile_path,circ_chr,circ_dict_sequence,a,ciri_dict_name,divided_sam_file_num,divided_sam_file_threshold,pair_count):
    CIGAR_set2 = CIGAR_set1 = ["M","N","D"]
    circ_dict_num = np.zeros(len(ciri_dict_name))##自己去找环形isoform bsj reads
    circ_dict = dict(zip(ciri_dict_name,circ_dict_num))
    circ_dict_reads_name = [""] * len(ciri_dict_name)
    circ_dict_reads = dict(zip(ciri_dict_name,circ_dict_reads_name))
    ##SAM文件后缀
    suffix = str(item)
    if (int(suffix) != int(divided_sam_file_num)):
        pair_count_part = divided_sam_file_threshold
    else:
        pair_count_part = pair_count
    in_sam = open(tempfile_path + "selected_flag_257_" + str(suffix) + ".sam", "r")
    i = 0
    sam_total_part = 0
    circ_match_num_part = 0
    circ_name_1 = set()
    circ_name_2 = set()
    not_sure_liner_or_circ = set()
    out_sam2 = open(tempfile_path + "liner_selected_" + str(suffix) +".sam", "a")
    my_find = []
    ciri_notstar = []
    while i != pair_count_part:
        line1 = in_sam.readline()
        line1 = line1.strip('\n')
        array1 = line1.split()
        line2 = in_sam.readline()
        line2 = line2.strip('\n')
        array2 = line2.split()
        i = i + 2
        if int(array1[3]) > int(array1[7]):
            line = line1
            line1 = line2
            line2 = line
        array1 = line1.split()
        array2 = line2.split()
        flag1 = re.findall("\D+", array1[5])
        flag2 = re.findall("\D+", array2[5])
        num1 = re.findall("\d+", array1[5])
        num2 = re.findall("\d+", array2[5])
        data1 = re.findall("\d+\D", array1[5])
        data2 = re.findall("\d+\D", array2[5])
        if array1[2][0:3] != "chr":
            chr_name = "chr" + array1[2]
        else:
            chr_name = array1[2]
        if (array1[5] == "*" and array2[5] == "*"):
            continue
        if array1[2] == array2[2]:
            if chr_name in circ_chr and array1[2] == array2[2]:
                k = circ_chr.index(array2[2])
                indicator_flag1 = 0
                indicator_flag2 = 0
                indicator_flag3 = 0
                ##情况1：左端reads跨越反向剪接位点两段reads中间插入的序列跨越反向剪接位点：
                if (array1[5] == "*"):
                    last_chr = array2[2]
                    s = 0
                    for j in data2:
                        if (j[-1] == "M" or j[-1] == "N"):
                            s = s + int(j[0:-1])
                    for n in range(len(a[k])):
                        if (int(array1[3]) + s) < (int(a[k][n][3]) + 3) and (int(array1[3]) > (int(a[k][n][2]) - 3)):
                            ##取环形RNA之间的外显子序列
                            exon_sequence = circ_dict_sequence[a[k][n][0]]
                            exon_sequence_len = len(exon_sequence)
                            if (exon_sequence_len > input_read_length):
                                exon_sequence = exon_sequence[exon_sequence_len - input_read_length:] + exon_sequence[:input_read_length]
                                thread = input_read_length
                            else:
                                thread = exon_sequence_len
                                exon_sequence += exon_sequence
                            window_length = 20
                            ##取左端的前20bp
                            sam_fasta = array1[9][:window_length]
                            sam_fasta_reverse = reverse(array1[9])[:window_length]
                            sam_fasta2 = array1[9][input_read_length - window_length:]
                            sam_fasta2_reverse = reverse(array1[9])[input_read_length - window_length:]
                            max_percent = 0
                            max_percent2 = 0
                            for o in range(len(exon_sequence) - window_length):
                                window = exon_sequence[o:o + window_length]
                                percent = SequenceMatcher(None, sam_fasta, window).ratio()
                                percent_reverse = SequenceMatcher(None, sam_fasta_reverse, window).ratio()
                                percent2 = SequenceMatcher(None, sam_fasta2, window).ratio()
                                percent_reverse2 = SequenceMatcher(None, sam_fasta2_reverse, window).ratio()
                                if (percent > max_percent):
                                    if (o < thread and int(o + input_read_length) > thread):
                                        max_percent = percent
                                if (percent_reverse > max_percent):
                                    if (o < thread and int(o + input_read_length) > thread):
                                        max_percent = percent_reverse
                                if (percent2 > max_percent2):
                                    max_percent2 = percent2
                                if (percent_reverse2 > max_percent2):
                                    max_percent2 = percent_reverse2
                                if (max_percent >= 0.95 and max_percent2 >= 0.95):
                                    break
                            if max_percent >= 0.95 and max_percent2 >= 0.95:
                                indicator_flag3 = 1
                                circ_name = a[k][n][0]
                                break
                    if (indicator_flag3 == 1):
                        circ_dict[circ_name] = int(circ_dict[circ_name] + 1)
                        circ_dict_reads[circ_name] = circ_dict_reads[circ_name] + str(array1[0]) + ","
                        circ_match_num_part += 1
                    continue
                ##情况2：两段reads中间插入的序列或右端跨越反向剪接位点
                if (array2[5] == "*"):
                    last_chr = array2[2]
                    s = 0
                    for j in data1:
                        if (j[-1] == "M" or j[-1] == "N"):
                            s = s + int(j[0:-1])
                    for n in range(len(a[k])):
                        if (int(array1[3]) + s) < (int(a[k][n][3]) + 3) and (int(array1[3]) > (int(a[k][n][2]) - 3)):
                            ##取环形RNA之间的外显子序列
                            exon_sequence = circ_dict_sequence[a[k][n][0]]
                            exon_sequence_len = len(exon_sequence)
                            if (exon_sequence_len > input_read_length):
                                exon_sequence = exon_sequence[exon_sequence_len - input_read_length:] + exon_sequence[:input_read_length]
                                thread = input_read_length
                            else:
                                thread = exon_sequence_len
                                exon_sequence += exon_sequence
                            window_length = 20
                            ##取右端的前20bp
                            sam_fasta = array2[9][:window_length]
                            sam_fasta_reverse = reverse(array2[9])[:window_length]
                            ##取右端的后20bp
                            sam_fasta2 = array2[9][input_read_length - window_length:]
                            sam_fasta2_reverse = reverse(array2[9])[input_read_length - window_length:]
                            max_percent = 0
                            max_percent2 = 0
                            for o in range(len(exon_sequence) - window_length):
                                window = exon_sequence[o:o + window_length]
                                percent = SequenceMatcher(None, sam_fasta, window).ratio()
                                percent_reverse = SequenceMatcher(None, sam_fasta_reverse, window).ratio()
                                percent2 = SequenceMatcher(None, sam_fasta2, window).ratio()
                                percent_reverse2 = SequenceMatcher(None, sam_fasta2_reverse, window).ratio()
                                if (percent > max_percent):
                                    if (o < thread and int(o + input_read_length) > thread):
                                        max_percent = percent
                                if (percent_reverse > max_percent):
                                    if (o < thread and int(o + input_read_length) > thread):
                                        max_percent = percent_reverse
                                if (percent2 > max_percent2):
                                    max_percent2 = percent2
                                if (percent_reverse2 > max_percent2):
                                    max_percent2 = percent_reverse2
                                if (max_percent >= 0.95 and max_percent2 >= 0.95):
                                    break
                            if max_percent >= 0.95 and max_percent2 >= 0.95:
                                indicator_flag3 = 1
                                circ_name = a[k][n][0]
                                break
                    if (indicator_flag3 == 1):
                        circ_dict[circ_name] = int(circ_dict[circ_name] + 1)
                        circ_dict_reads[circ_name] = circ_dict_reads[circ_name] + str(array1[0]) + ","
                        circ_match_num_part += 1
                    continue
                if "S" in flag1:
                    last_chr = array2[2]
                    idx = [i for i, x in enumerate(flag1) if x == "S"]
                    index = idx[0]
                    max_S = 0
                    for id in idx:
                        if int(num1[id]) >= int(num1[index]):
                            max_S = int(num1[id])
                            index = id
                    if max_S >= 2:
                        ##左端miss的bp数
                        idx_M = [i for i, x in enumerate(flag1) if x == "M"]
                        index_M = idx_M[0]
                        sum_M = 0
                        for id_M in idx_M:
                            if sum_M <= 6:
                                sum_M = sum_M + int(num1[id_M])
                                index_M = id_M
                        mark = 0
                        s = 0  ## s记录偏离值
                        for j in data1:
                            if mark < index_M and (j[-1] == "M" or j[-1] == "N"):
                                s = s + int(j[0:-1])
                            mark = mark + 1
                        ##右端跨越的bp数
                        s2 = 0
                        for j in data2:
                            if (j[-1] == "M" or j[-1] == "N"):
                                s2 = s2 + int(j[0:-1])
                        for n in range(len(a[k])):
                            if ((abs(int(array1[3]) - int(a[k][n][2])) <= 5) or (abs(int(array1[3]) + s - int(a[k][n][2])) <= 5)) and (int(array2[3]) >= int(a[k][n][2]) and (int(array2[3]) + s2) < int(a[k][n][3])):
                                ##取环形RNA之间的外显子序列
                                exon_sequence = circ_dict_sequence[a[k][n][0]]
                                exon_sequence_len = len(exon_sequence)
                                if (exon_sequence_len > input_read_length):
                                    exon_sequence = exon_sequence[exon_sequence_len - input_read_length:] + exon_sequence[:input_read_length]
                                    thread = input_read_length
                                else:
                                    thread = exon_sequence_len
                                    exon_sequence += exon_sequence
                                window_length = 20
                                ##取左端的前20bp
                                sam_fasta = array1[9][:window_length]
                                max_percent = 0
                                for o in range(len(exon_sequence) - window_length):
                                    window = exon_sequence[o:o + window_length]
                                    percent = SequenceMatcher(None, sam_fasta, window).ratio()
                                    if (percent > max_percent):
                                        if (o < thread and int(o + input_read_length) > thread):
                                            max_percent = percent
                                    if max_percent >= 0.95:
                                        break
                                if max_percent >= 0.9:
                                    indicator_flag1 = 1
                                    circ_name = a[k][n][0]
                                    break
                    if (indicator_flag1 == 1):
                        circ_dict[circ_name] = int(circ_dict[circ_name] + 1)
                        circ_dict_reads[circ_name] = circ_dict_reads[circ_name] + str(array1[0]) + ","
                        circ_match_num_part += 1
                        continue
                if "S" in flag2:
                    last_chr = array2[2]
                    idx = [i for i, x in enumerate(flag2) if x == "S"]
                    index = idx[0]
                    max_S = 0
                    for id in idx:
                        if int(num2[id]) >= int(num2[index]):
                            max_S = int(num2[id])
                            index = id
                    if max_S >= 2:
                        ##右端miss的bp数
                        idx_M = [i for i, x in enumerate(flag2) if x == "M"]
                        index_M = idx_M[-1]
                        sum_M = 0
                        for id_M in reversed(idx_M):
                            if sum_M <= 6:
                                sum_M = sum_M + int(num2[id_M])
                                index_M = id_M
                        mark = 0
                        s = 0
                        for j in data2:
                            if mark <= index_M and (j[-1] == "M" or j[-1] == "N"):
                                s = s + int(j[0:-1])
                            mark = mark + 1
                        ##左端跨越的bp数
                        s2 = 0
                        for j in data1:
                            if (j[-1] == "M" or j[-1] == "N"):
                                s2 = s2 + int(j[0:-1])
                        for n in range(len(a[k])):
                            if abs(int(array2[3]) + s - int(a[k][n][3])) <= 5 and (int(array1[3]) >= int(a[k][n][2]) and (int(array1[3]) + s2) < int(a[k][n][3])):
                                ##取环形RNA之间的外显子序列
                                exon_sequence = circ_dict_sequence[a[k][n][0]]
                                exon_sequence_len = len(exon_sequence)
                                if (exon_sequence_len > input_read_length):
                                    exon_sequence = exon_sequence[exon_sequence_len - input_read_length:] + exon_sequence[:input_read_length]
                                    thread = input_read_length
                                else:
                                    thread = exon_sequence_len
                                    exon_sequence += exon_sequence
                                window_length = 20
                                sam_fasta = array2[9][input_read_length - window_length:]
                                max_percent = 0
                                for o in range(len(exon_sequence) - window_length):
                                    window = exon_sequence[o:o + window_length]
                                    percent = SequenceMatcher(None, sam_fasta, window).ratio()
                                    if (percent > max_percent):
                                        if ((o + window_length > thread and o < thread) or (o > thread)):
                                            max_percent = percent
                                    if max_percent >= 0.95:
                                            break
                                if max_percent >= 0.9:
                                    circ_name = a[k][n][0]
                                    indicator_flag2 = 1
                                    break
                    if (indicator_flag2 == 1):
                        circ_dict[circ_name] = int(circ_dict[circ_name] + 1)
                        circ_dict_reads[circ_name] = circ_dict_reads[circ_name] + str(array1[0]) + ","
                        circ_match_num_part += 1
                        continue
            if (int(array1[1]) == 145 and int(array2[1]) == 97) or (int(array1[1]) == 81 and int(array2[1]) == 161):
                CIGAR = re.findall("\D+", array1[5]) + re.findall("\D+", array2[5])
                CIGAR = list(set(CIGAR))
                flag = set(CIGAR).issubset(set(CIGAR_set))
                if flag == True:
                    sam_total_part += 1
                    t11 = int(array1[3])
                    data1 = re.findall("\d+\D", array1[5])
                    s = 0
                    for j in data1:
                        if j[-1] in CIGAR_set2:
                            s = s + int(j[0:-1])
                    t12 = t11 + s - 1
                    t21 = int(array2[3])
                    data2 = re.findall("\d+\D", array2[5])
                    s = 0
                    for j in data2:
                        if j[-1] in CIGAR_set2:
                            s = s + int(j[0:-1])
                    t22 = t21 + s - 1
                    if t12 < t22 and array2[2] in circ_chr:
                        if (int(array2[3]) - int(array1[3]) > 2):
                            circ_name_2.add(array2[0])
                            line111 = convert_flag(array1)
                            line222 = convert_flag(array2)
                            out_samB = open(tempfile_path + "selected_" + chr_name + "_B" + str(suffix) + ".sam", "a")
                            out_samB.write(line111 + "\n")
                            out_samB.write(line222 + "\n")
                            out_samB.close()
                            continue
            if (int(array1[1]) == 99 and int(array2[1]) == 147) or (int(array1[1]) == 163 and int(array2[1]) == 83):
                CIGAR = re.findall("\D+", array1[5]) + re.findall("\D+", array2[5])
                CIGAR = list(set(CIGAR))
                flag = set(CIGAR).issubset(set(CIGAR_set))
                if flag == True:
                    sam_total_part += 1
                    t11 = int(array1[3])
                    data1 = re.findall("\d+\D", array1[5])
                    s = 0
                    for j in data1:
                        if j[-1] in CIGAR_set2:
                            s = s + int(j[0:-1])
                    t12 = t11 + s - 1
                    t21 = int(array2[3])
                    data2 = re.findall("\d+\D", array2[5])
                    s = 0
                    for j in data2:
                        if j[-1] in CIGAR_set2:
                            s = s + int(j[0:-1])
                    t22 = t21 + s - 1
                    if t12 < t22 and array2[2] in circ_chr:
                        out_samA = open(tempfile_path + "selected_" + chr_name + "_A" + str(suffix) + ".sam", "a")
                        out_samA.write(line1 + "\n")
                        out_samA.write(line2 + "\n")
                        out_samA.close()
                        k = circ_chr.index(array2[2])
                        for n in range(len(a[k])):
                            if (int(array1[3]) >= int(a[k][n][2]) - 2) and (t22 <= int(a[k][n][3]) + 2):
                                not_sure_liner_or_circ.add(array1[0])
                                out_sam5 = open(tempfile_path + "selected_circAST" + chr_name + "_A" + str(suffix) + ".sam", "a")
                                out_sam5.write(line1 + "\n")
                                out_sam5.write(line2 + "\n")
                                out_sam5.close()
                                out_sam2.write(line1 + "\n")
                                out_sam2.write(line2 + "\n")
                                break
            if (array1[0] not in circ_name_2) and (array1[0] not in not_sure_liner_or_circ):
                out_sam2.write(line1 + "\n")
                out_sam2.write(line2 + "\n")
    in_sam.close()
    out_sam2.close()
    return circ_dict,circ_dict_reads,sam_total_part,circ_match_num_part


if __name__ == "__main__":
    parse = optparse.OptionParser()
    parse.add_option('-G','--gtf',dest = 'gtf',action = 'store',metavar = 'gtf file name',help = 'please enter your annotation gtf file')
    parse.add_option('-I','--file',dest = 'file',action = 'store',metavar = 'input files',help = 'enter a SAM alignment of pair-end reads which is generated by tophat2 and should be sorted by name')
    parse.add_option('-F','--fasta',dest = 'fasta',action = 'store',metavar = 'input fasta file',help = 'please enter your fasta file path')
    parse.add_option('-J','--junction',dest = 'junction',action = 'store',metavar = 'circRNA junction list',help = 'please enter a circRNA junction list generated by CIRCexplorer2, CIRI and UROBORUS or  or other tools in a given format')
    parse.add_option('-T','--threshold',dest = 'threshold',action = 'store',metavar = 'define a threshold of reads number supporting back-splicing junction',help = 'please define a threshold of reads number supporting back-splicing junction')
    parse.add_option('-L','--length',dest = 'length',action = 'store',metavar = 'reads length',help = 'please enter the length of reads')
    parse.add_option('-S','--StringtieDir',dest = 'StringtieDir',action = 'store',metavar = 'StringtieDir',help = 'please enter your Stringtie location path')
    isTotal = 1
    (options,args) = parse.parse_args()
    #check input files
    for file in ([options.length,options.gtf,options.file,options.junction]):
        if not (file):
            print(sys.stderr, "\nError: Lack of input file!\n")
            parse.print_help()
            sys.exit(0)

    file_sorted = sorted_or_not(options.file)
    if file_sorted == 0:
        print(sys.stderr, "\nError: Input sam file must be pair-end reads and sorted by name!\n")
        parse.print_help()
        sys.exit(0)
    else:
        inPutFileName = options.file
    inputGtfName = options.gtf
    inputCircListNameCircAST2 = options.junction
    inputCircListName = options.junction
    input_read_length = int(options.length)
    inputfrag_len = input_read_length * 2 + 100
    StringTie_dir = options.StringtieDir
    input_fasta_name = options.fasta
    divided_sam_file_threshold = 20000000 ## 500000行SAM文件对应160M容量 6000000w大概2.5G  20000000
    MAX_WORKERS = 16


    for file in ([options.threshold]):
        if not (file):
            thresh = 1
        else:
            thresh = int(options.threshold)


    if os.path.exists("CircAST_temp/") == False:
        os.mkdir("CircAST_temp")
    tempfile_path = "CircAST_temp/"
    input_stringtie_gtf_name = tempfile_path + "stringtie_filtered.gtf"


    CIGAR_set = ["M","I","N","D","S"]
    CIGAR_set1 = ["M","N","D","S"]
    CIGAR_set2 = ["M","N","D"]


    ## Divide the GTF file by chromosome
    chr_name_gtf = []
    in_gtf = open(inputGtfName)  # split gtf files by chromosomes
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        if len(array[0]) < 6:
            if array[2] == "exon":
                chr_name_temp = array[0]
                out_gtf = open(tempfile_path + chr_name_temp + ".gtf", "a")
                out_gtf.write(line + "\n")
                out_gtf.close()
                if chr_name_temp not in chr_name_gtf:
                    chr_name_gtf.append(chr_name_temp)
    in_gtf.close()

    ## inputCircListNameCircAST2
    left = [None] * len(chr_name_gtf)
    right = [None] * len(chr_name_gtf)
    for i in range(0,len(chr_name_gtf)):
        in_gtf = open(tempfile_path + chr_name_gtf[i] + ".gtf","r")
        left[i] = []
        right[i] = []
        for line in in_gtf:
            line = line.strip('\n')
            array = line.split()
            if array[3] not in left[i]:
                left[i].append(array[3])
            if array[4] not in right[i]:
                right[i].append(array[4])
        in_gtf.close()

    list_class = 0
    # check inputCircListNameCircAST2
    in_txt = open(inputCircListNameCircAST2)
    line = in_txt.readline()
    line = in_txt.readline()
    line = line.strip('\n')
    array = line.split()
    in_txt.close()
    if len(array) > 9:
        if array[10] == "+" or array[10] == "-":  # CIRI's result
            list_class = 2
    if array[5] == "+" or array[5] == "-":  # circexplore's result
        list_class = 1
    if array[3] == "+" or array[3] == "-":  # UB's result
        list_class = 3
    if list_class == 0:
        print(sys.stderr, "\nError: Input circRNA list file must be from CIRCexplore2,CIRI OR UROBORUS!\n")
        parse.print_help()
        sys.exit(0)
    BSJ_left = 0
    BSJ_right = 0
    BSJ_name = 0

    ##########################################################################这里的最后一列circRNA reads数需要修改
    if list_class == 1:  # circexplore's result
        in_list = open(inputCircListNameCircAST2)
        for line in in_list:
            line = line.strip('\n')
            array = line.split()
            if array[0] in chr_name_gtf:
                k = chr_name_gtf.index(array[0])
                start = str(int(array[1]) + 1)
                end = array[2]
                temp = array[3].split("/")
                if start in left[k] and end in right[k]:
                    gene_name_L = findname1(array[0], start)
                    gene_name_R = findname2(array[0], end)
                    gene_name_same = list(set(gene_name_L).intersection(set(gene_name_R)))
                    if gene_name_same != []:
                        gene_name_same.sort()
                        out_list = open(tempfile_path + "circRNA_list_circAST2.txt", "a")
                        out_list.write(temp[0] + "\t" + array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[5] + "\t" + gene_name_same[0] + "\t" + str(int(array[2]) - int(array[1])) + "\t" + temp[1] + "\n")
                        out_list.close()
        inputCircListNameCircAST2 = tempfile_path + "circRNA_list_circAST2.txt"
        in_list.close()
    #################################################################################将匹配到的circRNA_junction文件中可以在gtf文件中找到开头结尾并且开头结尾匹配到同一基因的写出，为什么只取gene_name_same[0]这一个？？？？？？
    if list_class == 2:  # CIRI's result
        in_list = open(inputCircListNameCircAST2)
        in_list.readline()
        for line in in_list:
            line = line.strip('\n')
            array = line.split()
            if array[1][0:3] != "chr":
                array[1] = "chr" + array[1]
            if array[1] in chr_name_gtf:
                if array[8] == "exon":
                    k = chr_name_gtf.index(array[1])
                    start = array[2]
                    end = array[3]
                    if start in left[k] and end in right[k]:
                        gene_name_L = findname1(array[1], start)
                        gene_name_R = findname2(array[1], end)
                        gene_name_same = list(set(gene_name_L).intersection(set(gene_name_R)))
                        if gene_name_same != []:
                            gene_name_same.sort()
                            out_list = open(tempfile_path + "circRNA_list_circAST2.txt", "a")
                            out_list.write(array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[3] + "\t" + array[10] + "\t" + gene_name_same[0] + "\t" + str(int(array[3]) - int(array[2]) + 1) + "\t" + array[4] + "\n")
                            # out_list.write(
                                # array[1] + "\t" + str(int(array[2]) - 1) + "\t" + array[3] + "\t" + array[10] + "\t" +
                                # gene_name_same[0] + "\t" + str(int(array[3]) - int(array[2]) + 1) + "\t" + str(array[4]) + "\n")
                            out_list.close()
        inputCircListNameCircAST2 = tempfile_path + "circRNA_list_circAST2.txt"
        in_list.close()

    if list_class == 3:  # UB's result
        in_list = open(inputCircListNameCircAST2)
        for line in in_list:
            line = line.strip('\n')
            array = line.split()
            if array[0] in chr_name_gtf:
                k = chr_name_gtf.index(array[0])
                start = str(int(array[1]) + 1)
                end = array[2]
                if start in left[k] and end in right[k]:
                    gene_name_L = findname1(array[0], start)
                    gene_name_R = findname2(array[0], end)
                    gene_name_same = list(set(gene_name_L).intersection(set(gene_name_R)))
                    if gene_name_same != []:
                        gene_name_same.sort()
                        out_list = open(tempfile_path + "circRNA_list_circAST2.txt", "a")
                        out_list.write(str(array[0] + array[1] + array[2]) + "\t" + array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[3] + "\t" + gene_name_same[0] + "\t" + array[5] + "\t" + array[6] + "\n")
                        out_list.close()
        inputCircListNameCircAST2 = tempfile_path + "circRNA_list_circAST2.txt"
        in_list.close()


    circ_name_ciri = set()
    circ_chr = []
    ciri_dict_name = []
    in_txt = open(inputCircListNameCircAST2)
    for line in in_txt:
        line = line.strip('\n')
        array = line.split()
        array = line.split()
        if array[1][0:3] != "chr":
            array[1] = "chr" + array[1]
        if array[1] in chr_name_gtf:
            out_txt = open(tempfile_path + "ciri_" + array[1] + ".txt", "a")
            out_txt.write(line + "\n")
            out_txt.close()
            ciri_dict_name.append(array[0])
            # for j in range(int(array[4])):
                # circ_name_ciri.add(array[-1][j])
            if array[1] not in circ_chr:
                circ_chr.append(array[1])
    in_txt.close()
    ## sort circ TXT file to find bsj reads 
    ciri_sort = [None] * len(circ_chr)
    for i in range(0, len(circ_chr)):
        chr_name = circ_chr[i]
        m = 0
        in_ciri = open(tempfile_path + "ciri_" + chr_name + ".txt", "r")
        ciri_sort[i] = []
        for line in in_ciri:
            line = line.strip('\n')
            array = line.split()
            ciri_sort[i].append(array)
        in_ciri.close()
        ciri_sort[i].sort(key=lambda x: (int(x[3])-int(x[2]),int(x[3])+int(x[2])))
        out_ciri = open(tempfile_path + "ciri_" + chr_name + ".txt", "w")
        for j in range(len(ciri_sort[i])):
            out_ciri.write("\t".join(ciri_sort[i][j]) + "\n")
        out_ciri.close()

    circ_dict_num = np.zeros(len(ciri_dict_name))##自己去找环形isoform bsj reads
    circ_dict = dict(zip(ciri_dict_name,circ_dict_num))
    circ_dict_reads_name = [""] * len(ciri_dict_name)
    circ_dict_reads = dict(zip(ciri_dict_name,circ_dict_reads_name))
    circ_dict_sequence_list = [""] * len(ciri_dict_name)
    circ_dict_sequence = dict(zip(ciri_dict_name,circ_dict_sequence_list))



    ## Fasta information was recorded by chromosome
    fasta = [None] * len(circ_chr)
    len_fasta = 0
    for i in range(0, len(circ_chr)):
        chr_name = circ_chr[i]
        m = 0
        in_fasta = open(input_fasta_name + circ_chr[i] + ".fa", "r")
        fasta[i] = []
        for line in in_fasta:
            m = m + 1
            if m > 1:
                line = line.strip('\n')
                if len_fasta == 0:
                    len_fasta = len(line)
                fasta[i].append(line)
        in_fasta.close()

    ## CircRNA information was recorded by chromosome
    a = [None] * len(circ_chr)
    A = []
    for i in range(len(circ_chr)):
        a[i] = []
        chr_name = circ_chr[i]
        in_txt = open(tempfile_path + "ciri_" + chr_name + ".txt", "r")
        in_txt_read = in_txt.read()
        b = in_txt_read.split('\n')
        in_txt.close()
        for j in range(len(b) - 1):
            b[j] = b[j].split()
            a[i].append(b[j])
            A.append(b[j])

    # ##尝试把find_sequence放在外面
    # for n in range(len(A)):
        # exon_sequence = findExonSequence(A[n][1], circ_chr, fasta, int(A[n][2]), int(A[n][3]))
        # circ_dict_sequence[A[n][0]] = exon_sequence
    ##多线程findExonSequence
    start_time = time.time()
    partial_func_findExonSequence = partial(findExonSequence, circ_chr = circ_chr,fasta = fasta)
    with ThreadPoolExecutor() as pool:
        results = pool.map(partial_func_findExonSequence, A)
        resultsList = list(zip(A, results))
        for item, result in resultsList:
            circ_dict_sequence[item[0]] = result
    end_time = time.time()
    print("找环形RNA序列多线程运行成功，时间为" + str(end_time - start_time))

    ##----------------------------------------------------------------------------------------------------------------------
    ## Find the first line of the SAM file
    in_sam = open(inPutFileName)
    out_sam2 = open(tempfile_path + "liner_selected.sam", "a")
    i = 0
    for line in in_sam:
        line = line.strip('\n')
        array = line.split()
        if array[0][0] == '@':
            i = i + 1
            out_sam2.write(line + "\n")
        else:
            head_num = i
            break
    in_sam.close()
    out_sam2.close()
    # ##去掉最后一行\n
    # f = open(tempfile_path + "liner_selected.sam", 'rb+')
    # f.seek(-1, os.SEEK_END)
    # f.truncate()
    # f.close()

    for i in range(0, len(circ_chr)):
        chr_name = circ_chr[i]
        in_sam = open(inPutFileName)
        out_sam1 = open(tempfile_path + "selected_" + chr_name + "_A.sam", "a")
        out_sam2 = open(tempfile_path + "selected_" + chr_name + "_B.sam", "a")
        out_sam5 = open(tempfile_path + "selected_circAST" + chr_name + "_A.sam", "a")
        out_sam6 = open(tempfile_path + "selected_circAST" + chr_name + "_B.sam", "a")
        j = 0
        while j < head_num:
            line = in_sam.readline()
            line = line.strip('\n')
            out_sam1.write(line + "\n")
            out_sam2.write(line + "\n")
            out_sam5.write(line + "\n")
            out_sam6.write(line + "\n")
            j += 1
        in_sam.close()
        out_sam1.close()
        out_sam2.close()
        out_sam5.close()
        out_sam6.close()
        # f = open(tempfile_path + "selected_" + chr_name + "_A.sam", 'rb+')
        # f.seek(-1, os.SEEK_END)
        # f.truncate()
        # f.close()
        # f = open(tempfile_path + "selected_" + chr_name + "_B.sam", 'rb+')
        # f.seek(-1, os.SEEK_END)
        # f.truncate()
        # f.close()
        # f = open(tempfile_path + "selected_circAST" + chr_name + "_A.sam", 'rb+')
        # f.seek(-1, os.SEEK_END)
        # f.truncate()
        # f.close()
        # f = open(tempfile_path + "selected_circAST" + chr_name + "_B.sam", 'rb+')
        # f.seek(-1, os.SEEK_END)
        # f.truncate()
        # f.close()

    ## Select primary alignment reads
    i = 0
    pair_count = 0
    chr_name_sam = []
    in_sam = open(inPutFileName)
    out_sam_divide = 0
    divided_sam_file_num = 0
    divided_sam_file_num_list = []
    divided_sam_file_num_list.append(divided_sam_file_num)
    out_sam = open(tempfile_path + "selected_flag_257_" + str(divided_sam_file_num) + ".sam", "w")
    for line in in_sam:
        i = i + 1
        line = line.strip('\n')
        if i >= head_num + 1:
            array = line.split()
            flag = int(array[1])
            if flag < 257 and out_sam_divide < divided_sam_file_threshold:
                out_sam.write(line + "\n")
                pair_count += 1
                out_sam_divide += 1
                continue
            if flag < 257 and out_sam_divide >= divided_sam_file_threshold:
                out_sam.close()
                pair_count = 0
                out_sam_divide = 0
                divided_sam_file_num += 1
                divided_sam_file_num_list.append(divided_sam_file_num)
                out_sam = open(tempfile_path + "selected_flag_257_" + str(divided_sam_file_num) + ".sam", "w")
                out_sam.write(line + "\n")
                out_sam_divide += 1
                pair_count += 1
    out_sam.close()
    in_sam.close()

    print("开始proceccSamFile")
    sam_total = 0
    circ_match_num = 0
    start_time = time.time()
    partial_func_proceccSamFile = partial(proceccSamFile, tempfile_path=tempfile_path, circ_chr=circ_chr, circ_dict_sequence=circ_dict_sequence, a=a, ciri_dict_name=ciri_dict_name, divided_sam_file_num=divided_sam_file_num, divided_sam_file_threshold=divided_sam_file_threshold, pair_count=pair_count)
    with ProcessPoolExecutor() as pool:
        results = pool.map(partial_func_proceccSamFile, divided_sam_file_num_list)
        resultsList = list(zip(divided_sam_file_num_list, results))
        for item, result in resultsList:
            print(str(item) + "\t" + str(result))
            for key, value in result[0].items():
                if key in circ_dict:
                    circ_dict[key] += int(float(value))
            for key, value in result[1].items():
                if key in circ_dict_reads:
                    circ_dict_reads[key] += value
            sam_total += int(result[2])
            circ_match_num += int(result[3])
    end_time = time.time()
    ### return circ_dict,circ_dict_reads,sam_total_part,circ_match_num_part
    print("proceccSamFile多线程运行成功，时间为" + str(end_time - start_time))
    ## shutil.copyfile(tempfile_path + gene_id + "_A.sam", tempfile_path + gene_id + ".sam")
    for i in range(divided_sam_file_num + 1):
        ##Removes the last "\n" of the file
        # f = open(tempfile_path + "liner_selected_" + str(i) +".sam", 'rb+')
        # f.seek(-1, os.SEEK_END)
        # f.truncate()
        # f.close()
        os.system("cat " + tempfile_path + "liner_selected_" + str(i) +".sam >> "+ tempfile_path + "liner_selected.sam")
        os.remove(tempfile_path + "liner_selected_" + str(i) +".sam")
        for j in range(0, len(circ_chr)):
            chr_name = circ_chr[j]
            if(os.path.exists(tempfile_path + "selected_" + chr_name + "_A" + str(i) + ".sam")):
                # f = open(tempfile_path + "selected_" + chr_name + "_A" + str(i) + ".sam", 'rb+')
                # f.seek(-1, os.SEEK_END)
                # f.truncate()
                # f.close()
                os.system("cat " + tempfile_path + "selected_" + chr_name + "_A" + str(i) + ".sam >> " + tempfile_path + "selected_" + chr_name + "_A.sam")
                os.remove(tempfile_path + "selected_" + chr_name + "_A" + str(i) + ".sam")
            if (os.path.exists(tempfile_path + "selected_" + chr_name + "_B" + str(i) + ".sam")):
                # f = open(tempfile_path + "selected_" + chr_name + "_B" + str(i) + ".sam", 'rb+')
                # f.seek(-1, os.SEEK_END)
                # f.truncate()
                # f.close()
                os.system("cat " + tempfile_path + "selected_" + chr_name + "_B" + str(i) + ".sam >> " + tempfile_path + "selected_" + chr_name + "_B.sam")
                os.remove(tempfile_path + "selected_" + chr_name + "_B" + str(i) + ".sam")
            if (os.path.exists(tempfile_path + "selected_circAST" + chr_name + "_A" + str(i) + ".sam")):
                # f = open(tempfile_path + "selected_circAST" + chr_name + "_A" + str(i) + ".sam", 'rb+')
                # f.seek(-1, os.SEEK_END)
                # f.truncate()
                # f.close()
                os.system("cat " + tempfile_path + "selected_circAST" + chr_name + "_A" + str(i) + ".sam >> " + tempfile_path + "selected_circAST" + chr_name + "_A.sam")
                os.remove(tempfile_path + "selected_circAST" + chr_name + "_A" + str(i) + ".sam")
            if (os.path.exists(tempfile_path + "selected_circAST" + chr_name + "_B" + str(i) + ".sam")):
                # f = open(tempfile_path + "selected_circAST" + chr_name + "_B" + str(i) + ".sam", 'rb+')
                # f.seek(-1, os.SEEK_END)
                # f.truncate()
                # f.close()
                os.system("cat " + tempfile_path + "selected_circAST" + chr_name + "_B" + str(i) + ".sam >> " + tempfile_path + "selected_circAST" + chr_name + "_B.sam")
                os.remove(tempfile_path + "selected_circAST" + chr_name + "_B" + str(i) + ".sam")


    # ##----------------------------------------------------------------------------------------------------------------------
    out_BSJ_reads = open(tempfile_path + "out_bsj_reads.txt", "w")
    for i in range(0, len(ciri_dict_name)):
        out_BSJ_reads.write(str(ciri_dict_name[i]) + "\t" + circ_dict_reads[ciri_dict_name[i]].strip(",") + "\n")
    out_BSJ_reads.close()



    ##call stringtie
    os.system('samtools sort -o ' + tempfile_path +'liner_selected.bam ' + tempfile_path +'liner_selected.sam')
    os.system(StringTie_dir + ' -p 8 -G ' + inputGtfName +' -e -o ' +tempfile_path +'test.gtf -l test ' + tempfile_path +'liner_selected.bam')
    # os.system(StringTie_dir + ' -p 8 -G ' + inputGtfName +' -o ' +tempfile_path +'test.gtf -l test ' + tempfile_path +'liner_selected.bam')

    ##filter low expression transcripts
    in_gtf = open(tempfile_path +'test.gtf')
    out_gtf = open(tempfile_path + "stringtie_filtered.gtf","w")
    canWrite = 0
    i = 0
    for line in in_gtf:
        i += 1
        if i <= 2:
            continue
        line = line.strip('\n')
        array = line.split('\t')
        if len(array[0]) < 6:
            if array[2] == "transcript" and float(array[8].split('"')[-6]) > 10:
                out_gtf.write(line + "\n")
                canWrite = 1
            if array[2] == "transcript" and float(array[8].split('"')[-6]) < 10:
                canWrite = 0
            if array[2] == "exon" and canWrite == 1:
                out_gtf.write(line + "\n")
    in_gtf.close()
    out_gtf.close()

    # # Call CircAST
    # CircAST(inputGtfName,tempfile_path + "circAST_input_pe.sam",inputCircListNameCircAST2,1,input_read_length)

    # #Split circAST files by chromosomes
    # in_txt = open(input_circAST_txt_name)
    # for line in in_txt:
        # line = line.strip('\n')
        # array = line.split()
        # if len(array[0]) < 6:
            # chr_name = array[0]
            # out_gtf = open(tempfile_path + chr_name + "_circAST.txt","a")
            # out_gtf.write(line + "\n")
            # out_gtf.close()
    # in_txt.close()

    #Split stringtie files by chromosomes
    chr_name_stringtie_gtf = []
    in_gtf = open(input_stringtie_gtf_name)
    for line in in_gtf:
        line = line.strip('\n')
        array = line.split()
        if len(array[0]) < 6:
            if array[2] == "exon" or array[2] == "transcript":
                chr_name = array[0]
                out_gtf = open(tempfile_path + chr_name + "_stringtie.gtf","a")
                out_gtf.write(line + "\n")
                out_gtf.close()
                if chr_name not in chr_name_stringtie_gtf:
                    chr_name_stringtie_gtf.append(chr_name)
    in_gtf.close()


    ##########################################################################这里的最后一列circRNA reads数需要修改
    if list_class == 1:  # circexplore's result
        in_list = open(inputCircListName)
        for line in in_list:
            line = line.strip('\n')
            array = line.split()
            if array[0] in chr_name_gtf:
                k = chr_name_gtf.index(array[0])
                start = str(int(array[1]) + 1)
                end = array[2]
                temp = array[3].split("/")
                if start in left[k] and end in right[k]:
                    gene_name_L = findname1(array[0], start)
                    gene_name_R = findname2(array[0], end)
                    gene_name_same = list(set(gene_name_L).intersection(set(gene_name_R)))
                    if gene_name_same != [] and isTotal == 1:
                        gene_name_same.sort()
                        out_list = open(tempfile_path + "circRNA_list.txt", "a")
                        if int(float(circ_dict[temp[0]])) != 0:
                            out_list.write(
                                array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[5] + "\t" + gene_name_same[
                                    0] + "\t" + str(int(array[2]) - int(array[1])) + "\t" + str(circ_dict[temp[0]]) + "\n")
                        if int(float(circ_dict[temp[0]])) == 0:
                                                    out_list.write(
                            array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[5] + "\t" + gene_name_same[
                                0] + "\t" + str(int(array[2]) - int(array[1])) + "\t" + temp[1] + "\n")
                        out_list.close()
                    if gene_name_same != [] and isTotal == 0:
                        gene_name_same.sort()
                        out_list = open(tempfile_path + "circRNA_list.txt", "a")
                        out_list.write(
                            array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[5] + "\t" + gene_name_same[
                                0] + "\t" + str(int(array[2]) - int(array[1])) + "\t" + temp[1] + "\n")
                        out_list.close()
        inputCircListName = tempfile_path + "circRNA_list.txt"
        in_list.close()
    #################################################################################将匹配到的circRNA_junction文件中可以在gtf文件中找到开头结尾并且开头结尾匹配到同一基因的写出，为什么只取gene_name_same[0]这一个？？？？？？
    if list_class == 2:  # CIRI's result
        in_list = open(inputCircListName)
        in_list.readline()
        for line in in_list:
            line = line.strip('\n')
            array = line.split()
            if array[1][0:3] != "chr":
                array[1] = "chr" + array[1]
            if array[1] in chr_name_gtf:
                if array[8] == "exon":
                    k = chr_name_gtf.index(array[1])
                    start = array[2]
                    end = array[3]
                    if start in left[k] and end in right[k]:
                        gene_name_L = findname1(array[1], start)
                        gene_name_R = findname2(array[1], end)
                        gene_name_same = list(set(gene_name_L).intersection(set(gene_name_R)))
                        if gene_name_same != [] and isTotal == 1:
                            gene_name_same.sort()
                            out_list = open(tempfile_path + "circRNA_list.txt", "a")
                            if int(float(circ_dict[array[0]])) != 0:
                                out_list.write(
                                    array[1] + "\t" + str(int(array[2]) - 1) + "\t" + array[3] + "\t" + array[10] + "\t" +
                                    gene_name_same[0] + "\t" + str(int(array[3]) - int(array[2]) + 1) + "\t" + str(circ_dict[array[0]]) + "\n")
                            if int(float(circ_dict[array[0]])) == 0:
                                out_list.write(
                                    array[1] + "\t" + str(int(array[2]) - 1) + "\t" + array[3] + "\t" + array[10] + "\t" +
                                    gene_name_same[0] + "\t" + str(int(array[3]) - int(array[2]) + 1) + "\t" + str(array[4]) + "\n")
                                out_list.close()
                        if gene_name_same != [] and isTotal == 0:
                            gene_name_same.sort()
                            out_list = open(tempfile_path + "circRNA_list.txt", "a")
                            out_list.write(
                                array[1] + "\t" + str(int(array[2]) - 1) + "\t" + array[3] + "\t" + array[10] + "\t" +
                                gene_name_same[0] + "\t" + str(int(array[3]) - int(array[2]) + 1) + "\t" + array[4] + "\n")
                            # out_list.write(
                                # array[1] + "\t" + str(int(array[2]) - 1) + "\t" + array[3] + "\t" + array[10] + "\t" +
                                # gene_name_same[0] + "\t" + str(int(array[3]) - int(array[2]) + 1) + "\t" + str(array[4]) + "\n")
                            out_list.close()
        inputCircListName = tempfile_path + "circRNA_list.txt"
        in_list.close()

    if list_class == 3:  # UB's result
        in_list = open(inputCircListName)
        for line in in_list:
            line = line.strip('\n')
            array = line.split()
            if array[0] in chr_name_gtf:
                k = chr_name_gtf.index(array[0])
                start = str(int(array[1]) + 1)
                end = array[2]
                if start in left[k] and end in right[k]:
                    gene_name_L = findname1(array[0], start)
                    gene_name_R = findname2(array[0], end)
                    gene_name_same = list(set(gene_name_L).intersection(set(gene_name_R)))
                    if gene_name_same != [] and isTotal == 1:
                        gene_name_same.sort()
                        out_list = open(tempfile_path + "circRNA_list.txt", "a")
                        if int(float(circ_dict[str(array[0] + array[1] + array[2])])) != 0:
                            out_list.write(
                                array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[3] + "\t" + gene_name_same[
                                    0] + "\t" + array[5] + "\t" + str(circ_dict[str(array[0] + array[1] + array[2])]) + "\n")
                        if int(float(circ_dict[str(array[0] + array[1] + array[2])])) == 0:
                            out_list.write(
                            array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[3] + "\t" + gene_name_same[
                                0] + "\t" + array[5] + "\t" + array[6] + "\n")
                        out_list.close()
                    if gene_name_same != [] and isTotal == 0:
                        gene_name_same.sort()
                        out_list = open(tempfile_path + "circRNA_list.txt", "a")
                        out_list.write(
                            array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + array[3] + "\t" + gene_name_same[
                                0] + "\t" + array[5] + "\t" + array[6] + "\n")
                        out_list.close()
        inputCircListName = tempfile_path + "circRNA_list.txt"
        in_list.close()
    ## get gene_name and chr_no for future use
    inputCircListName = tempfile_path + "circRNA_list.txt"
    chr_set = []
    gene_name = []
    chr_no = []
    in_txt = open(inputCircListName)
    out_txt = open(tempfile_path + "circRNA_list_selected.txt","w")
    for line in in_txt:
        line = line.strip('\n')
        array = line.split()
        if int(float(array[6])) >= thresh:
            out_txt.write(line+ "\n")
            if array[0] not in chr_set:
                chr_set.append(array[0])
            if array[4] not in gene_name:
                gene_name.append(array[4])
                chr_no.append(array[0])
            else:######################################这步是考虑到一个基因名对应到多个染色体的情况
                chr_temp = []
                for i in range(0,len(gene_name)):
                    if gene_name[i] == array[4]:
                        chr_temp.append(chr_no[i])
                if array[0] not in chr_temp:
                    gene_name.append(array[4])
                    chr_no.append(array[0])
    in_txt.close()
    out_txt.close()

    match_num_total = sam_total + circ_match_num


    if (isTotal == 1):
        circ_reads_name = []
        partial_func_calculate_mixed = partial(calculate_mixed, gene_name = gene_name,chr_no = chr_no,match_num_total = match_num_total,head_num = head_num,inputfrag_len = inputfrag_len,tempfile_path = tempfile_path)
        with ProcessPoolExecutor() as pool:
            results = pool.map(partial_func_calculate_mixed, gene_name)
            resultsList = list(zip(gene_name, results))
            for gene, circ_reads in resultsList:
                if (circ_reads != 0):
                    circ_reads_name = circ_reads_name + circ_reads

    if (isTotal == 0):
        partial_func_calculate = partial(calculate, gene_name = gene_name,chr_no = chr_no,match_num_total = match_num_total,head_num = head_num,inputfrag_len = inputfrag_len,tempfile_path = tempfile_path)
        with ProcessPoolExecutor() as pool:
            pool.map(partial_func_calculate, gene_name)
    gene = []
    BSJ = []
    BSJ_iso = []
    gene_BSJ = []
    FPKM_total = 0
    in_table = open("result_table.txt", "r")
    for line in in_table:
        line = line.strip('\n')
        array = line.split("\t")
        circ = array[0] + "*" + array[1] + "*" + array[2] + "*" + array[3]
        if array[3] not in gene:
            gene.append(array[3])
        if circ not in BSJ:
            BSJ.append(circ)
            gene_BSJ.append(array[3])
        BSJ_iso.append(circ)
        FPKM_total = FPKM_total + float(array[9])
    in_table.close()

    BSJ_num = [0] * len(gene)
    for i in range(0, len(gene)):
        BSJ_num[i] = gene_BSJ.count(gene[i])
    iso_num = [0] * len(BSJ)
    for i in range(0, len(BSJ)):
        iso_num[i] = BSJ_iso.count(BSJ[i])

    iso_name = []
    for i in gene:
        for j in range(0, BSJ_num[0]):
            for k in range(0, iso_num[0]):
                iso_name.append("circ" + i + "-" + str(j + 1) + "-" + str(k + 1))
            iso_num.pop(0)
        BSJ_num.pop(0)

    in_table = open("result_table.txt", "r")
    out_table = open("result_circ.txt", "w")
    for line in in_table:
        line = line.strip('\n')
        array = line.split("\t")
        TPM = float(array[9]) * 10 ** 6 / (FPKM_total * 1.0)
        TPM = float('%.4f' % TPM)
        out_table.write(array[0] + "\t" + array[1] + "\t" + array[2] + "\t" + iso_name[0] + "\t" + array[3] + "\t" + array[4] + "\t" + array[5] + "\t" + array[6] + "\t" + array[7] + "\t" + array[8] + "\t" + array[9] + "\t" + str(TPM)  + "\n")
        iso_name.pop(0)
    in_table.close()
    out_table.close()


    os.remove("result_table.txt")