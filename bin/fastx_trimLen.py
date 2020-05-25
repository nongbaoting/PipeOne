#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/12 21:46
# @Author  : Baoting Nong
# @Email   : 523135753@qq.com
import fire,subprocess,os
from Bio import SeqIO
from statistics import median
def get_len(fq):
    mylen = []
    count = 0
    for record in SeqIO.parse(fq,'fastq'):
        mylen.append(len(record))
        count +=1
        if count >100:
            break
    this_len = round(median( mylen)) - 10
    print(this_len )
    return 0

if __name__ == '__main__':
    fire.Fire(get_len)