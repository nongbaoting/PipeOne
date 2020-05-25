#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/9/29 18:58
# @Author  : Baoting Nong
import fire
from collections import defaultdict

class Venn:

    def venn(self, file1, file2,out,out_rev, col1=0,col2=0,header1=0,header2=0):
        lines1 = open(file1,'r').read().strip().split('\n')
        lines2 = open(file2, 'r').read().strip().split('\n')

        fo = open(out, 'w')
        fo_rev = open(out_rev, 'w')
        start1, start2 = 0,0
        if header1:
            fo.write(lines1[0] + '\n')
            fo_rev.write(lines1[0] + '\n')
            start1 = 1
        if header2:
            start2 = 1
        myd =[]
        for li in lines2[start2 :]:
            cell = li.split('\t')
            myd.append(cell[col2] )
        myd =set(myd)
        for l in lines1[start1 :]:
            cl = l.split('\t')

            if cl[col1].strip() in myd:
                fo.write( l + '\n')
            else:
                fo_rev.write(l + '\n')
        fo.close()

    def circ_bed(self,bed,circ,out):
        circlist = open(circ, 'r').read().strip().split('\n')
        circlist = set(circlist)
        with open(bed, 'r') as fb, open(out, 'w') as fo:
            for line in fb:
                cell = line.split('\t')
                circid = f"{cell[0]}:{cell[1]}|{cell[2]}"
                if circid in circlist:
                    fo.write(line)

    def sub_table(self, out, mytable, info):
        infoDict = defaultdict(int)
        with open(info, 'r') as f:
            for li in f:
                infoDict[ li.strip().split('\t', 1)[0] ] =1
        with open(mytable, 'r') as ft, open(out, 'w') as fo:
            fo.write(next(ft))
            for li in ft:
                prob = li.strip().split('\t', 1)[0].strip('"')
                if prob in infoDict:
                    fo.write(li)

if __name__ == '__main__':
    fire.Fire(Venn)