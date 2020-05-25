#!/usr/bin/env python3

import fire
from collections import defaultdict

class Venn:

    def venn(self, file1,file2,out,col1=0,col2=0,header1=0,header2=0):
        lines1 = open(file1,'r').read().strip().split('\n')
        lines2 = open(file2, 'r').read().strip().split('\n')

        fo = open(out, 'w')
        start1, start2 = 0,0
        if header1:
            fo.write(lines1[0] + '\n')
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
        fo.close()


if __name__ == '__main__':
    fire.Fire(Venn)