#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/2/20 23:29
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, re, os
from  collections import defaultdict

def bed2dict(bed, mydict_int ):
    mydict = defaultdict(str)
    with open(bed, 'r') as f:
        for li in f:
            chrom, start, end, circ, score, strand = li.split('\t')[0:6]
            circ_id = f'{chrom}:{start}|{end}|{strand}'
            if not circ_id in mydict:
                mydict[ circ_id ] = li
                mydict_int[circ_id] += 1

    return mydict, mydict_int

class MYRUN:

    def ciri2(self, circ_bed, circ_res):
        fo_bed = open(circ_bed, 'w')
        lines = open(circ_res, 'r').read().strip().split('\n')
        for line in lines[1:]:
            cell = line.split('\t')
            circrna_id_1base = cell[0]
            chrom, pos = circrna_id_1base.split(':')
            start, end = pos.split('|')
            start = int(start) - 1
            strand = cell[10]
            circrna_id = f'{chrom}:{start}|{end}|{strand}'
            junc_read = cell[4]
            circrna_id = re.sub('aaAaa', '__', circrna_id)
            fo_bed.write('\t'.join([chrom, str(start), end, circrna_id, junc_read, strand]) + '\n')
        fo_bed.close()

    def circexplorer2(self,  circ_bed, circ_res):
        fo_bed = open(circ_bed, 'w')
        lines = open(circ_res, 'r').read().strip().split('\n')
        for line in lines[1:]:
            cell = line.split('\t')
            chrom, start, end, jun = cell[0:4]
            strand = cell[5]
            #start = int(start) - 1
            circrna_id = f'{chrom}:{start}|{end}|{strand}'
            junc_read = jun.split('/')[-1]
            circrna_id = re.sub('aaAaa', '__', circrna_id)
            fo_bed.write('\t'.join([chrom, str(start), end, circrna_id, junc_read, strand]) + '\n')
        fo_bed.close()

    def dcc(self,  circ_bed, circ_res):
        fo_bed = open(circ_bed, 'w')
        lines = open(circ_res, 'r').read().strip().split('\n')
        for line in lines[1:]:
            chrom,start, end, Gene, JunctionType, strand,  =line.split('\t')[0:6]
            circrna_id = f'{chrom}:{start}|{end}|{strand}'
            circrna_id = re.sub('aaAaa', '__', circrna_id)
            fo_bed.write('\t'.join([chrom, str(start), end, circrna_id, '.', strand ]) + '\n')
        fo_bed.close()

    def filter(self,  outbed,  ciri2_bed, ce2_bed, dcc_bed):
        cNum = defaultdict(int)
        ciri2_dict, cNum = bed2dict(ciri2_bed, cNum)
        ce2_dict, cNum   = bed2dict(ce2_bed, cNum)
        dcc_dict, cNum   = bed2dict(dcc_bed, cNum)

        with open(outbed, 'w') as fo:
            for circ_id, num in cNum.items():
                if num >= 2:
                    if circ_id in ciri2_dict:
                        fo.write(ciri2_dict[circ_id] )
                    elif circ_id in ce2_dict:
                        fo.write(ce2_dict[circ_id] )
                    elif circ_id in dcc_dict:
                        fo.write(dcc_dict[circ_id] )
                    else:
                        print('some thing wrong ! ' + outbed )
                        exit(1)

if __name__ == '__main__':
    fire.Fire(MYRUN)
