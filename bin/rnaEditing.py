#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/12/17 20:38
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, os
import  re
from collections import defaultdict

class RUN:

    def cat(self,output, A2I_editing_out,  mydir = "results"):
        reg = re.compile('.SPRINT_identified_all.res')
        fo = open(output, 'w')
        fo_a2i = open(A2I_editing_out, 'w')

        fo.write('\t'.join(["chrom", "start", "end", "type", "supporting_reads", "strand", "AD:DP", "AD", "DP", "sample"]) + '\n')
        fo_a2i.write('\t'.join(["chrom", "start", "end", "type", "supporting_reads", "strand", "AD:DP", "AD", "DP", "sample"]) + '\n')

        for entry in os.scandir(mydir):

            if os.path.isdir(entry.path):continue
            sampleName = reg.sub('', entry.name)
            with open(entry.path, 'r') as f:
                next(f)
                for li in f:
                    li = li.strip()
                    chrom, start, end, Type, supporting_reads, strand,  AD_DP = li.split('\t')
                    AD,DP = AD_DP.split(':')
                    fo.write( li + '\t' + '\t'.join([AD, DP, sampleName] ) + '\n' )
                    if Type == "AG" or Type == "TC":
                        fo_a2i.write(li + '\t' + '\t'.join([AD, DP, sampleName] ) + '\n' )

    def toTable(self, res, out, cutoff=5):
        samples = []
        mydict = defaultdict(dict)

        with open(res, 'r') as f:
            next(f)
            for li in f:
                cell = li.strip().split('\t')
                chrom, start, end, Type, supporting_reads, strand, AD_DP, editing_rate, sampleName = cell
                if int(supporting_reads) >= cutoff:
                    editing_site = ':'.join([chrom, start, end, Type, strand])
                    mydict[editing_site][sampleName] = editing_rate
                if sampleName not in samples:
                    samples.append(sampleName)

        samples = sorted(set(samples))
        fo_tab = open(out, 'w')
        fo_tab.write("RES\t" + '\t'.join(samples) + '\n')
        for editing_site in mydict:
            vals = [editing_site]
            for sp in samples:
                value = 'NA'
                if sp in mydict[editing_site]:
                    value = mydict[editing_site][sp]
                vals.append(value)
            fo_tab.write('\t'.join(vals) + '\n')

    def filter_tab(self, inbed, tab, out ):
        mydict = defaultdict(int)

        with open(inbed, 'r') as f:
            for li in f:
                chrom ,start, end = li.strip().split('\t')
                mydict[f'{chrom}:{start}:{end}'] = 1

        with open(tab, 'r') as fr, open(out, 'w') as fo:
            fo.write(next(fr))
            for li in fr:
                cell = li.strip().split('\t')
                chrom,start,end = cell[0].split(':')[0:3]
                if f'{chrom}:{start}:{end}' in mydict:
                    continue
                fo.write(li )

    def filter2(self, inbed, tab, out ):
        mydict = defaultdict(int)

        with open(tab, 'r') as f:
            for li in f:
                chrom ,start, end = li.strip().split('\t')[1:3]
                mydict[f'{chrom}:{start}:{end}'] = li

        with open(inbed, 'r') as fr, open(out, 'w') as fo:
            fo.write(next(fr))
            for li in fo:
                cell = li.strip().split('\t')
                chrom,start,end,*_ = cell[0].split(':')
                if f'{chrom}:{start}:{end}' in mydict:
                    continue
                fo.write(li )



    def filter_res(self, inbed, tab, out ):
        mydict = defaultdict(int)

        with open(inbed, 'r') as f:
            for li in f:
                chrom ,start, end = li.strip().split('\t')
                mydict[f'{chrom}:{start}:{end}'] = 1

        with open(tab, 'r') as fr, open(out, 'w') as fo:
            fo.write(next(fr))
            for li in fr:
                cell = li.strip().split('\t')
                chrom,start,end= cell[0:3]
                if f'{chrom}:{start}:{end}' in mydict:
                    continue
                fo.write(li)



if __name__ == '__main__':
    fire.Fire(RUN)