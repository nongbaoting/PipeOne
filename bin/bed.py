#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/10/3 17:06
# @Author  : Baoting Nong
# @Email   : 523135753@qq.com
import fire,os,re
from collections import defaultdict

class BED6:
    def __init__(self,line):
        cell = line.split('\t')
        self.chrom = cell[0]
        self.start = int(cell[1])
        self.end   = int(cell[2])
        self.feature= cell[3]
        self.score = int(cell[4])
        self.strand = cell[5]

        self.id = f'{self.chrom}:{self.start}:{self.end}:{self.strand}'
        self.length = self.end - self.start
        self.circRNA_id = f'{self.chrom}:{self.start}|{self.end}'

def id_count(bedfile,bedDict):
    lines = open(bedfile, 'r').read().strip().split('\n')
    for li in lines:
        bed = BED6(li)
        #print(bed.id)
        bedDict[bed.id] +=1

class Bed:

    def select_common(self, freq,outfile, *beds):
        bedDict = defaultdict(int)
        for bedfile in beds:
            id_count(bedfile, bedDict)
        fo =  open(outfile, 'w')
        for bed_id, val in bedDict.items():
            if val >= freq:
                chrom,start,end,strand = bed_id.split(':')
                fo.write('\t'.join([chrom,start, end, bed_id, '0', strand]) + '\n')
    def venn(self,bed, circRNA,outfile):
        ids = set([li for li in open(circRNA, 'r').read().strip().split('\n') ])
        with open(bed, 'r') as fb, open(outfile, 'w') as fo:
            for bl in fb:
                bed = BED6(bl)
                if bed.circRNA_id in ids:
                    fo.write(bl)
    def annotate_byHommer(self, bed,out,genome='hg19'):
        hommer = '~/bin/homer/bin/annotatePeaks.pl'
        import tempfile
        tf_bed = open('bed_tmp', 'w')
        lines = open(bed, 'r').read().strip().split('\n')
        for li in  lines:
            bb = BED6(li)
            li_1 = '\t'.join([bb.chrom, str(bb.start),str(bb.start), bb.feature + ':start', str(bb.score), bb.strand])
            li_2 = '\t'.join([bb.chrom, str(bb.end),str(bb.end), bb.feature + ':end', str(bb.score), bb.strand])
            ne_li = li_1 + '\n' + li_2 + '\n'
            tf_bed.write(ne_li)
        tf_bed.close()
        os.system(f'{hommer} bed_tmp hg19 > result_tmp')
        annDict = defaultdict(list)
        with open('result_tmp', 'r') as fr, open(out, 'w') as fo:
            reg_rep = re.compile(':(start|end)$')
            next(fr)
            for li in fr:
                cell = li.strip('\n').split('\t')
                ids = cell[0]
                id = re.sub(reg_rep, '', ids)
                annDict[id].append(cell)
            header = ['chrom', 'start','end','id','score','strand','circType','distance_TSS','n_promoter','n_entrez',
                      'n_unigene', 'n_refseq','n_ensembl', 'n_geneName', 'n_geneAlias', 'n_genedis', 'n_geneType']
            fo.write('\t'.join(header) + '\n')
            for id in annDict:
                ce = annDict[id]
                ce = sorted(ce,key=lambda a: int(a[2]) )
                ids, chrom, start, end, strand, pscore, ratio, anno, detail_anno, distance_TSS, \
                n_promoter, n_entrez, n_unigene, n_refseq, n_ensembl, n_geneName, n_geneAlias, n_genedis, n_geneType = ce[0]

                ids_2, chrom_2, start_2, end_2 = ce[1][0:4]
                anno_2, detail_anno_2,distance_TSS_2 =ce[1][7:10]

                distance = distance_TSS
                if int(distance_TSS_2) < int(distance_TSS_2):
                    distance = distance_TSS_2

                circType = anno.split(' ')[0] + '-' + anno_2.split(' ')[0]

                circ_line = '\t'.join([ chrom, end, end_2, id, '0', strand, circType,distance_TSS,n_promoter, n_entrez,
                                        n_unigene, n_refseq, n_ensembl,n_geneName, n_geneAlias, n_genedis, n_geneType ])
                fo.write(circ_line + '\n')


        os.system('rm *tmp')

if __name__ == '__main__':
    fire.Fire(Bed)

