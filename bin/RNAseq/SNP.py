#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/2/20 0:07
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'

import fire, re, os
from  collections import defaultdict

def scanAndFind_pattern(mydir, mypattern):
    wantdir = []
    for entry in os.scandir(mydir):
        if mypattern.search(entry.name):
            wantdir.append(entry)
        elif entry.is_dir() :
            wantdir.extend( scanAndFind_pattern(entry.path, mypattern) )

    return wantdir


class MYRUN:

    def gene_base_table(self, outfile, annovar_dir):
        reg_anno = re.compile("multianno.txt$")
        myfiles = scanAndFind_pattern(annovar_dir, reg_anno)
        non_silence = ['frameshift deletion', 'frameshift insertion', 'nonframeshift deletion',
                       'nonframeshift insertion', 'nonsynonymous SNV', 'stopgain', 'stoploss']
        mydict = defaultdict(dict)
        allSam = []
        for entry in myfiles:
            print(entry.path)
            sample = entry.name.split(".")[0]
            allSam.append(sample)
            with open(entry.path, 'r') as f:
                next(f)
                for li in f:
                    cell = li.strip().split('\t')
                    func = cell[8]
                    gene_name = cell[6]
                    if func in non_silence:
                        mydict[gene_name][sample] = 1

        allSam = sorted(allSam)

        with open(outfile, 'w') as fo:
            fo.write('gene_name\t' + '\t'.join(allSam) + '\n')
            for gene_name in mydict:
                cell = [gene_name]
                for sam  in allSam:
                    if sam in mydict[gene_name]:
                        cell.append('1')
                    else:
                        cell.append('0')
                fo.write('\t'.join(cell) + '\n')


    def annovar_res(self, res, outfile):
        fo = open(outfile, 'w')
        f  = open(res, 'r')
        head = next(f)
        fo_geneName = open('feature_geneName.tsv', 'w')
        fo.write('\t'.join(["a2i_id", 'gene_name', 'func_gene', 'exonFunc' ] ) + '\n' )
        for li in f:
            cell = li.strip().split('\t')
            exonFunc = cell[8]
            chrom, start, end,ref, alt, func_gene, gene_name = cell[0:7]
            a2i_id = ':'.join([chrom, start, end ])
            fo.write('\t'.join([a2i_id, gene_name, func_gene, exonFunc] ) + '\n' )
            fo_geneName.write('\t'.join([a2i_id, gene_name]) + '\n')

        fo.close()
        f.close()
        fo_geneName.close()



if __name__ == '__main__':
    fire.Fire(MYRUN)
