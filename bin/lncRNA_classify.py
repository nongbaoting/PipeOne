#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/7/29 9:16
__author__ = 'Baoting Nong'
__email__ = '523135753@qq.com'
import fire,os
from collections import defaultdict

def get_lnc_cls(closest):
    lnccls = defaultdict(list)
    with open(closest, 'r') as f:
        for li in f:
            cell = li.strip().split('\t')
            lnc = cell[3]
            cls_protein = cell[9]
            distance  = int(cell[-1])
            strand_1 = cell[5]
            strand_2 = cell[11]
            if lnc in lnccls:
                if abs(distance) < abs(lnccls[lnc][1]):
                    lnccls[lnc] = [cls_protein, distance, strand_1, strand_2 ]
            else:
                lnccls[lnc] = [cls_protein, distance, strand_1, strand_2]
    return lnccls

def getProteinInfo(protein_info):
    geneName = defaultdict(str)
    with open(protein_info, 'r') as f:
        for li in f:
            cell = li.strip().split('\t')
            geneName[cell[1] ] =cell[2]
    return geneName

class RUN:

    def final(self, lncClass, closest, closest_sense, lnc_info, protein_info, mydist = 2000):
        clsDict = get_lnc_cls(closest)
        cls_senseDict = get_lnc_cls(closest_sense)
        geneNameDict = getProteinInfo(protein_info)

        lines = open(lncClass, 'r').read().strip().split('\n')
        lncDict = defaultdict(str)
        for li in lines:
            cell = li.split('\t')
            lncDict[cell[0] ] = cell[1]

        with open(lnc_info, 'r') as f, open('lncRNA_class_closest_PCG.tsv', 'w') as fo:
            head =next(f).strip().split('\t')
            head_add = [ "lncRNA_type_4", "lncRNA_type", "closest_PCG","closest_PCG_Name", "strand_PCG", "distance_PCG",
                                "closest_PCG_sense", "closest_PCG_Name_sense", "strand_PCG_sense", "distance_PCG_sense"]
            head.extend(head_add)
            fo.write('\t'.join(head) + '\n')
            for li in f:
                cell = li.strip().split('\t')
                lnc  = cell[1]
                lnc_class =lncDict[lnc]
                cls_pcg, distance, strand_1, strand_2  = clsDict[lnc]
                cls_pcg_name = geneNameDict[cls_pcg]
                cls_pcg_sence, distance_sence, strand_1_sense, strand_2_sense = cls_senseDict[lnc]
                cls_pcg_sence_name = geneNameDict[cls_pcg_sence]
                class_2k = lnc_class
                if lnc_class == "Intergenic" and abs(distance) < mydist:
                    if strand_1 == strand_2:
                        class_2k = "cis-regulatory"
                    else:
                        class_2k = 'bidirectional'
                cell_add = [lnc_class, class_2k , cls_pcg,  cls_pcg_name, strand_2, str(distance),
                            cls_pcg_sence, cls_pcg_sence_name, strand_2_sense, str(distance_sence)]
                cell.extend(cell_add)

                fo.write('\t'.join(cell) + '\n')

if __name__ == '__main__':
    fire.Fire(RUN)