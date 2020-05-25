#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/5/19 11:13
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, os, csv, gzip
from collections import  defaultdict

def lncRNA_mRNA(tpm, lnc_info):
    mylnc = []
    with open(lnc_info, 'r')  as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter="\t")
        for item in csv_reader:
            mylnc.append(item['gene_id'])

    fo_lnc = gzip.open("s1_lncR_mRNA/lncR_gene.tpm.csv.gz", 'wt')
    fo_prot = gzip.open("s1_lncR_mRNA/prot_gene.tpm.csv.gz", 'wt')
    spamwriter = csv.writer(fo_lnc)
    spamwriter_prot = csv.writer(fo_prot)
    with  open(tpm, 'r') as csv_file_2:
        csv_tpm = csv.reader(csv_file_2, delimiter = "\t" )
        header = next(csv_tpm)
        spamwriter.writerow(header)
        spamwriter_prot.writerow(header)
        for row in csv_tpm:
            if row[0] in mylnc:
                spamwriter.writerow(row)
            else:
                spamwriter_prot.writerow(row)

    fo_lnc.close()
    fo_prot.close()

class RUN:

    def merge(self,):
        home = os.path.dirname(__file__)
        ## 1,2 lncRNA mRNA
        # lncRNA_mRNA("../s1_lncRNA/results/salmon/salmon_gene_tpm.tsv", "../s1_lncRNA/results/novel_lncRNA/all_lncRNA_info.tsv")

        ## s3 APA
        # os.system(f"Rscript {home}/apa_3utr_filter.R")
        ## s4 retro
        os.system("mkdir -p s4_retrotranscriptome")
        os.system(f"cp ../s4_retrotranscriptome/results/telescope/telescope.FPKM-divide_totalMapReads.csv s4_retrotranscriptome/FPKM-divide_totalMapReads.csv")

        ## s5 Fusion
        # os.system("mkdir -p s5_fusion")
        # os.system(f"cp ../s5_fusion/results/arriba/fusion_arriba_out.tsv s5_fusion/fusion_arriba_out.tsv")

        ## s6_rnaEditing
        os.system("mkdir -p s6_rnaEditing")
        os.system(f"cp ../s6_rnaEditing/results/sprint/Merge/SPrint_A2I_table.tsv s6_rnaEditing/SPrint_A2I_table.tsv")

        ## s7
        os.system("mkdir -p s7_alternative_splicing")
        os.system(f"cp ../s7_alternative_splicing/ s7_alternative_splicing/")

        ## s8_SNP
        os.system("mkdir -p s8_SNP")
        os.system(f"cp ../s8_SNP/ s8_SNP/")


if __name__ == '__main__':
    fire.Fire(RUN)
