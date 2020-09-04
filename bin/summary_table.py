#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/5/19 11:13
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, os, csv, re
from collections import  defaultdict
import pandas as pd
import numpy as np

def chck_dir(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

def lncRNA_mRNA(tpm, lnc_info):
    mylnc = []
    with open(lnc_info, 'r')  as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter="\t")
        for item in csv_reader:
            mylnc.append(item['gene_id'])

    chck_dir("s1.1_lncR_mRNA")
    fo_lnc  = open("s1.1_lncR_mRNA/lncR_gene.tpm.tsv", 'w')
    fo_prot = open("s1.1_lncR_mRNA/prot_gene.tpm.tsv", 'w')
    with  open(tpm, 'r') as csv_file_2:
        csv_tpm = csv.reader(csv_file_2, delimiter = "\t" )
        header = next(csv_tpm)
        header.insert(0,'gene_id')
        fo_lnc.write('\t'.join( header) + '\n')
        fo_prot.write('\t'.join(header) + '\n')
        for row in csv_tpm:
            if row[0] in mylnc:
                fo_lnc.write('\t'.join(row) + '\n')
            else:
                fo_prot.write('\t'.join(row) + '\n')

    fo_lnc.close()
    fo_prot.close()

def add_marks(infile, outifle, marks, sep="\t" ):
    df = pd.read_csv(infile, sep=sep)
    if df.shape[0] >0:
        df.iloc[:,0] = [f'{marks}_{i}' for i in list(df.iloc[:,0])]
        df.to_csv(outifle, header=True, index=False,sep=",")

def rename_col(infile, outfile, sampleDict):
    df_raw = pd.read_csv(infile)
    df_raw.rename(columns=sampleDict, inplace=True)
    df_raw.to_csv(outfile)

class RUN:

    def check_tables(self, library = 'polyA'):
        home = os.path.dirname(__file__)
        ## 1,2 lncRNA mRNA
        lncRNA_mRNA("../s1.1_lncRNA/results/salmon/salmon_gene_tpm.tsv", "../s1.1_lncRNA/results/novel_lncRNA/all_lncRNA_info.tsv")
        
        ## s2 circRNA
        if library == "total":
            chck_dir("s1.2_circRNA")
            os.system("cp ../s1.2_circRNA/results/CIRIquant/circRNA_cpm.csv s1.2_circRNA/circRNA_cpm.csv")

        ## s3 APA
        os.system(f"Rscript --vanilla {home}/apa_3utr_filter.R")
        ## s4 retro
        chck_dir("s1.4_retrotranscriptome")
        os.system(f"cp ../s1.4_retrotranscriptome/results/telescope/telescope.FPKM-divide_totalMapReads.csv s1.4_retrotranscriptome/FPKM-divide_totalMapReads.csv")

        ## s5 Fusion
        chck_dir("s1.5_fusion")
        os.system(f"cp ../s1.5_fusion/results/arriba/fusion_arriba_out.tsv s1.5_fusion/fusion_arriba_out.tsv")

        ## s6_rnaEditing
        chck_dir("s1.6_rnaEditing")
        os.system(f"cp ../s1.6_rnaEditing/results/sprint/Merge/SPrint_A2I_table.annovar.csv s1.6_rnaEditing/SPrint_A2I_table.annovar.csv")

        ## s7 AS
        chck_dir("s1.7_alternative_splicing")
        os.system(f"cp ../s1.7_alternative_splicing/results/spladder_out_table/*confirmed.psi.txt.gz s1.7_alternative_splicing/")

        ## s8_SNP
        chck_dir("s1.8_SNP")
        os.system(f"cp ../s1.8_SNP/results/annovar_table/snp.geneBase.tsv s1.8_SNP/snp.geneBase.tsv")

    def mark_feature(self, library = 'polyA'):
        chck_dir('00_rawdata')
        add_marks("s1.1_lncR_mRNA/prot_gene.tpm.tsv", "00_rawdata/prot_gene.tpm.csv", "mRNA")
        add_marks("s1.1_lncR_mRNA/lncR_gene.tpm.tsv", "00_rawdata/lncR_gene.tpm.csv", "lncRNA")
        if library == "total":
            add_marks("s1.2_circRNA/circRNA_cpm.csv", "00_rawdata/circRNA_cpm.csv", "circRNA")
        add_marks("s1.3_APA-3TUR/pau_results.filterPau-distal-proximal.txt", "00_rawdata/APA_pau-distal-proximal.csv", "APA")
        add_marks("s1.4_retrotranscriptome/FPKM-divide_totalMapReads.csv", "00_rawdata/retro-FPKM-divide_totalMapReads.csv", "Retro", sep=",")
        add_marks("s1.5_fusion/fusion_arriba_out.tsv", "00_rawdata/fusion_arriba_out.csv", "Fusion")
        add_marks("s1.6_rnaEditing/SPrint_A2I_table.annovar.csv", "00_rawdata/RNA-editing-rate.csv", "RNAEditing",sep=",")

        for entry in os.scandir("s1.7_alternative_splicing"):
            if re.search(".*confirmed.psi.txt.gz", entry.name ):
                new_name = re.sub("txt.gz$", 'csv', entry.name)
                add_marks(entry.path , f"00_rawdata/{new_name}", "AS")

        add_marks("s1.8_SNP/snp.geneBase.tsv", "00_rawdata/snp.geneBase.csv", "SNP")

    def rename_sample(self,sample_sheet):
        outdir = "./00_rawdata_tcga"
        chck_dir(outdir)

        df_s = pd.read_csv(sample_sheet, sep="\t")

        for entry in os.scandir("00_rawdata"):
            print(entry.name)
            rename_col(entry.path, f'{outdir}/{entry.name}', sampleDict )

    def tcga_sample_sheet(self):
        pass


if __name__ == '__main__':
    fire.Fire(RUN)
