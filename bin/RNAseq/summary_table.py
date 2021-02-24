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
    fo_lnc  = open("lncR_gene.tpm.tsv", 'w')
    fo_prot = open("prot_gene.tpm.tsv", 'w')
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
    if df.shape[0] > 1:
        df.iloc[:,0] = [f'{marks}_{i}' for i in list(df.iloc[:,0])]
        df.to_csv(outifle, header=True, index=False,sep=",")
    


def rename_col(infile, outfile, sampleDict):
    df_raw = pd.read_csv(infile)
    df_raw.rename(columns=sampleDict, inplace=True)
    df_raw.to_csv(outfile)

class RUN:

    def sep_lnc_mRNA(self, tpm, lnc_info):
        lncRNA_mRNA(tpm, lnc_info)

    def mark_feature(self, infile, outifle, marks, sep="\t"):
        add_marks(infile, outifle, marks, sep="\t")
    
    def mark_feature_spladder(self, mydir):
        for entry in os.scandir("mydir"):
            if re.search(".*confirmed.psi.txt.gz", entry.name ):
                new_name = re.sub("txt.gz$", 'csv', entry.name)
                add_marks(entry.path , f"{new_name}", "AS")

if __name__ == '__main__':
    fire.Fire(RUN)
