#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/6/30 14:35
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, os,re
import pandas as pd
from collections import defaultdict

def get_geneName(feature_name, gdict):
    ff = feature_name.split('_')
    dataType = ff[0]
    gene_name = 'NA'
    if dataType == "Fusion":
        gene_name =  ff[1]
    elif dataType == "APA":
        gene_id = ff[1].split('.')[0]
        gene_name = gdict[gene_id]
    elif dataType == "lncRNA" or dataType == "mRNA":
        gene_id = ff[1].split('.')[0]
        gene_name = gdict[gene_id]
    elif dataType == "AS":
        prefix, gene_id = feature_name.split(":")
        gene_name = gdict[ gene_id.split('.')[0] ]
    elif dataType == "Retro":
        gene_name = feature_name.split('_',1)[1]
    elif dataType == "RNAEditing":
        editing_site = '_'.join(ff[1:])
    elif dataType == "SNP":
        gene_name = ff[1]
    return gene_name


class myResult_summary:

    def feature(self, rf_res_fi, ginfo_fi ):

        #rf_res_fi = "/home/nbt2/proj/2020-TCGA/proc/KIRP/00_tables/pipeOne/FeatureSelection/feature(all)_importance.csv"
        #ginfo_fi = "/home/nbt2/proj/2020-TCGA/proc/KIRP/s1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv"
        rf_fi_base = os.path.basename(rf_res_fi)
        rf_fi_dir  = os.path.dirname(rf_res_fi)
        rf_fi_base = re.sub('.csv$', '-addName.csv',rf_fi_base )
        out_fi = os.path.join(rf_fi_dir, rf_fi_base)
        rf_res = pd.read_csv(rf_res_fi)
        ginfo = pd.read_csv(ginfo_fi, sep="\t")
        gg = ginfo[['gene_id', 'gene_name']]
        
        gdict = defaultdict(str)
        for i in range(len(gg)):
            gene_id = gg.loc[i, 'gene_id']
            gene_id = gene_id.split('.')[0]
            gene_name = gg.loc[i, 'gene_name']
            gdict[gene_id] = gene_name
        
        rf_res.columns = ['feature', 'weight']

        imp = rf_res[rf_res['weight'] > 0]

        imp['data_type'] = imp['feature'].str.split('_', 1).str.get(0)

        imp['gene_name'] = imp['feature'].apply(get_geneName, args=(gdict,))
        imp.sort_values(by = ['weight'], ascending=False).to_csv(out_fi, index=False)
    
    def filter_chrM(self, rawdir, geneInfo, outdir, rRNA_Fi = "/dat1/dat/ref/hg38/hg38+gencode.v32/rRNA/rRNA.info"):
        ginfo = pd.read_csv(geneInfo, sep="\t")
        ginfo_rRNA = pd.read_csv(rRNA_Fi, sep="\t")
        ginfo_chrM = ginfo[ginfo.chrom == "chrM"]
        ginfo_chrM = ginfo_chrM[['gene_id', 'gene_name']]
        ginfo_rRNA = ginfo_rRNA[['gene_id', 'gene_name']]
        ginfo_chrM = ginfo_chrM.append(ginfo_rRNA, ignore_index=True)

        gg = ginfo[['gene_id', 'gene_name']]
        gdict = defaultdict(str)
        for i in range(len(gg)):
            gene_id = gg.loc[i, 'gene_id']
            gene_id = gene_id.split('.')[0]
            gene_name = gg.loc[i, 'gene_name']
            gdict[gene_id] = gene_name

        os.system('mkdir -p %s' % outdir) 
        for entry in os.scandir(rawdir):
            print(entry.name)
            dat = pd.read_csv(entry.path, header=0)
            dat['gene_name2'] = dat.iloc[:, 0].apply( get_geneName, args=(gdict,) ) 
            fi_name_out = f'{outdir}/{entry.name}'
            dat_nochrM = dat[dat.gene_name2.isin( ginfo_chrM.gene_name ) == False ]
            dat_nochrM = dat_nochrM.drop( columns = ['gene_name2'] )
            dat_nochrM.to_csv(fi_name_out, index=False )


if __name__ == '__main__':
    fire.Fire(myResult_summary)