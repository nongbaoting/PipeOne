#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/6/29 17:14
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, re, os
import pandas as pd

def chck_dir(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

class MYRUN_subset_sample:

    def subset(self,indir, outdir,sample_info):
        sinfo = pd.read_csv(sample_info)
        chck_dir(outdir)
        for entry in os.scandir(indir):
            print(entry.name)
            dat = pd.read_csv(entry.path , index_col=0)
            outfi = os.path.join(outdir, entry.name)
            #sub_dat = dat[dat.iloc[:,0].isin(sinfo.iloc[:,0]) ]
            sub_dat = dat.reindex(index=sinfo.iloc[:,0] ,fill_value=0)
            sub_dat.to_csv(outfi) 
        os.system(f"cp {sample_info} {outdir}/dat_sample_clustering.csv")

    def sub_feature_samples(self, indir, train_dir,  outdir, sample_info):
        sinfo = pd.read_csv(sample_info)
        chck_dir(outdir)
        for entry in os.scandir(indir):
            print(entry.name)
            dat = pd.read_csv(entry.path, header=0, index_col=0 )
            datT = dat.T
            fi_name_out = f'proc_{entry.name}'
            train_fi = os.path.join(train_dir, fi_name_out )
            dat_train = pd.read_csv(train_fi, header=0, index_col=0 )

            sub_dat = datT.reindex(index= sinfo.iloc[:,0], columns=dat_train.columns, fill_value=0)
            outfi = os.path.join(outdir, fi_name_out)
            sub_dat.to_csv(outfi, index=True, header=True)
        os.system(f"cp {sample_info} {outdir}/dat_sample_clustering.csv")

if __name__ == '__main__':
    fire.Fire(MYRUN_subset_sample)

