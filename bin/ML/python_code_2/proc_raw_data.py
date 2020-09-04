
import numpy as np
import pandas as pd
# from scipy.io import savemat
import os, fire, re
from sklearn.model_selection import train_test_split
from collections import Counter

def chck_dir(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

def fill_na(df):
    # df: sample x feature
    df_fill = df.fillna(0)
    return df_fill

def select_top_k_var(df, topk):
    # df: sample x feature
    var = df.var(axis=0, skipna=True)
    var = var.fillna(0)
    sort_id = np.argsort(-var.values)
    df_out = df.iloc[:, sort_id[0: topk]]
    return df_out

def select_top_k_nonzero_ratio(df, topk):
    # df: sample x feature
    # no nan in df
    sample_num = df.shape[0]
    non_zero_count = df.astype(bool).sum(axis=0)
    non_zero_count_ratio = non_zero_count / sample_num
    sort_id = np.argsort(-non_zero_count_ratio.values)
    df_out = df.iloc[:, sort_id[: topk]]
    return df_out

# preprocessing A
# step 1: fill nan with zeros
# step 2: select top k var
# step 3: rearrange samples

# preprocessing B
# step 1: fill nan with zeros
# step 2: drop too many zeros
# step 3: rearrange samples

def get_sample_head_idx(df):
    sample_head = re.compile("Run|Sample")
    for i, c in enumerate(df.columns):
        if sample_head.match(c):
            return i


def sub_feature_samples(indir, train_dir, sample_info, outdir):
    sinfo = pd.read_csv(sample_info)
    chck_dir(outdir)
    for entry in os.scandir(indir):
        print(entry.name)
        if re.search("clustering", entry.name): continue
        dat = pd.read_csv(entry.path, header=0, index_col=0)
        datT = dat.T
        fi_name_out = f'proc_{entry.name}'
        train_fi = os.path.join(train_dir, fi_name_out)
        dat_train = pd.read_csv(train_fi, header=0, index_col=0)

        sub_dat = datT.reindex(index = sinfo.iloc[:, 0], columns=dat_train.columns, fill_value=0)
        outfi = os.path.join(outdir, fi_name_out)
        sub_dat.to_csv(outfi, index=True, header=True)
    os.system(f"cp {sample_info} {outdir}/dat_sample_clustering.csv")


def samllSample_train_test_split(sample_df, test_size=0.3, random_state=2):
    train, test = pd.DataFrame(), pd.DataFrame()
    for group_ in set(sample_df["Group"]):
        the_df = sample_df[ sample_df["Group"] == group_ ]
        train_, test_ = train_test_split(the_df, test_size = test_size, random_state= random_state)
        train = train.append(train_, ignore_index=True )
        test = test.append(test_, ignore_index=True)
    return train, test


class MYRUN_procRaw:

    def proc(self, rawdir ,  sample_info, var_topk = 1000, tdir = "./data/proc/" ):
        #tdir = "./data/proc/"
        chck_dir(tdir)
        for entry in os.scandir(rawdir):
            if not re.search("csv$", entry.path): continue
            fileA = entry.name
            fiA = entry.path
            # fiA = "./data/rawdata/" + fileA + ".csv" #### 可以改成os.scandir
            dfA = pd.read_csv(fiA, header=0, index_col=0)
            dfA = dfA.T
            # rearrange
            clinical_fi = sample_info
            clinical_df = pd.read_csv(clinical_fi, header=0, index_col=None)

            sample_idx = get_sample_head_idx(clinical_df)
            #sample_id = clinical_df['Run']
            sample_id = clinical_df.iloc[:,sample_idx] ### 之后修改,用固定位置
            #print(sample_id)

            dfA = dfA.reindex(sample_id, fill_value=0)
            dfA = fill_na(dfA)
            # var_topk = 1000 ##### 可以改
            dfA_nonzero = select_top_k_var(dfA, var_topk)

            foutA = os.path.join( tdir, f'proc_{fileA}' )  #### 输出
            dfA_nonzero.to_csv(foutA, header=True, index=True)

    def train_test_split(self, rawdir, sample_info, test_size=0.3, random_state=2):
        sinfo = pd.read_csv(sample_info)
        train_dir, test_dir = "./data/train_dir/", "./data/test_dir/"
        chck_dir(train_dir); chck_dir(test_dir)
        sinfo_train, sinfo_test = samllSample_train_test_split(sinfo,test_size, random_state)
        sinfo_train.to_csv("./data/train_dir/train_sample_clustering.csv", index=False)
        sinfo_test.to_csv("./data/test_dir/test_sample_clustering.csv", index=False)

        for entry in os.scandir(rawdir):
            print(entry.name)
            dat = pd.read_csv(entry.path, header=0, index_col=0)
            fi_name_out = f'proc_{entry.name}'
            train_fi = os.path.join(train_dir, fi_name_out)
            test_fi  = os.path.join(test_dir, fi_name_out)

            sub_train = dat.reindex(index=sinfo_train.loc[:,'Sample'], fill_value=0)
            sub_test  = dat.reindex(index=sinfo_test.loc[:, 'Sample'], fill_value=0)

            sub_train.to_csv(train_fi)
            sub_test.to_csv(test_fi)


    def proc_sep_train_test(self, rawdir, sample_fi, var_topk =1000 ,test_size = 0.3 ):
        sample_df = pd.read_csv(sample_fi).rename(str.capitalize, axis=1)

        if "Group" not in sample_df.columns :
            print( f"Your provide file {sample_fi} must contain header: Group")
            exit(1)

        train_df, test_df = samllSample_train_test_split(sample_df, test_size)

        print(Counter(train_df.loc[:, "Group"]))
        print(Counter(test_df.loc[:, "Group"]))
        train_df.to_csv("train_sample_clustering.csv", index=False)
        test_df.to_csv("test_sample_clustering.csv", index=False)

        myrun = MYRUN_procRaw()

        train_dir  = "./data/proc/"
        test_dir   = "./data/test_proc/"
        #os.system(f'rm -rf {train_dir} {test_dir}')
        myrun.proc(rawdir ,  "train_sample_clustering.csv",  var_topk, train_dir )
        sub_feature_samples(rawdir, train_dir, "test_sample_clustering.csv", test_dir)
        os.system(f"cp train_sample_clustering.csv {train_dir}")


if __name__ == '__main__':
    fire.Fire(MYRUN_procRaw)