#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/2/2 17:46
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, os, gzip
import re
from collections import  defaultdict

def get_geneIDName(info):
    md = defaultdict(str)
    with open(info, 'r') as f:
        heads = next(f).strip().split('\t')
        index_gene_id = heads.index('gene_id')
        index_gene_name = heads.index('gene_name')

        for li in f:
            cell = li.strip().split('\t')
            gene_id = cell[index_gene_id].split('.')[0]
            gene_name = cell[index_gene_name]
            md[gene_id] = gene_name
    return md

def scanAndFind_pattern(mydir, pattern):
    wantdir = []
    for entry in os.scandir(mydir):
        if pattern.search(entry.name):
            wantdir.append(entry)
        elif entry.is_dir() :
            wantdir.extend( scanAndFind_pattern(entry.path, pattern) )
    return wantdir

class RUN:

    def spladder_PSI_tab(self, infile):
        infile_basename = os.path.basename(infile)
        outfile  = re.sub("txt.gz","psi.txt.gz", infile_basename)

        f  = gzip.open(infile, 'rt')
        fo = gzip.open(outfile, 'wt')
        heads = next(f).strip().split('\t')
        event_id_index = heads.index('event_id')
        gene_name_index = heads.index('gene_name')
        psi_indices = [i for i, s in enumerate(heads) if ':psi' in s]
        new_heads = [heads[i].split(':')[0] for i in  psi_indices]
        new_heads.insert(0, 'event_id:gene_name')
        fo.write('\t'.join( new_heads ) + '\n')
        for li in f:
            cell = li.strip().split('\t')
            event_id = cell[event_id_index]
            gene_name = cell[gene_name_index]
            new_id = f'{event_id}:{gene_name}'
            new_cell = [cell[i] for i in psi_indices ]
            new_cell.insert(0, new_id)
            fo.write('\t'.join(new_cell) + '\n' )
        f.close()
        fo.close()

    def feature_geneName(self, info, outfile, pattern = "merge.*txt.gz"):
        gdict = get_geneIDName(info)
        fo = open(outfile, 'w')
        fo.write('\t'.join(['feature_id', 'gene_name']) + '\n')
        for entry in scanAndFind_pattern('.', re.compile(pattern)):
            f = gzip.open(entry.path, 'rt')
            for li in f:
                feature_id = li.split('\t')[0]
                gene_id = feature_id.split(':')[1].split('.')[0]
                if gene_id in gdict:
                    fo.write('\t'.join([feature_id, gdict[gene_id]] ) + '\n')

            f.close()
        fo.close()

if __name__ == '__main__':
    fire.Fire(RUN)
