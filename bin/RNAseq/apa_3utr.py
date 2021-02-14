#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/12/31 21:17
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, re, os
from Bio import SeqIO
from collections import defaultdict

def scanAndFind(mydir, myfile_name):
    wantdir = []
    for entry in os.scandir(mydir):
        if entry.name == myfile_name:
            wantdir.append(entry)
        elif entry.is_dir():
            wantdir.extend(scanAndFind(entry.path, myfile_name ) )

    return wantdir

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

class RUN:

    def replace_SalmonIndex_ID(self, salmon_index_id,  mydir, myfile = "quant.sf"):

        mydict = defaultdict(str)
        with open(salmon_index_id, 'r') as fs:
            for li in fs:
                pseudo_id, rel_id = li.strip().split('\t')[0:2]
                mydict[pseudo_id] = rel_id

        for entry in scanAndFind(mydir, myfile):
            new_file = os.path.join(os.path.dirname(entry.path), "quant_replace.sf")
            with open(entry.path, 'r') as f, open(new_file, 'w') as fo:
                fo.write(next(f) )
                for li in f:
                    cell = li.strip('\n').split('\t')
                    cell[0] = mydict[cell[0] ]
                    fo.write('\t'.join(cell) + '\n')

    def shorten_id(self, fa, outFa, sep = '_'):
        count = 0
        reg_sp = re.compile(sep)
        fo_tab = open('fa_id_relate.txt', 'w')
        fo = open(outFa, 'w')
        for record in SeqIO.parse(fa, 'fasta'):
            new_id = record.id
            if len(record.id) >200:
                dess = reg_sp.split(record.id)
                if len(dess) > 1 and len( dess[0] ) < 200:
                    new_id = f'{ dess[0] }-{count}'
                else:
                    new_id = f'new_id-{count}'
                count +=1
            print(new_id, record.id, sep='\t', file=fo_tab)
            record.id = new_id
            SeqIO.write(record, fo, 'fasta')

    def relate_to_geneName(self, infile, gene_info, outfile):
        gdict = get_geneIDName(gene_info)
        reg_sp = re.compile('_|\.')
        with open(infile, 'r') as f, open(outfile, 'w') as fo:
            fo.write('\t'.join(['feature_id','gene_name']) + '\n')
            for li in f:
                feature_id = li.split('\t')[0]
                gene_id = reg_sp.split(feature_id)[0]
                if gene_id in gdict:
                    gene_name = gdict[gene_id ]
                    fo.write('\t'.join([feature_id, gene_name]) + '\n')


if __name__ == '__main__':
    fire.Fire(RUN)