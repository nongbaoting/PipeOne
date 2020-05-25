#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/3/25 9:43
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire, os, re
from collections import defaultdict

def gtf_attr(attr):
    attributes = defaultdict(str)
    ats = attr.strip(' ').strip(';').split('; ')
    for at in ats:
        k,v = at.split(' ')[0:2]
        attributes[k] = v.strip('"')
    return attributes

def gtf_attrs_cat(attrsDict):
    it = []
    for key,val in attrsDict.items():
        it.append( '{} "{}"'.format(key, val) )
    attr = '; '.join(it) + ';'
    return attr

class MYRUN:
    def get_txID_geneID(self,gtf,out):
        ids = []
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if cell[2] != 'exon': continue
                attr = gtf_attr(cell[8])
                tx_gene = attr['transcript_id'] + '\t' + attr['gene_id']
                ids.append(tx_gene)
        ids = list(set(ids))
        fo = open(out, 'w')
        fo.write('txID\tgeneID\n')
        fo.write('\n'.join(sorted(ids)) + '\n')
        fo.close()

    def cpat(self, res_cpat, species, out_non="CPAT.noncoding", out_coding = "CPAT.coding"):
        cutoffset = {'human' : 0.364,
                     'mouse' : 0.44,
                     'zebrafish' : 0.381
                    }

        cutoff = float( cutoffset[species])
        print(f"species: {species}, cutoff: {cutoff}")
        fo_coding = open(out_coding, 'w')
        with open(res_cpat, 'r') as f, open(out_non,'w') as fo:
            next(f)
            for li in f:
                cell = li.strip().split('\t')
                tx_id =cell[0]
                coding_prob = float(cell[-1] )
                if coding_prob < cutoff:
                    fo.write(tx_id + "\n")
                else:
                    fo_coding.write(tx_id + '\n')
        fo_coding.close()

    def combine_coding_prediction(self, tx_gene,  CPAT_res, PLEK_res, CPPred_res,
                                  species='human', noncoding_ls ="novel_lncRNA.list"):
        cutoffset = {'human': 0.364,
                     'mouse': 0.44,
                     'zebrafish': 0.381
                     }
        CPAT_cutoff = cutoffset[species]

        d_cpat  = defaultdict(str)
        d_plek  = defaultdict(str)
        d_cppred = defaultdict(str)
        d_txGene = defaultdict(str)

        with open(tx_gene) as ftg:
            for li in ftg:
                tx_id,gene_id = li.split('\t')[0:2]
                d_txGene[tx_id ] = gene_id

        txs = []
        with open(CPAT_res, 'r') as f:
            next(f)
            for li in f:
                cell = li.strip().split('\t')
                tx_id =cell[0]
                coding_prob = float(cell[-1] )
                txs.append(tx_id)
                if coding_prob < CPAT_cutoff:
                    d_cpat[tx_id] = "noncoding"
                else:
                    d_cpat[tx_id] = "coding"

        with open(PLEK_res, 'r') as fp:
            for li in fp:
                cell = li.strip().split('\t')
                tx_id = cell[2].split(' ')[0].replace('>', '')
                if cell[0] == "Non-coding":
                    d_plek[tx_id] = "noncoding"
                else:
                    d_plek[tx_id] = "coding"

        with open(CPPred_res) as fc:
            next(fc)
            for li in fc:
                cell = li.strip().split('\t')
                tx_id = cell[0]
                d_cppred[tx_id] = cell[-2]

        fo_coding = open("TUCP.list", 'w')
        coding_locus = []
        with open("coding_potential_sum.tsv", 'w') as fo, open(noncoding_ls, 'w') as fl:
            fo.write('tx_id\tCPAT\tPLEK\tCPPred\n')
            for tx_id in txs:
                fo.write('\t'.join([tx_id, d_cpat[tx_id],d_plek[tx_id], d_cppred[tx_id]  ]) + '\n')
                if not ( d_cpat[tx_id] == d_cppred[tx_id] == d_plek[tx_id] == "noncoding"):
                    coding_locus.append(d_txGene[tx_id])

            for tx_id in txs:
                gene_id = d_txGene[tx_id ]
                if gene_id in coding_locus:
                    fo_coding.write(tx_id + '\n')
                else:
                    fl.write(tx_id + '\n')

        fo_coding.close()


    def gffcompare_tmap(self, tmap, out):
        mytx = defaultdict(list)
        mygene_not = defaultdict(int)
        mygene_len  = defaultdict(int)
        mygene_exon = defaultdict(int)

        with open(tmap) as f, open(out, 'w') as fo:
            next(f)
            for li in f:
                cell = li.strip().split('\t')
                class_code, qry_gene_id, qry_id, num_exons = cell[2:6]
                tx_len = cell[9]
                if class_code not in [ "u" ,"i" , "x"] :
                    mygene_not[qry_gene_id] = 1

                tx_len = int(tx_len)
                num_exons = int(num_exons)

                mygene_len[qry_gene_id] = max(mygene_len[qry_gene_id], tx_len)
                mygene_exon[qry_gene_id] = max(mygene_exon[qry_gene_id], num_exons )

                mytx[qry_id] = [qry_gene_id, class_code, int(num_exons), int(tx_len) ]

            for tx_id in mytx:
                gene_id, class_code, num_exons, tx_len = mytx[tx_id]
                if gene_id in mygene_not:
                    continue
                if tx_len <= 200: continue

                max_exon = mygene_exon[gene_id]
                max_len  = mygene_len[gene_id]
                ok = 0
                if max_exon >1:
                    ok = 1

                else:
                    if max_len > 2000 and class_code == "u":
                        ok =1

                if ok :
                    fo.write('\t'.join( [tx_id, gene_id, class_code,
                                        str(num_exons), str(tx_len) ]) + '\n')






if __name__ == '__main__':
    fire.Fire(MYRUN)
