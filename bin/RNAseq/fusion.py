#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/12/23 10:46
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'


import fire,os,re
from collections import defaultdict

def getline(myfile,mydict):
    lines = open(myfile,'r').read().strip().split('\n')
    for li in  lines:
        cell = li.split('\t')
        #circ_id = ':'.join(cell[0:3]) + ':' + cell[5]
        circ_id = f'{cell[0]}:{cell[1]}|{cell[2]}'
        mydict[circ_id] = cell[:]
    return mydict

def mytrim( mname ):
    cell= mname.split(',')
    ag = []
    for g in cell:
        gene = g.split('(')[0]
        ag.append(gene)
    new = '^'.join(ag)
    return new

def tumorSample(datainfo):
    lines = open(datainfo, 'r').read().strip().split('\n')
    tumor = []
    for li in lines[1:]:
        cell = li.split('\t')
        if cell[3] =='Tumor':
            tumor.append(cell[1])
    return tumor

def scanAndFind(mydir, myfile):
    wantdir = []
    for entry in os.scandir(mydir):
        if entry.name == myfile:
            wantdir.append(entry)
        elif entry.is_dir() :
            wantdir.extend( scanAndFind(entry.path, myfile) )
    return wantdir

def scanAndFind_pattern(mydir, pattern):
    wantdir = []
    for entry in os.scandir(mydir):
        if pattern.search(entry.name):
            wantdir.append(entry)
        elif entry.is_dir() :
            wantdir.extend( scanAndFind_pattern(entry.path, pattern) )
    return wantdir

def write_to_bed(all_fusion, out_prefix = 'all_fusion'):
    fo_bed1 = open(f'{out_prefix}.1.bed', 'w')
    fo_bed2 = open(f'{out_prefix}.2.bed', 'w')
    header = 'chr\tstart\tend\tfusion\tSample\n'
    fo_bed1.write(header)
    fo_bed2.write(header)
    for li in all_fusion:
        cell = li.strip('\n').split('\t')
        fusion, mysample,gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, type, direction1, direction2, \
        split_reads1, split_reads2, discordant_mates, *rest = cell
        chrom1, end1 = breakpoint1.split(':')
        chrom2, end2 = breakpoint2.split(':')
        chrom1 = 'chr' + chrom1
        chrom2 = 'chr' + chrom2
        start1 = str(int(end1) - 1)
        start2 = str(int(end2) - 1)
        fo_bed1.write('\t'.join([chrom1, start1, end1, fusion, mysample]) + '\n')
        fo_bed2.write('\t'.join([chrom2, start2, end2, fusion, mysample]) + '\n')

    fo_bed1.close()
    fo_bed2.close()

def sp_fusion_gene(gene, sample , myd):
    reg_gene = re.compile('\(\d+?\)')
    new_gene = reg_gene.sub('', gene)
    sp_gene = new_gene.split(',')
    if len(sp_gene) > 1:
        for gg in sp_gene:
            if gg != 'NA':
                myd[gg][sample] = 2
    else:
        myd[ sp_gene[0] ][sample] = 1
    return myd

class RUN:

    def arriba(self, outfile, *mydirs):
        tumor = tumorSample('datainfo_4.txt')
        reg_file = re.compile('.fusions.tsv')
        allfile = []
        for myd in mydirs:
            for entry in os.scandir(myd):
                if re.search(reg_file, entry.name):
                    allfile.append(entry)

        fusionDict = defaultdict(dict)
        samples = []
        freq = defaultdict(int)

        all_new, all_head = [], []
        for entry in allfile:
            mysample = re.sub(reg_file, '', entry.name)
            mysample = re.sub('-', '.', mysample)
            lines = open(entry.path, 'r').read().strip().split('\n')
            all_head = lines[0]
            for li in lines[1:]:
                cell = li.split('\t')
                gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, type, direction1, direction2, \
                split_reads1, split_reads2, discordant_mates, confidence, closest_genomic_breakpoint1, \
                closest_genomic_breakpoint2, filters, fusion_transcript, read_identifiers = cell
                # if confidence == 'low': continue
                if confidence != 'high': continue
                if site1 == 'intergenic':
                    gene1 = mytrim(gene1)
                if site2 == 'intergenic':
                    gene2 = mytrim(gene2)
                fusion = f'{gene1}__{gene2}'
                fusionDict[fusion][mysample] = 1
                samples.append(mysample)
                new_line = [fusion, mysample]
                new_line.extend(cell)
                new_line_cat = '\t'.join(new_line)
                all_new.append(new_line_cat)

        write_to_bed(all_new)

        gene_freq = defaultdict(int)
        fusion_freq = defaultdict(int)
        samples = sorted(list(set(samples)))
        fo = open(outfile, 'w')
        fo.write('\t'.join(samples) + '\n')
        reg_hat = re.compile('\^')
        for fusion in sorted(fusionDict.keys()):
            all_val = [fusion]
            for s in samples:
                val = '0'
                if s in fusionDict[fusion]:
                    fusion_freq[fusion] += 1
                    val = '1'
                    ##counts
                    gene1, gene2 = fusion.split('__')
                    if s in tumor:
                        for g1 in reg_hat.split(gene1):
                            gene_freq[g1] += 1
                        for g2 in reg_hat.split(gene2):
                            gene_freq[g2] += 1
                all_val.append(val)

            fo.write('\t'.join(all_val) + '\n')

        fo.close()

        with open('gene_freq.txt', 'w') as foc:
            foc.write('gene\tcounts\tfusion\tfusion_counts\n')
            for fusion in fusion_freq:
                fusion_count = fusion_freq[fusion]
                gene1, gene2 = fusion.split('__')
                for g1 in reg_hat.split(gene1):
                    foc.write(f'{g1}\t{gene_freq[g1]}\t{fusion}\t{fusion_count}\n')
                for g2 in reg_hat.split(gene2):
                    foc.write(f'{g2}\t{gene_freq[g2]}\t{fusion}\t{fusion_count }\n')

        with open('arriba_all.txt', 'w') as fo_all:
            fo_all.write('fusion\tsample\t' + all_head + '\n')
            fo_all.write('\n'.join(all_new) + '\n')


    def arriba_table(self, outfile,   mydir):
        reg_file = re.compile('.fusions.tsv')
        resFiles = scanAndFind_pattern(mydir, re.compile(r'.*.fusions.tsv') )


        myd = defaultdict(dict)
        samples = []



        for entry in resFiles:
            print(entry.name)
            mysample = reg_file.sub( '', entry.name)
            samples.append(mysample)
            lines = open(entry.path, 'r').read().strip().split('\n')
            all_head = lines[0]
            for li in lines[1:]:
                cell = li.split('\t')
                gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, type, direction1, direction2, \
                split_reads1, split_reads2, discordant_mates, confidence = cell[0:15]
                if confidence == 'low': continue
                myd = sp_fusion_gene(gene1, mysample, myd)
                myd = sp_fusion_gene(gene2, mysample, myd)

        samples = sorted(list(set(samples) ) )
        with open(outfile, 'w') as fo:
            fo.write('gene_name\t' + '\t'.join(samples) + '\n')
            for myg in myd:
                myvalues = [myg ]
                for mys in samples:
                    if mys in myd[myg]:
                        myvalues.append( str( myd[myg][mys] ) )
                    else:
                        myvalues.append( "0" )
                fo.write('\t'.join(myvalues) + '\n')


    def star_fusion(self):
        pass


if __name__ == '__main__':
    fire.Fire(RUN)
