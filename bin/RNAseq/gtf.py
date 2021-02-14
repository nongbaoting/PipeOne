#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/9/29 18:58
# @Author  : Baoting Nong
__author__ = 'Baoting Nong'
__email__ = '523135753@qq.com'
#import pysam
from collections import defaultdict
import fire, os, re, intervals
import operator

class GTF:
    def __init__(self,gtffile,attr_type):
        gtfDict = defaultdict(list)
        with open(gtffile, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                line = line.strip()
                cell = line.split('\t')
                if cell[2] != 'exon': continue
                attrs = gtf_attr(cell[8])
                gene_id = attrs[attr_type]
                gene_id = gene_id.split('.')[0]
                gtfDict[gene_id].append(line)
        self.gtfDict = gtfDict

    def get_start_end(self,attr_id):
        start = -1
        end = -1
        for li in self.gtfDict[attr_id]:
            exon = EXON(li)
            if  start == -1:
                start = exon.start
            elif start > exon.start:
                start = exon.start
            if end == -1:
                end =exon.end
            elif end < exon.end:
                 end = exon.end
        return start, end

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

def to_gff_attr(attrsDict):
    '''
    only exon
    '''
    attr = f'Parent={ attrsDict["transcript_id"] };gene_id={ attrsDict["gene_id"] };transcript_id={ attrsDict["transcript_id"]};'
    return attr

def to_gtf_line(chrom,source,feature,start,end,score,strand,frame,attrDict):
    attr_cat = gtf_attrs_cat(attrDict)
    line = '\t'.join([chrom,source,feature,str(start),str(end),score,strand,frame,attr_cat])
    return line

class EXON:
    def __init__(self,cell):
        self.chrom   = cell[0]
        self.source  = cell[1]
        self.feature = cell[2]
        self.start   = int(float( cell[3]))
        self.end     = int(float( cell[4]))
        self.score   = cell[5]
        self.strand  = cell[6]
        self.frame   = cell[7]
        self.attr    = cell[8]
        self.gene_type = 'NA'
        attrDict     = gtf_attr(self.attr)
        self.gene_id       = attrDict['gene_id']
        self.gene_name = self.gene_id
        if 'gene_name' in attrDict:
            self.gene_name = attrDict['gene_name']
        if 'transcript_id' in attrDict:
            self.transcript_id = attrDict['transcript_id']
        if 'gene_type' in attrDict:
            self.gene_type = attrDict['gene_type']
        if 'transcript_type' in attrDict:
            self.transcript_type = attrDict['transcript_type' ]

from pprint import pprint
class GTF_circ:
    def __init__(self, tabixfile, chrom, start, end):
        self.tabixfile = tabixfile
        self.chrom     = chrom
        self.start     = start
        self.end       = end
        self.rows = tabixfile.fetch(chrom, start, end, parser = pysam.asTuple() )
        self.gene_id = ''
        self.gene_name = ''
        self.transcript_id = ''
        for row in self.rows:
            exon = EXON(row)
            self.gene_id = exon.gene_id
            self.gene_name = exon.gene_name
            break

    def check_exons(self):
        up_exons = self.tabixfile.fetch(self.chrom, self.start,self.start + 1,parser=pysam.asTuple() )
        dw_exons = self.tabixfile.fetch(self.chrom, self.end -1 ,  self.end ,  parser=pysam.asTuple() )

        pos1,pos2=0,0
        if up_exons :
            for up in up_exons:
                if up[2] == 'exon':
                    pos1 = 1
        if dw_exons:
            for dw in dw_exons:
                if dw[2] == 'exon':
                    pos2 = 1
        circType = ''
        if pos1 and pos2:
            circType = 'exon'
        elif pos1 or pos2:
            circType ='exon-intron'
        else:
            'intron'
        return circType

class GTF_exe:

    def add_gene_name(self, gtf, out_gtf):
        fo = open(out_gtf, 'w')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                attr = gtf_attr(cell[8])
                if 'gene_name' not in attr:
                    attr['gene_name'] = attr['gene_id']
                attr_cat = gtf_attrs_cat(attr)
                cell[8] = attr_cat
                new_line = '\t'.join(cell)
                fo.write(new_line + '\n')
        fo.close()

    def simple(self, gtf, out):
        fo = open(out, 'w')
        mydict = defaultdict(list)
        mytype = ['exon', 'CDS', 'start_codon', 'stop_codon']
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if int(cell[3]) >= int(cell[4]): continue
                if cell[2] not in mytype:continue
                attr = gtf_attr(cell[8])
                #new_dict =defaultdict(int)
                #new_dict['gene_id'] = attr['gene_id']
                tx_id = attr['transcript_id']
                #attr_cat = gtf_attrs_cat(new_dict)
                #cell[8] = attr_cat
                mydict[ tx_id].append(cell)

        for tx_id in sorted( mydict.keys() ):
            allExon = mydict[tx_id]
            for cell in sorted(allExon,key=lambda k:int(k[3])):
                new_line = '\t'.join(cell)
                fo.write(new_line + '\n')
        fo.close()

    def add_quotes(self,gtf,out):
        fo = open(out, 'w')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line):continue
                if re.match('track', line):continue
                cell = line.strip().split('\t')
                if int(cell[3]) >= int(cell[4]):continue
                attr =  gtf_attr(cell[8])
                attr_cat = gtf_attrs_cat(attr)
                cell[8]  = attr_cat
                new_line = '\t'.join(cell)
                fo.write(new_line + '\n')
        fo.close()

    def gene_id2name(self, gtf, out):
        fo = open(out, 'w')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if int(cell[3]) >= int(cell[4]): continue
                attr = gtf_attr(cell[8])
                attr['gene_id'] = attr['gene_name']
                attr_cat = gtf_attrs_cat(attr)
                cell[8] = attr_cat

                new_line = '\t'.join(cell)
                fo.write(new_line + '\n')
        fo.close()

    def get_by_trans_id(self, gtf, info, out, info_column=0):
        reg_sp = re.compile('\s')
        fo = open(out, 'w')
        lines = open(info, 'r').read().strip().split('\n')
        ids = [reg_sp.split(line)[info_column]  for line in lines]
        ids = set(ids)

        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line):  continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                exon = EXON(cell)
                if exon.transcript_id in ids:
                    fo.write(line)
        fo.close()

    def get_by_gene_id(self, gtf, info, out, info_column=0):
        reg_sp = re.compile('\s')
        fo = open(out, 'w')
        lines = open(info, 'r').read().strip().split('\n')
        ids = [reg_sp.split(line)[info_column]  for line in lines]
        ids = set(ids)

        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line):  continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                exon = EXON(cell)
                if exon.gene_id in ids:
                    fo.write(line)
        fo.close()

    def get_by_type(self,gtf, attr_type,attr_value, out):
        fo = open(out, 'w')
        vals = attr_value.split(',')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line):  continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                attr = gtf_attr(cell[8])
                attr_cat = gtf_attrs_cat(attr)
                if attr[attr_type] in vals:
                    fo.write(line)
        fo.close()

    def get_id(self,gtf, attr_type, out):
        fo = open(out, 'w')
        wanted = []
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                attr = gtf_attr(cell[8])
                val =attr[attr_type]
                wanted.append(val)
        wanted = set(wanted)
        fo.write('\n'.join(sorted(wanted)) + '\n')
        fo.close()

    def start_must_less_than_end(self,gtf,out):
        fo = open(out, 'w')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if int(cell[3]) >= int(cell[4]): continue
                if cell[2] != 'exon': continue
                attr = gtf_attr(cell[8])
                new_attr = defaultdict(str)
                new_attr['transcript_id'] = attr['transcript_id']
                new_attr['gene_id'] = 'G' + attr['gene_id']
                attr_cat = gtf_attrs_cat(new_attr)
                cell[8] = attr_cat
                new_line = '\t'.join(cell)
                fo.write(new_line + '\n')
        fo.close()

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

    def get_txID_geneID_geneName(self,gtf,out):
        ids = []
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if cell[2] != 'exon': continue
                attr = gtf_attr(cell[8])
                tx_gene = attr['transcript_id'] + '\t' + attr['gene_id'] + '\t' + attr['gene_name']
                ids.append(tx_gene)
        ids = list(set(ids))
        fo = open(out, 'w')
        fo.write('txID\tgeneID\tgeneName\n')
        fo.write('\n'.join(sorted(ids)) + '\n')
        fo.close()

    def get_exonLen(self, gtf, out):
        fo = open(out, 'w')
        fo.write('gene_id\ttx_id\texon_id\texonLen\n')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if cell[2] != 'exon': continue
                if int(cell[3]) >= int(cell[4]): continue
                exonLen = int( cell[4] ) - int(cell[3]) + 1
                attr = gtf_attr(cell[8])
                tx_id = attr['transcript_id']
                gene_id = attr['gene_id']
                exon_id = f'{tx_id}.{cell[3]}.{cell[4]}'
                new_line  = '\t'.join([gene_id, tx_id, exon_id, str(exonLen)])
                fo.write(new_line + '\n')
        fo.close()

    def to_info(self, gtf, out_tsv):
        f = open(gtf, 'r')
        genedict = defaultdict(list)
        for line in f:
            cell = line.strip().split('\t')
            if re.match('#', line): continue
            if not re.match('exon', cell[2]): continue
            # cell[8]='gene_id "ENSG00000237375"; transcript_id "ENST00000327822"; exon_number "11"; gene_name "BX072566.1"; gene_biotype "protein_coding";'
            attributes = gtf_attr(cell[8])
            start = int(float(cell[3])) - 1
            end = int(float(cell[4]))
            strand = cell[6]
            exon_len = end - start
            gene_type = 'NA'
            if 'gene_type' in attributes:
                gene_type = attributes['gene_type']
            elif 'gene_biotype' in attributes:
                gene_type = attributes['gene_biotype']
            gene_name = attributes['gene_id']
            if 'gene_name' in attributes:
                gene_name = attributes['gene_name']

            genedict[attributes['transcript_id']].append( [attributes['gene_id'], gene_name, cell[0], start, end,strand,
                                               gene_type,exon_len, cell[1] ])
        fo = open(out_tsv, 'w')
        fo.write('tx_id\tgene_id\tgene_name\tchrom\tstart\tend\tstrand\texons\tsplice_len\tgene_type\tsource\n')
        for transcript_id in genedict:
            exons = genedict[transcript_id]
            exons_sorted = sorted(exons,key= lambda k : k[3])

            splice_len = 0
            tx_start, tx_end = exons_sorted[0][3:5]
            for ex in exons_sorted:
                gene_id, gene_name, chrom, start, end, strand, gene_type, exon_len, gene_source = ex
                splice_len += exon_len
                if tx_start > start: tx_start = start
                if tx_end <  end :   tx_end   = end

            gene_id, gene_name, chrom, start, end,strand, gene_type, exon_len, gene_source = exons[0]
            exon_num = len(exons)
            genomic_len = str(tx_end - tx_start + 1)
            #print(genomic_len)
            fo.write('\t'.join([transcript_id, gene_id, gene_name, chrom, str(tx_start),
                               str(tx_end), strand,str(exon_num), str(splice_len), gene_type, gene_source ]) + '\n')
        fo.close()

    def to_TSS_bed(self,gtf,out):
        f = open(gtf, 'r')
        genedict = defaultdict(list)
        for line in f:
            cell = line.strip().split('\t')
            if re.match('#', line): continue
            if not re.match('exon', cell[2]): continue
            exon = EXON(cell)
            genedict[exon.transcript_id ].append( exon )

        fo_tss_200 = open(f'{out}.200.bed', 'w')
        fo_tss_first = open(f'{out}.firstExon.bed', 'w')
        fo_tss_1500 = open(f'{out}.1500.bed', 'w')
        fo_tss_2k = open(f'{out}.2k.bed', 'w')
        fo_tss_3k = open(f'{out}.3k.bed', 'w')

        for transcript_id in genedict:
            exons = genedict[transcript_id]
            exons.sort(key = lambda k : k.start)
            first_exon = exons[0]
            tss = first_exon.start
            if first_exon.strand =='-':
                first_exon = exons[-1]
                tss = first_exon.end

            fo_tss_first.write('\t'.join([first_exon.chrom, str(first_exon.start -1), str(first_exon.end), transcript_id, '0',
                                          first_exon.strand]) + '\n')
            fo_tss_200.write('\t'.join([first_exon.chrom, str( (tss - 200) if tss >200 else 0), str(tss + 200), transcript_id, '0',
                                        first_exon.strand]) + '\n')
            fo_tss_1500.write('\t'.join([first_exon.chrom, str(  (tss - 1500) if tss >1500 else 0 ), str(tss + 1500), transcript_id, '0',
                                         first_exon.strand]) + '\n')
            fo_tss_2k.write('\t'.join(
                [first_exon.chrom, str( (tss - 2000) if tss > 2000 else 0), str(tss + 2000), transcript_id, '0',
                 first_exon.strand]) + '\n')

            fo_tss_3k.write('\t'.join(
                [first_exon.chrom, str( (tss - 3000) if tss > 3000 else 0), str(tss + 3000), transcript_id, '0',
                 first_exon.strand]) + '\n')

        fo_tss_200.close()
        fo_tss_first.close()
        fo_tss_1500.close()
        fo_tss_2k.close()
        fo_tss_3k.close()

    def to_extend_bed(self, gtf, out, length = 100000):
        f = open(gtf, 'r')
        genedict = defaultdict(list)
        for line in f:
            cell = line.strip().split('\t')
            if re.match('#', line): continue
            if not re.match('exon', cell[2]): continue
            exon = EXON(cell)
            genedict[exon.transcript_id].append(exon)

        fo = open(out, 'w')
        for transcript_id in genedict:
            exons = genedict[transcript_id]
            exons.sort(key=lambda k: k.start)
            first_exon = exons[0]
            last_exon  = exons[-1]

            fo.write(
                '\t'.join([first_exon.chrom, str((first_exon.start  - length ) if first_exon.start > length else 0),
                           str(last_exon.end + length), transcript_id, '0',
                           first_exon.strand]) + '\n')
        fo.close()

    def get_tx_frame(self, gencode_gtf,out):
        with open(gencode_gtf, 'r') as f, open(out, 'w') as fo:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if not cell[2] =='transcript':continue
                attr = gtf_attr(cell[8])
                txID = attr['transcript_id']
                frame = cell[7]
                print(txID, frame)

    def to_gff(self,gtf,out):
        '''
        only exon
        '''
        fo = open(out, 'w')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if cell[2] != 'exon': continue
                attr = gtf_attr(cell[8])
                attr_cat = to_gff_attr(attr)
                cell[8] = attr_cat
                new_line = '\t'.join(cell)
                fo.write(new_line + '\n')
        fo.close()

    def annotate_circRNA(self, gtf, all_annote_ciri2_ce2, out):
        if not os.path.exists('gtf.sort.gz'):
            os.system(f'grep -v ^"#" {gtf} | sort -k1,1 -k4,4n | bgzip > gtf.sort.gz')
            os.system(f'tabix -p gff gtf.sort.gz')
        tabixfile = pysam.TabixFile('gtf.sort.gz')
        fo = open(out, 'w')
        lines = open(all_annote_ciri2_ce2, 'r').read().strip().split('\n')
        for li in lines[1:]:
            circ_id, geneName, txid, circType = li.split('\t')[0:4]
            chrom,start,end,strand = circ_id.split(':')
            start = int(start)
            end   = int(end)
            if circType != 'intergenic':
                circO = GTF_circ(tabixfile,chrom,start,end)
                if circO.gene_id:
                    circType = circO.check_exons()
            fo.write('\t'.join([circ_id, geneName, txid, circType]) + '\n')

    def anno_circID(self, gtf, circRNA, out):
        sorted_gtf = os.path.basename(gtf) + '.sort.gz'
        if not os.path.exists( sorted_gtf ):
            os.system(f'grep -v ^"#" {gtf} | sort -k1,1 -k4,4n | bgzip >{sorted_gtf}')
            os.system(f'tabix -p gff {sorted_gtf}')
        tabixfile = pysam.TabixFile(sorted_gtf)
        fo = open(out, 'w')
        f =  open(circRNA, 'r')
        reg_chrom = re.compile(r'_|\.')
        fo.write('circ_id_len\tid_relate\tgene_name\tcircType\n')
        for li in f:
            li = li.strip('\n')
            circ_id = li.split('\t')[0]
            chrom,pos = circ_id.split(':')
            if re.search(reg_chrom, chrom):continue
            start,end = pos.split('|')[0:2]
            start = int(start)
            end   = int(end)
            circO = GTF_circ(tabixfile, chrom, start, end)
            circType = 'intergenic'
            gene_name = 'NA'
            if circO.gene_id:
                circType = circO.check_exons()
                gene_name = circO.gene_name
            fo.write('\t'.join([li, gene_name, circType]) + '\n')
        fo.close()
        f.close()

    def circ_possible_tx(self, gtf, circRNA, out):
        sorted_gtf = os.path.basename(gtf) + '.sort.gz'
        if not os.path.exists(sorted_gtf):
            os.system(f'grep -v ^"#" {gtf} | sort -k1,1 -k4,4n | bgzip >{sorted_gtf}')
            os.system(f'tabix -p gff {sorted_gtf}')
        tabixfile = pysam.TabixFile(sorted_gtf)
        fo = open(out, 'w')
        f = open(circRNA, 'r')
        reg_chrom = re.compile(r'_|\.')
        fo.write('circ_id_len\tid_relate\tgene_name\tcircType\n')
        for li in f:
            li = li.strip('\n')
            circ_id = li.split('\t')[0]
            chrom, pos = circ_id.split(':')
            if re.search(reg_chrom, chrom): continue
            start, end = pos.split('|')[0:2]
            start = int(start)
            end = int(end)
            circO = GTF_circ(tabixfile, chrom, start, end)
            circType = 'intergenic'
            gene_name = 'NA'
            if circO.gene_id:
                circType = circO.check_exons()
                gene_name = circO.gene_name
            fo.write('\t'.join([li, gene_name, circType]) + '\n')
        fo.close()
        f.close()

    def circ_possible_tx2(self, gtf, circRNA, out):
        cir = open(circRNA, 'r').read().strip().split('\n')
        mygtf = GTF(gtf, 'gene_id')
        fo = open(out, 'w')
        for circ_line in cir[1:]:
            cell = circ_line.split('\t')
            genename = cell[1]
            #chrom, start, end, strand = cell[0].split(':')
            chrom, pos = cell[0].split(':')
            start, end = pos.split('|')
            start = int(start); end = int(end)
            circInterval  = intervals.closedopen(start, end)
            if not genename in mygtf.gtfDict: continue
            exonLines = mygtf.gtfDict[genename]
            mytempDict = defaultdict(list)
            for el in exonLines:
                myexon = EXON(el.split('\t'))
                mytempDict[myexon.transcript_id].append(myexon)

            mycircDict = defaultdict(list)
            for mytx in mytempDict:
                circExons = []
                exons = mytempDict[mytx]
                exons_sort = sorted(exons, key=lambda k: k.start)
                start_loc,end_loc = 0,0
                for index, this_exon in enumerate(exons_sort):
                    # interval from bed
                    exonInterval = intervals.closedopen( this_exon.start - 1, this_exon.end )
                    if start in exonInterval:
                        start_loc = 1
                    if end -1 in exonInterval:
                        end_loc = 1
                    newExonInterval = exonInterval & circInterval
                    if newExonInterval:
                        circExons.append(newExonInterval)
                    if index + 1 < len(exons_sort):
                        next_exon = exons_sort[index + 1]
                        intronInterval = intervals.closedopen(this_exon.end + 1, next_exon.start - 1)
                        newIntronInterval = intronInterval & circInterval
                        if newIntronInterval:
                            circExons.append(intronInterval)
                mytype = start_loc + end_loc
                tx_circ_id = f'{exons_sort[0].transcript_id}^{cell[0]}'
                attrDict = {'gene_id': exons_sort[0].gene_id, 'transcript_id':  tx_circ_id,}

                new_circ_line_list =[]
                for new_exon in circExons:
                    new_exon = new_exon.to_atomic()
                    new_circ_line = to_gtf_line(exons_sort[0].chrom, exons_sort[0].source, exons_sort[0].feature,
                                                new_exon.lower + 1, new_exon.upper , exons_sort[0].score, exons_sort[0].strand,
                                                exons_sort[0].frame, attrDict)
                    #print(new_circ_line, mytype)
                    new_circ_line_list.append(new_circ_line)
                mycircDict[tx_circ_id] = [mytype,new_circ_line_list]
            mycon = 0
            for mytxcirc,value in sorted(mycircDict.items(), key=lambda x: x[1][0], reverse=True):
                mytype,new_circ_line_list = value
                if mytype == 2:
                    fo.write('\n'.join(new_circ_line_list) + '\n')
                    print(mytype)
                    mycon = 1
                elif mytype ==1 and mycon == 0:
                    fo.write('\n'.join(new_circ_line_list) + '\n')
                    print(mytype)
                    mycon = 1
                elif mytype ==0  and mycon == 0:
                    fo.write('\n'.join(new_circ_line_list) + '\n')
                    print(mytype)
                    mycon = 1

        fo.close()

    def protein_gene_name(self, gtf, out):
        fo = open(out, 'w')
        pro = []
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line):  continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                attr = gtf_attr(cell[8])
                if attr['gene_type'] ==  "protein_coding":
                    if attr['protein_id'] in pro:continue
                    pro.append(attr['protein_id'] )
                    fo.write('\t'.join([attr['protein_id'], attr['gene_name'], attr['gene_id']]) + '\n')
        fo.close()

    def to_bed(self, gtf, out, only_exon=1, name = 'gene_id'):
        '''
        :param gtf:
        :param out:
        :return:
        '''
        fo = open(out, 'w')
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if cell[2] != 'exon' and only_exon: continue
                exon = EXON(cell)
                feature_id = exon.transcript_id
                if name == 'gene_id': 
                    feature_id = exon.gene_id
                elif name == 'transcript_id':  
                    feature_id = exon.transcript_id
                elif name == 'gene_name':
                    feature_id = exon.gene_name

                bed = [exon.chrom, str(exon.start - 1), str(exon.end), feature_id, exon.score, exon.strand]
                fo.write('\t'.join(bed) + '\n')
        fo.close()


    def to_bed_PLAIDOH(self, out, gtf, expr):
        exprDict = defaultdict(list)
        header = ['#CHR', 'START', 'STOP', 'NAME', 'TYPE']
        fo = open(out, 'w')
        with open(expr, 'r') as f:
            head = next(f).strip('\n').strip().split('\t')
            headss = [f'sample{ i+ 1}' for i in range(len(head) )]
            header.extend(headss)
            fo.write('\t'.join(header) + '\n')
            #exit()
            for li in f:
                cell = li.strip().split('\t')
                exprDict[cell[0] ] = cell[1:]

        itemss = []
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                if cell[2] != 'exon':continue
                exon = EXON(cell)
                gene_type = 'lncRNA'
                if exon.gene_type == 'protein_coding':
                    gene_type = 'protein_coding'
                item = f'{exon.chrom}-{exon.start}-{exon.end}'
                if exon.gene_id in exprDict and item not in itemss:
                    itemss.append(item )
                    bed = [exon.chrom, str(exon.start), str(exon.end), exon.gene_name,gene_type]
                    bed.extend(exprDict[exon.gene_id] )
                    fo.write('\t'.join(bed) + '\n')
        fo.close()

    def to_bedGene_PLAIDOH(self, out, gtf, expr):
        exprDict = defaultdict(list)
        header = ['#CHR', 'START', 'STOP', 'NAME', 'TYPE']
        with open(expr, 'r') as f:
            head = next(f).strip('\n').strip().split('\t')
            headss = [f'sample{ i+ 1}' for i in range(len(head) )]
            header.extend(headss)
            #exit()
            for li in f:
                cell = li.strip().split('\t')
                exprDict[cell[0] ] = cell[1:]

        geneDict = defaultdict(list)
        fo_all = open(f"{out}.all.tmp" , 'w')
        fo_all.write('\t'.join(header) + '\n')

        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                exon = EXON(cell)
                gene_type = 'lncRNA'
                if exon.gene_type == 'protein_coding':
                    gene_type = 'protein_coding'
                if exon.gene_id in exprDict:
                    bed = [exon.chrom, str(exon.start), str(exon.end), exon.gene_name,gene_type]
                    bed.extend(exprDict[ exon.gene_id ])
                    fo_all.write('\t'.join(bed) + '\n')


                if exon.gene_id in geneDict:
                    start = min(geneDict[exon.gene_id][1], exon.start)
                    end = max(geneDict[exon.gene_id][2], exon.end)
                    geneDict[exon.gene_id] = [exon.chrom, start, end, exon.gene_name,gene_type]
                else:
                    geneDict[exon.gene_id] = [exon.chrom, exon.start, exon.end, exon.gene_name,gene_type]

        fo = open(f'{out}.tmp', 'w')
        fo.write('\t'.join(header) + '\n')
        for gene_id in geneDict:
            if gene_id in exprDict:
                bed = [ str(i) for i in geneDict[gene_id] ]
                bed.extend(exprDict[gene_id])
                fo.write('\t'.join(bed) + '\n')
        fo.close()
        fo_all.close()

        os.system(f" sort -k1,1 -k2,2n {out}.tmp >{out}")
        os.system(f" sort -k1,1 -k2,2n {out}.all.tmp >{out}.all")

    def to_gene_bed(self, gtf, out):
        geneDict = defaultdict(list)
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                exon = EXON(cell)
                if exon.gene_id in geneDict:
                    start = min(geneDict[exon.gene_id][1], exon.start - 1)
                    end = max(geneDict[exon.gene_id][2], exon.end)
                    geneDict[exon.gene_id] = [exon.chrom, start, end, exon.gene_id, '.', exon.strand]
                else:
                    geneDict[exon.gene_id] = [exon.chrom, exon.start - 1, exon.end, exon.gene_id, '.', exon.strand]

        with open(out, 'w') as fo:
            for cell in geneDict.values():
                fo.write('\t'.join([str(i) for i in cell]) + '\n' )

    def to_gene_bedName(self, gtf, out):
        geneDict = defaultdict(list)
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue
                cell = line.strip().split('\t')
                exon = EXON(cell)
                if exon.gene_name in geneDict:
                    start = min(geneDict[exon.gene_name][1], exon.start - 1)
                    end = max(geneDict[exon.gene_name][2], exon.end)
                    geneDict[exon.gene_name] = [exon.chrom, start, end, exon.gene_name, '.', exon.strand]
                else:
                    geneDict[exon.gene_name] = [exon.chrom, exon.start - 1, exon.end, exon.gene_name, '.', exon.strand]

        with open('bed_tmp.txt', 'w') as fo:
            for cell in geneDict.values():
                fo.write('\t'.join([str(i) for i in cell]) + '\n' )

        os.system(f'sort -k1,1 -k2,2n bed_tmp.txt > {out}; rm bed_tmp.txt' )

    def to_gencode2biomark5(self,gtf, out):
        f = open(gtf, 'r')
        fo = open(out, 'w')
        fo.write("Gene stable ID\tTranscript stable ID\tGene type\tTranscript type\tGene name\n")
        txList = []
        for line in f:
            cell = line.strip().split('\t')
            if re.match('#', line): continue
            if not re.match('exon', cell[2]): continue
            exon = EXON(cell)
            tx_id = exon.transcript_id
            if tx_id not in txList:
                gene_id_stable = exon.gene_id.split('.')[0]
                tx_id_stable = exon.transcript_id.split('.')[0]
                fo.write('\t'.join([gene_id_stable, tx_id_stable, exon.gene_type, exon.transcript_type , exon.gene_name] ) + '\n')

    def to_txDict(self, gtf, out):
        geneDict = defaultdict(list)
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                if re.match('track', line): continue

                cell = line.strip().split('\t')
                if cell[2] != 'exon': continue
                exon = EXON(cell)
                if exon.transcript_id in geneDict:
                    start = min(geneDict[exon.transcript_id][3], exon.start - 1)
                    end = max(geneDict[exon.transcript_id][4], exon.end)
                    geneDict[exon.transcript_id] = [exon.chrom, exon.source, exon.feature, start, end, exon.strand,exon.score,  exon.frame, exon.gene_name]
                else:
                    geneDict[exon.transcript_id] = [exon.chrom, exon.source, exon.feature, exon.start, exon.end, exon.score, exon.strand,  exon.frame,  exon.gene_name]

        gtfDict  = defaultdict(list)
        for tx_id in geneDict:
            chrom, source, feature, start, end, score, strand,  frame,  gene_name =geneDict[tx_id] 
            gtfDict["Chromosome"].append(chrom)
            gtfDict["Source"].append(source)
            gtfDict["Feature"].append(feature)
            gtfDict["Start"].append(start)
            gtfDict["End"].append(end)
            gtfDict["Score"].append(score)
            gtfDict["Strand"].append(strand)
            gtfDict["Frame"].append(frame)
            gtfDict["gene_name"].append(gene_name)
            gtfDict["transcript_id"].append(tx_id)
        
if __name__ == '__main__':
    fire.Fire(GTF_exe)
