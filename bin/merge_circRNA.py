#!/usr/bin/env python3
import os, fire, re
from collections import defaultdict


def rowSums(cell, start, end, cutoff=1):
    num = 0
    for i in cell[start:end]:
        if int(i) >= cutoff:
            num += 1
    return str(num)


def rowSums_lines(lines, start, end, cutoff=1):
    num = 0
    newLines = []
    for line in lines:
        cell = line.strip('\n').split('\t')
        num = rowSums(cell, start, end, cutoff)
        newLines.append(line.strip() + '\t' + num)
    return newLines


def to_dict(lines):
    mydict = defaultdict(dict)
    header = lines[0].split('\t')
    for line in lines[1:]:
        cell = line.split('\t')
        for i in range(1, len(cell)):
            mydict[cell[0]][header[i - 1]] = cell[i]
    return mydict


def to_dict_liens(lines, soft, oldLines):
    header = lines[0].split('\t')
    new_header = '\t'.join([soft + '_' + i for i in header]) + '\t' + soft + '.rowSums'
    for line in lines[1:]:
        cell = line.split('\t')
        num = rowSums(cell, 1, len(cell))
        val = '\t'.join(cell[1:]) + '\t' + str(num)
        for i in range(1, len(cell)):
            oldLines[cell[0]][soft] = val
    oldLines['header'][soft] = new_header
    return oldLines


class ResultCircRNA:
    def ciri2(self, outfile='ciri2.table.txt', *files):
        fo_bed = open('ciri2.bed', 'w')
        read_table = defaultdict(dict);
        all_sample = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            for line in lines[1:]:
                cell = line.split('\t')

                circrna_id_1base = cell[0]
                chrom, pos = circrna_id_1base.split(':')
                start, end = pos.split('|')
                start = int(start) - 1
                strand = cell[10]
                circrna_id = f'{chrom}:{start}|{end}'
                # circrna_id = f'{chrom}:{start}|{end}:{strand}'
                junc_read = cell[4]
                circrna_id = re.sub('aaAaa', '__', circrna_id)
                read_table[circrna_id][sample] = junc_read

                fo_bed.write('\t'.join([chrom, str(start), end, 'CIRI2', '0', strand]) + '\n')
        fo_bed.close()

        out_lines = []
        out_lines.append('\t'.join(all_sample))
        for circrna_id in read_table.keys():
            new_line = circrna_id
            for sample in all_sample:
                count = '0'
                if sample in read_table[circrna_id].keys():
                    count = read_table[circrna_id][sample]
                new_line += '\t' + count
            out_lines.append(new_line)

        with open(outfile, 'w') as f:
            f.write('\n'.join(out_lines) + '\n')

        newLines_sum = rowSums_lines(out_lines[1:], 1, len(out_lines[0].strip().split('\t')) + 2)
        with open(outfile + '.rowSums', 'w') as f:
            f.write('circRNA_id\t' + out_lines[0].strip() + '\t' + 'ciri2_rowSums' + '\n')
            f.write('\n'.join(newLines_sum) + '\n')

    def circRNA_finder(self, outfile='circRNA_finder.table.txt', *files):
        read_table = defaultdict(dict);
        all_sample = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            for line in lines[:]:
                cell = line.split('\t')
                # cell[1] = str( int(cell[1]) + 1)
                circrna_id = ':'.join(cell[0:2]) + '|' + cell[2]
                junc_read = cell[4]
                circrna_id = re.sub('aaAaa', '__', circrna_id)
                read_table[circrna_id][sample] = junc_read
        # pprint(read_table)
        out_lines = []
        out_lines.append('\t'.join(all_sample))
        for circrna_id in read_table.keys():
            new_line = circrna_id
            for sample in all_sample:
                count = '0'
                if sample in read_table[circrna_id].keys():
                    count = read_table[circrna_id][sample]
                new_line += '\t' + count
            out_lines.append(new_line)

        with open(outfile, 'w') as f:
            f.write('\n'.join(out_lines) + '\n')

        newLines_sum = rowSums_lines(out_lines[1:], 1, len(out_lines[0].strip().split('\t')) + 2)
        with open(outfile + '.rowSums', 'w') as f:
            f.write('circRNA_id\t' + out_lines[0].strip() + '\t' + 'circRNA_finder_rowSums' + '\n')
            f.write('\n'.join(newLines_sum) + '\n')

    def find_circ(self, outfile='find_circ.table.txt', *files):
        read_table = defaultdict(dict)
        all_sample = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            if len(lines[0]) == 0:
                print('warining: find_circ did not found result in sample: ' + sample)
                continue
            for line in lines[:]:
                cell = line.split('\t')
                # cell[1] = str( int(cell[1]) + 1)
                circrna_id = ':'.join(cell[0:2]) + '|' + cell[2]
                junc_read = cell[4]
                circrna_id = re.sub('aaAaa', '__', circrna_id)
                read_table[circrna_id][sample] = junc_read
        # pprint(read_table)
        out_lines = []
        out_lines.append('\t'.join(all_sample))
        for circrna_id in read_table.keys():
            new_line = circrna_id
            for sample in all_sample:
                count = '0'
                if sample in read_table[circrna_id].keys():
                    count = read_table[circrna_id][sample]
                new_line += '\t' + count
            out_lines.append(new_line)

        with open(outfile, 'w') as f:
            f.write('\n'.join(out_lines) + '\n')

        newLines_sum = rowSums_lines(out_lines[1:], 1, len(out_lines[0].strip().split('\t')) + 2)
        with open(outfile + '.rowSums', 'w') as f:
            f.write('circRNA_id\t' + out_lines[0].strip() + '\t' + 'find_circ_rowSums' + '\n')
            f.write('\n'.join(newLines_sum) + '\n')

    def circexplorer2(self, outfile='circexplorer2.table.txt', *files):
        read_table = defaultdict(dict)
        all_sample = []
        fo_bed = open('circexplorer2.bed', 'w')
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            for line in lines[1:]:
                cell = line.split('\t')
                chrom, start, end, jun = cell[:4]
                strand = cell[5]
                tran = cell[15]
                circrna_id = f'{chrom}:{start}|{end}'
                #circrna_id = f'{chrom}:{start}|{end}:{strand}:{tran}'
                junc_read = jun.split('/')[-1]
                circrna_id = re.sub('aaAaa', '__', circrna_id)
                read_table[circrna_id][sample] = junc_read

                circrna_id_tran = f'{chrom}:{start}|{end}:{strand}:{tran}'
                fo_bed.write('\t'.join([chrom, str(start), end, 'circexplorer2', '0', strand]) + '\n')
        fo_bed.close()
        out_lines = []
        out_lines.append('\t'.join(all_sample))
        for circrna_id in read_table.keys():
            new_line = circrna_id
            for sample in all_sample:
                count = '0'
                if sample in read_table[circrna_id].keys():
                    count = read_table[circrna_id][sample]
                new_line += '\t' + count
            out_lines.append(new_line)

        with open(outfile, 'w') as f:
            f.write('\n'.join(out_lines) + '\n')

        newLines_sum = rowSums_lines(out_lines[1:], 1, len(out_lines[0].strip().split('\t')) + 2)
        with open(outfile + '.rowSums', 'w') as f:
            f.write('circRNA_id\t' + out_lines[0].strip() + '\t' + 'ciri2_rowSums' + '\n')
            f.write('\n'.join(newLines_sum) + '\n')

    def dcc(self, outfile='dcc.table.txt', *files):
        read_table = defaultdict(dict);
        all_sample = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            for line in lines[1:]:
                chrom,start, end, junc_read =line.split('\t')
                circrna_id = f'{chrom}:{start}|{end}'
                # circrna_id = f'{chrom}:{start}|{end}:{strand}'
                circrna_id = re.sub('aaAaa', '__', circrna_id)
                read_table[circrna_id][sample] = junc_read
        out_lines = []
        out_lines.append('\t'.join(all_sample))
        for circrna_id in read_table.keys():
            new_line = circrna_id
            for sample in all_sample:
                count = '0'
                if sample in read_table[circrna_id].keys():
                    count = read_table[circrna_id][sample]
                new_line += '\t' + count
            out_lines.append(new_line)
        with open(outfile, 'w') as f:
            f.write('\n'.join(out_lines) + '\n')
        newLines_sum = rowSums_lines(out_lines[1:], 1, len(out_lines[0].strip().split('\t')) + 2)
        with open(outfile + '.rowSums', 'w') as f:
            f.write('circRNA_id\t' + out_lines[0].strip() + '\t' + 'ciri2_rowSums' + '\n')
            f.write('\n'.join(newLines_sum) + '\n')

    def dcc_bed(self, *files):
        fo_bed = open('dcc.bed', 'w')
        ids = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            for line in lines[1:]:
                chrom, start, end, _, _, strand,*more =line.split('\t')
                circrna_id = f'{chrom}:{start}|{end}'
                # circrna_id = f'{chrom}:{start}|{end}:{strand}'
                circrna_id = re.sub('aaAaa', '__', circrna_id)
                if circrna_id in ids:continue
                fo_bed.write('\t'.join([chrom, str(start), end, 'CIRI2', '0', strand]) + '\n')
        fo_bed.close()

    def merge_tables(self, outtable, *circTable):
        newDict = defaultdict(dict)
        oldLines = defaultdict(dict)
        new_head = [];
        softs = []
        for myfile in circTable:
            lineS = open(myfile, 'r').read().strip().split('\n')
            mydict = to_dict(lineS)
            header = lineS[0].split('\t')
            soft = re.split(r'\.|_table', myfile)[0]
            softs.append(soft)
            oldLines = to_dict_liens(lineS, soft, oldLines)
            print(soft)
            for circRNA_id in mydict.keys():
                if 'soft' in newDict[circRNA_id]:
                    newDict[circRNA_id]['soft'] += ',' + soft
                else:
                    newDict[circRNA_id]['soft'] = soft
                for sample in header:
                    if sample in newDict[circRNA_id]:
                        newDict[circRNA_id][sample] = max(newDict[circRNA_id][sample], int(mydict[circRNA_id][sample]))
                    else:
                        newDict[circRNA_id][sample] = int(mydict[circRNA_id][sample])
            new_head = header

        newLines = []
        new_header = ''

        with open(outtable, 'w') as f:
            new_header = 'circRNA_algrithm\t' + '\t'.join(new_head)
            f.write(new_header + '\n')
            for circRNA_id in newDict.keys():
                vals = []
                for sample in new_head:
                    vals.append(str(newDict[circRNA_id][sample]))
                new_line = circRNA_id + '\t' + newDict[circRNA_id]['soft'] + '\t' + '\t'.join(vals)
                f.write(new_line + '\n')
                num = rowSums(vals, 0, len(vals))
                new_line += '\t' + num

                for soft in softs:
                    soft_val = ''
                    if soft in oldLines[circRNA_id]:
                        soft_val = oldLines[circRNA_id][soft]
                    else:
                        soft_val = '\t'.join(['0' for i in range(len(header) + 1)])
                    new_line += '\t' + soft_val
                newLines.append(new_line)

                chrom, pos = circRNA_id.split(':')
                start, end = pos.split('|')

        with open(outtable + '.rowSums', 'w') as f:
            new_header += '\t' + 'merge_rowSums'
            for soft in softs:
                new_header += '\t' + oldLines['header'][soft]
            f.write(new_header + '\n')
            f.write('\n'.join(newLines) + '\n')

    def merge_bed(self, *beds):
        bedDict = defaultdict(dict)
        fo = open('all_circRNAs.bed', 'w')
        fo_hc = open('all_circRNAs.hc.bed', 'w')

        for myfile in beds:
            lines = open(myfile, 'r').read().strip().split('\n')
            for li in lines:
                chrom, start, end, soft, score, strand = li.split('\t')
                circrna_id = f'{chrom}:{start}|{end}:{strand}'
                bedDict[circrna_id][soft] = li.split('\t')
        for circrna_id in bedDict.keys():
            softs = list(bedDict[circrna_id].keys())

            softs_cat = ','.join(softs)
            chrom, start, end, soft, score, strand = bedDict[circrna_id][softs[0]]
            for so in softs:
                if so in bedDict[circrna_id]:
                    cell = bedDict[circrna_id][so]
                    if cell[-1] != strand:
                        print("strand error! " + circrna_id + ' ' + cell[-1] + ' ' + strand)
                        print(so, cell[3], soft)
                        print(softs)
            fo.write('\t'.join([chrom, start, end, circrna_id, score, strand]) + '\n')
            if len(softs) > 1:
                fo_hc.write('\t'.join([chrom, start, end, circrna_id, score, strand]) + '\n')

        fo.close()
        fo_hc.close()

    def venn(self, ciri2, ce2, outfile):
        li_ce2 = open(ce2, 'r').read().strip().split('\n')
        li_ci2 = open(ciri2, 'r').read().strip().split('\n')
        set_ci2 = ['_'.join(li.split('\t')[0:3]) for li in li_ci2]
        set_ci2 = set(set_ci2)
        fo = open(outfile, 'w')
        for li in li_ce2:
            myid = '_'.join(li.split('\t')[0:3])
            if myid in set_ci2:
                fo.write(li + '\n')
        fo.close()


if __name__ == '__main__':
    fire.Fire(ResultCircRNA)
