#!/usr/bin/env python3
from collections import defaultdict
import fire,os,re

class Table:


    def sailfish(self, outfile,col , *files):
        read_table = defaultdict(dict)
        all_sample = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            for line in lines[1:]:
                cell = line.split('\t')
                name, Length, EffectiveLength, TPM, numReads = cell
                val = cell[col]
                read_table[name][sample] = val
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


    def kallisto(self, outfile, col , *files):
        read_table = defaultdict(dict)
        all_sample = []
        for myfile in files:
            lines = open(myfile, 'r').read().strip().split('\n')
            sample = os.path.basename(myfile).split('.')[0]
            all_sample.append(sample)
            for line in lines[1:]:
                cell = line.split('\t')
                target_id,length,eff_length,est_counts,tpm = cell
                val = cell[col]
                read_table[target_id][sample] = val
        out_lines = []
        #header = 'target_id\t' + '\t'.join(all_sample)
        header =  '\t'.join(all_sample)
        out_lines.append(header)
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



if __name__ == '__main__':
    fire.Fire(Table)
