#!/usr/bin/env python3
import os, sys,re
from collections import defaultdict
def get_gtf_attr(gtf_line):
    if re.match('^#', gtf_line):
        return ('none ')
    chr, source, feature, start, end, score, strand, frame, atrr = gtf_line.strip().split('\t')
    atrr_dict = defaultdict(list)
    for field in atrr.strip(';').split(';'):
        field_name, field_value = field.strip(' ').split(' ')
        atrr_dict[field_name] = field_value.strip('"')
    return (atrr_dict)

if __name__ == '__main__':

    gtf = sys.argv[1]
    attr = sys.argv[2]
    info = sys.argv[3]
    info_id = open(info, 'r').read().strip().split('\n')
    info_id = set(info_id)
    with open(gtf,'rU') as f:
        for line in f:
            if re.match('^#', line):
                continue
            cell = line.strip().split('\t')
            if cell[2] == 'exon':
                attr_dict = get_gtf_attr(line)
                if attr_dict[attr] in info_id:
                    print(line, end='')

