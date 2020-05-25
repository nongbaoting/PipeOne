#!/usr/bin/env python3
import fire
from Bio import SeqIO
from collections import defaultdict
def extract(fa,info,out):

    lines = open(info, 'r').read().strip().split('\n')
    ids = [line.split(' ')[0]  for line in lines]
    fo = open(out, 'w')
    for record in SeqIO.parse(fa, 'fasta'):
        if record.id in ids:
            SeqIO.write(record,fo,'fasta')
    fo.close()
if __name__ == '__main__':
    fire.Fire(extract)
