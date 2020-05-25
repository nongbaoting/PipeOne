#!/usr/bin/env python3
__author__ = "Baoting Nong"
__data__   = "18/09/15"

import fire
from Bio import SeqIO
from collections import defaultdict

class Fasta:

    def get_seq_by_id(self, fa, info, out, id_col =0):
        """
        :param fa: input fasta file
        :param info: id file in tab delimeter
        :param out: output fasta file
        :param id_col: id column number, default 0
        """
        lines = open(info, 'r').read().strip().split('\n')
        ids = [line.split('\t')[id_col] for line in lines]
        fo = open(out, 'w')
        for record in SeqIO.parse(fa, 'fasta'):
            if record.id in ids:
                SeqIO.write(record, fo, 'fasta')
        fo.close()

if __name__ == '__main__':
    fire.Fire(Fasta)