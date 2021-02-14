#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/12/19 13:18
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire,re, os
from  collections import defaultdict

class RUN:

    def merge(self, out, mydir, pattern = ".telescope_report.tsv"):
        mydict = defaultdict(dict)
        ss = []
        reg_pat = re.compile(pattern)
        for entry in os.scandir(mydir):
            if os.path.isfile(entry.path):
                sample = reg_pat.sub('', entry.name )
                print(entry.path)
                print(entry.name)
                print(sample)
                with open(entry.path, 'r') as f:
                    for li in f:
                        if re.match("##|transcript", li):
                            continue
                        cell = li.strip().split('\t')
                        tx_id = cell[0]
                        mydict[tx_id][sample] = cell[2]
                ss.append(sample)

        with open(out, 'w') as fo:
            ss_sorted = sorted(ss)
            fo.write('tx_id\t' + '\t'.join(ss_sorted) + '\n')
            for tx_id in mydict:
                items = [ tx_id ]
                for sample in ss_sorted:
                    if sample in mydict[tx_id]:
                        items.append( mydict[tx_id][sample] )
                    else:
                        items.append('0')
                fo.write('\t'.join(items) + '\n')

if __name__ == '__main__':
    fire.Fire(RUN)