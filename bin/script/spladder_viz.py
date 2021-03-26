#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : spladder_viz.py
# @Date            : 2021/03/08 19:26:56
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, random

class MYRUN:

    def cmd(self, res_dir, bam_dir, info_fi, the_id, outfi = 'out.pdf' ):
        tumor, normal = [], []
        with open(info_fi) as f:
            head = next(f).strip()
            for li in f:
                cell = li.strip().split(',')
                bam = f'{bam_dir}/{cell[0]}.bam'
                if cell[1] == '0':
                    normal.append(bam) 
                else:
                    tumor.append(bam)
        num = len(normal)
        cmd = f"spladder viz --track coverage tumor:{','.join(tumor[:]) } --track coverage adjacent_tumor:{','.join(normal) } \
            --track event {the_id} --track splicegraph ENSG00000154175 \
            -o {res_dir}  \
            -O {outfi}"
        print(cmd)

    def list_bam(self,  bam_dir, info_fi, outFi):
        tumor, normal = [], []
        with open(info_fi) as f, open(outFi, 'w') as fo:
            head = next(f).strip()
            for li in f:
                cell = li.strip().split(',')
                bam = f'{bam_dir}/{cell[0]}.bam'
                id_ = cell[0].split(".")[0]
                if not os.path.exists(bam):
                    continue
                if cell[1] == '0':
                    normal.append([id_, bam, 'Adjacent Tumor'] )
                    # fo.write('\t'.join([id_, bam, 'Adjacent Tumor']) + '\n' )
                else:
                    tumor.append([id_, bam, 'Tumor'])
                    # fo.write('\t'.join([id_, bam, 'Tumor']) + '\n' )

            # max_num = len(normal)
            # random.sample(tumor,max_num):
            for row in tumor:
                fo.write('\t'.join(row ) + '\n' )
            for row in normal:
                fo.write('\t'.join(row ) + '\n' )



if __name__ == '__main__':
    fire.Fire(MYRUN)