#!/usr/bin/env python3
import  os, sys

def meta_to_coding_noncoding(myfile):
    f = open(myfile,'r')
    head = next(f)
    coding,noncoding =[],[]
    for line in f:
        cell = line.strip().split('\t')
        annotation=cell[9]
        category = cell[10]
        ref_tran_id = cell[15]
        if category =='protein_coding':
            coding.extend( ref_tran_id.split(','))
        elif(category=='lncrna' and annotation=='annotated'):
            noncoding.extend(ref_tran_id.split(','))
    return (coding,noncoding)
if __name__ == '__main__':
    met = sys.argv[1]
    gtf = sys.argv[2]
    novo_lncrna = sys.argv[3]
    met_dir = os.path.dirname(met)
    novo = open(novo_lncrna, 'r').read().strip().split('\n')
    coding,noncoding = meta_to_coding_noncoding(met)
    noncoding.extend(novo)
    open( os.path.join(met_dir, 'coding.txt'), 'w').write('\n'.join(coding) + '\n')
    open(os.path.join(met_dir, 'noncoding.txt'), 'w').write('\n'.join(noncoding) + '\n')
    home = os.path.dirname( os.path.realpath(__file__))
    cmd = []
    cmd.append('python3 ' + home + '/subset_gtf_by_attr.py ' + gtf + ' transcript_id ' + met_dir + '/noncoding.txt > ' + met_dir + '/noncoding_all.gtf')
    cmd.append('cat ' + gtf + '|grep protein_coding > ' + met_dir + '/coding_all.gtf' )

    for c in cmd:
        print(c)
        os.system(c)


