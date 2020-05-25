import fire
from collections import defaultdict
import intervals

def annotate_type(strand_q, strand):
    if strand_q == strand:
        return 'cis-regulatory'
    else:
        return 'bidirectional'
def closest2bed(lines):
    clsExonDict = defaultdict(list)
    for li in lines:
        cell = li.split('\t')

        qry_id = cell[3]
        pro_id = cell[9]
        strand_q = cell[5]
        strand = cell[11]
        dist = int(cell[12])
        if dist == 0:
            start_lnc, end_lnc = int(cell[1]), int(cell[2])
            start_pro, end_pro = int(cell[7]), int(cell[8])
            region_1 = intervals.openclosed(start_lnc, end_lnc)
            region_2 = intervals.openclosed(start_pro, end_pro)
            if not region_1 & region_2:
                dist =1
        if qry_id in clsExonDict:
            if dist < clsExonDict[qry_id][1]:
                clsExonDict[qry_id] = [pro_id, dist, strand_q, strand]
        else:
            clsExonDict[qry_id] = [pro_id, dist, strand_q, strand]
    return clsExonDict

class Info:

    def add_gffcompre(self, info, gffcompre, out):
        gff_line = open(gffcompre, 'r').read().strip().split('\n')
        gffDict = defaultdict(list)
        cc = {
            '=': 'complete_match',
            'i': 'intronic',
            'u': 'intergenic',
            'x': 'antisense'
        }
        for li in gff_line[1:]:
            cell = li.split('\t')[0:6]
            qry_id = cell[4]
            gffDict[qry_id] = cell

        info_line = open(info, 'r').read().strip().split('\n')
        header = info_line[0]
        new_lines = [ header + '\t' + '\t'.join(['lncRNA_type','closest_PCT','closest_PCG','distance_PCG', 'strand_PCT']) ]
        uncomapre = 0
        for li in info_line[1:]:
            cell = li.split('\t')
            tx_id = cell[0]
            lnc_type, ref_id, ref_gene_id,dist, strand ='NA','NA','NA','NA','NA'
            if not tx_id in gffDict:
                print(tx_id)
                uncomapre +=1
            else:
                ref_gene_id, ref_id, class_code, qry_gene_id, qry_id, num_exons = gffDict[ tx_id ]

            if class_code in cc:
                lnc_type = cc[class_code]
                if class_code in ['=','i','x']:
                    dist = '0'
            else:
                lnc_type = class_code

            add_info = [lnc_type,ref_id,ref_gene_id,dist]

            new_lines.append( li + '\t' + '\t'.join(add_info))

        with open(out, 'w') as f:
            f.write('\n'.join(new_lines) + '\n')

    def add_bedtool_closest(self,info,closest,out):
        info_line = open(info, 'r').read().strip().split('\n')
        cls = open(closest, 'r').read().strip().split('\n')
        clsDict = defaultdict(list)
        for li in cls:
            cell = li.split('\t')
            qry_id = cell[3]
            pro_id = cell[9]
            strand_q = cell[5]
            strand = cell[11]
            dist = int( cell[12] )
            if qry_id in clsDict:
                if dist > clsDict[qry_id][1]:
                    dist = clsDict[qry_id][1]
                    strand = clsDict[qry_id][3]
            clsDict[qry_id] =[pro_id,dist,strand_q,strand]
        print(len(clsDict))
        new_line = [info_line[0]]

        test_type = ['antisense', 'intronic']
        for li in info_line[1:]:
            cl = li.split('\t')
            gene_id = cl[1]
            new_l = li

            if gene_id in clsDict:
                pro_id, dist, strand_q, strand = clsDict[ gene_id ]
                if dist == 0 and strand != strand_q:
                    cl[10] = 'antisense'
                elif  0 < dist <= 2000 and cl[ 10 ] not in test_type :
                    cl[10] = annotate_type(strand_q, strand)
                elif dist >2000:
                    cl[10] = 'intergenic'
                cl[12] = pro_id
                cl[13] = str(dist)
                cl.append(strand)
            new_l = '\t'.join(cl)
            new_line.append(new_l)
        with open(out, 'w') as f:
            f.write('\n'.join(new_line) +'\n')

    def add_2_bedtool_closest(self,info,closest_gene, closest_exon,out):
        info_line = open(info, 'r').read().strip().split('\n')
        clsExon = open(closest_exon, 'r').read().strip().split('\n')
        clsGene = open(closest_gene, 'r').read().strip().split('\n')
        clsExonDict = closest2bed(lines = clsExon)
        clsGeneDict = closest2bed(lines = clsGene)

        new_line = [info_line[0]]

        for li in info_line[1:]:
            cl = li.split('\t')
            gene_id = cl[1]
            new_l = li
            if gene_id in clsGeneDict:
                pro_id, dist, strand_q, strand = clsGeneDict[ gene_id ]
                if dist == 0 and strand != strand_q:
                    cl[10] = 'antisense'
                elif dist == 0 :
                    pro_id_exon, dist_exon, strand_q_exon, strand_exon = clsExonDict[gene_id]
                    if strand_exon != strand_q_exon:
                        cl[10] = 'antisense'
                    elif dist_exon > 0 :
                        cl[10] = 'intronic'
                    else:
                        cl[10] = "sense_overlapping"
                        print('distance error ! ',gene_id, clsExonDict[gene_id], clsGeneDict[gene_id], sep='\t')
                elif  0 < dist <= 2000 :
                    cl[10] = annotate_type(strand_q, strand)
                elif dist >2000:
                    cl[10] = 'intergenic'
                cl[12] = pro_id
                cl[13] = str(dist)
                cl.append(strand)
            else:
                print('not found distance value: ' + gene_id)
            new_l = '\t'.join(cl)
            new_line.append(new_l)
        with open(out, 'w') as f:
            f.write('\n'.join(new_line) +'\n')

    def add_geneName(self,info,info_protein,out):
        info_line    = open(info, 'r').read().strip().split('\n')
        protein_line = open(info_protein, 'r').read().strip().split('\n')
        proDict = defaultdict(list)
        for li in protein_line[1:]:
            cell = li.split('\t')
            proDict[ cell[1] ] = cell[2]

        new_line = [ info_line[0] + '\tclosest_PCG_Name']
        for li in info_line[1:]:
            cl = li.split('\t')
            gene_id = cl[12]
            geneName = gene_id
            if gene_id in proDict:
                geneName = proDict[gene_id]
            new_l = li + '\t' + geneName
            new_line.append(new_l)
        with open(out, 'w') as f:
            f.write('\n'.join(new_line) + '\n')
    def info_and_expr(self,info,expr,out,col = 0):
        info_line = open(info, 'r').read().strip().split('\n')
        expr_line = open(expr, 'r').read().strip().split('\n')
        expDict = defaultdict(list)
        expHead = '\t'.join( expr_line[0].split('\t')[1:] )
        for li in expr_line[1:]:
            cell = li.split('\t')
            expDict[ cell[col] ] = '\t'.join( cell[1:] )

        new_line = [info_line[0] + '\t' + expHead]
        for li in info_line[1:]:
            cl = li.split('\t')
            gene_id = cl[1]
            if gene_id in expDict:
                new_l = li + '\t' +expDict[gene_id]
                new_line.append(new_l)
        with open(out, 'w') as f:
            f.write('\n'.join(new_line) + '\n')

if __name__ == '__main__':
    fire.Fire(Info)