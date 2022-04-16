

// to get known annotate genes/lncRNA and novel lncRNAs
process gffcompare{
	publishDir "${params.outdir_sub}/gffcompare/", mode: 'link',
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	path fasta
	path fasta_fai
	path gtf
	path lncpedia_gtf
	path "assembly.gtf" 
	
	
	output:
	tuple  path("novel_lncRNA_candidate.gtf"), path("novel_lncRNA_candidate.fa") , emit: nov_lnc_candidate
	path "novel_lncRNA_candidate.fa" , emit: nov_lnc_candidate_fa
	path "*"
	path "gffcmp.assembly.gtf.tmap" , emit: gffcompare_tmap
	
	"""
	set +u; source activate pipeOne_py3; set -u
	#!/bin/bash
	cat ${gtf} ${lncpedia_gtf} >all.gtf
	gffcompare -r all.gtf assembly.gtf 

	## get multi-exons lncRNAs or single exon length >= 2000	
	python3 ${baseDir}/bin/RNAseq/lncRNA.py gffcompare_tmap gffcmp.assembly.gtf.tmap novel_lncRNA_candidate.list ${params.novel_lnc_single_exon_len}
	
	python3 ${baseDir}/bin/RNAseq/gtf.py get_by_trans_id assembly.gtf novel_lncRNA_candidate.list novel_lncRNA_candidate.gtf
	gffread novel_lncRNA_candidate.gtf -g ${fasta} -w novel_lncRNA_candidate.fa -W
	"""
}

process prepare_ref_gtf {
    publishDir "${params.outdir_sub}/reference_gtf_info/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles  ? "$filename" : null }
	
	input:
	path "gencode_annotation.gtf" 
	path "lncpedia.gtf" 
	path genecode_lncRNA_gtf
    
    output:
    tuple path("gencode_protein_coding.gtf"), path("known_lncRNA.gtf"), path("known_lncRNA.list") , emit: known_gtf_lst
    path "known_lncRNA.gtf" , emit: known_lncRNA_ch
    tuple path("known_lncRNA.list"), path("nonOverlap_lncpedia.list"),path("gencode_protein_coding_tran_id.list")
    
    """
    #!/bin/bash
	set +u; source activate pipeOne_py3; set -u
	gffcompare -r ${genecode_lncRNA_gtf} lncpedia.gtf 
	awk '\$3 =="x"||\$3=="u"||\$3=="i"{print}' gffcmp.lncpedia.gtf.tmap > nonOverlap.gtf.tmap    
	cut -f5 nonOverlap.gtf.tmap > nonOverlap_lncpedia.list
	python3 ${baseDir}/bin/RNAseq/gtf.py get_by_trans_id lncpedia.gtf nonOverlap_lncpedia.list  nonOverlap_lncpedia.gtf
    cat gencode_annotation.gtf  |grep "protein_coding" > gencode_protein_coding.gtf
    cat ${genecode_lncRNA_gtf} nonOverlap_lncpedia.gtf > known_lncRNA.gtf
    python3 ${baseDir}/bin/RNAseq/gtf.py get_id gencode_protein_coding.gtf transcript_id gencode_protein_coding_tran_id.list
    python3 ${baseDir}/bin/RNAseq/gtf.py get_id known_lncRNA.gtf transcript_id known_lncRNA.list
    """
}


process filter_coding_potentail{

    publishDir "${params.outdir_sub}/annotations_and_fasta/", mode: 'link'
		

	publishDir "${params.outdir_sub}/novel_lncRNA/", mode: 'link',
		saveAs: { filename-> 
			if( filename =~ /novel_lncRNA/) filename
			}
	
	publishDir "${params.outdir_sub}/novel_lncRNA/tmp/", mode: 'link',
		saveAs: { filename ->  params.saveIntermediateFiles  ? "$filename" : null }

	publishDir "${params.outdir_sub}/coding_potential/", mode: 'link',
		saveAs: { filename ->  filename ==  "coding_potential_sum.tsv" ? "$filename" : null }

			
    input:
    path gtf
	path fasta
	path fasta_fai
	path "CPAT.out" 
    path "PLEK.out" 
	path "CPPred.out"
    tuple path("gencode_protein_coding.gtf"), path("known_lncRNA.gtf"),path("known_lncRNA.list") 
	tuple path("novel_lncRNA_candidate.gtf"), path("novel_lncRNA_candidate.fa") 
   
    
    output:
    path "protein_coding_and_all_lncRNA.gtf" , emit: cal_expr_gtf
    path "protein_coding_and_all_lncRNA.fa" , emit: cal_expr_fa
    tuple path("all_lncRNA.gtf"),path("all_lncRNA.list"),path("novel_lncRNA.fa"),path("novel_lncRNA.gtf"), path("novel_lncRNA.list") , emit: nov_lst
	path "all_lncRNA.gtf" , emit: all_lnc_gtf
    path "TUCP*"
	path "all_lncRNA.list" , emit: all_lnc_lst
	path "coding_potential_sum.tsv"
	path "novel_lncRNA.gtf" , emit: nov_lnc_gtf

    """
	set +u; source activate pipeOne_py3; set -u;
	echo "filter coding potential!"
	python3 ${baseDir}/bin/RNAseq/gtf.py to_info novel_lncRNA_candidate.gtf novel_lncRNA_candidate.info
	python3 ${baseDir}/bin/RNAseq/lncRNA.py combine_coding_prediction novel_lncRNA_candidate.info CPAT.out PLEK.out CPPred.out ${params.species}  novel_lncRNA.list
	
	# all non-coding RNA
	python3 ${baseDir}/bin/RNAseq/get_seq_by_id.py novel_lncRNA_candidate.fa novel_lncRNA.list novel_lncRNA.fa
	python3 ${baseDir}/bin/RNAseq/gtf.py get_by_trans_id novel_lncRNA_candidate.gtf novel_lncRNA.list  novel_lncRNA.gtf
	
	# TUCP
	python3 ${baseDir}/bin/RNAseq/get_seq_by_id.py novel_lncRNA_candidate.fa TUCP.list TUCP.fa
	python3 ${baseDir}/bin/RNAseq/gtf.py get_by_trans_id novel_lncRNA_candidate.gtf TUCP.list  TUCP.gtf
	
	# for expression calculate
	cat novel_lncRNA.list known_lncRNA.list >all_lncRNA.list
	cat novel_lncRNA.gtf  known_lncRNA.gtf  >all_lncRNA.gtf
    cat gencode_protein_coding.gtf novel_lncRNA.gtf known_lncRNA.gtf > protein_coding_and_all_lncRNA.tmp.gtf
	python3 ${baseDir}/bin/RNAseq/gtf.py add_gene_name protein_coding_and_all_lncRNA.tmp.gtf protein_coding_and_all_lncRNA.gtf
	gffread all_lncRNA.gtf -g ${fasta} -w all_lncRNA.fa -W
	gffread protein_coding_and_all_lncRNA.gtf -g ${fasta} -w protein_coding_and_all_lncRNA.fa -W
    """
}

process format_lncRNA_info{
	publishDir "${params.outdir_sub}/novel_lncRNA/", mode: 'link'

	input:
	
	tuple path("all_lncRNA.gtf"),path("all_lncRNA.list"),path("novel_lncRNA.fa"),path("novel_lncRNA.gtf"), path("novel_lncRNA.list")
	tuple path("gencode_protein_coding.gtf"), path("known_lncRNA.gtf"), path("known_lncRNA.list")
	path "protein_coding_and_all_lncRNA.gtf"
	
	output:
	tuple path("all_lncRNA_info.tsv"), path("novel_lncRNA_info.tsv") , emit: lncRNA_res
	path "protein_coding_and_all_lncRNA.info.tsv" 
	path "all_lncRNA_info.tsv", emit: all_lnc_info
	
	shell:
	'''
	set +u; source activate pipeOne_py3; set -u
    gffcompare -r gencode_protein_coding.gtf all_lncRNA.gtf
	
	#conver to bed and colsest to bed
	awk '$3=="exon"{print}'  gencode_protein_coding.gtf >gencode_protein_coding.exon.gtf
	awk '$3=="exon"{print}'  all_lncRNA.gtf             >all_lncRNA.exon.gtf
	gtf2bed < gencode_protein_coding.exon.gtf |cut -f1-6 >gencode_protein_coding.exon.gtf.bed
	gtf2bed < all_lncRNA.exon.gtf |cut -f1-6 >all_lncRNA.exon.gtf.bed
	bedtools closest -a all_lncRNA.exon.gtf.bed -b gencode_protein_coding.exon.gtf.bed  -d |sort -u > cloest.txt
	
	#gtf to info
	python3  !{baseDir}/bin/RNAseq/gtf.py to-info  all_lncRNA.gtf all_lncRNA.info
	python3  !{baseDir}/bin/RNAseq/gtf.py to-info  gencode_protein_coding.gtf gencode_protein_coding.info.tsv
	python3  !{baseDir}/bin/RNAseq/gtf.py to-info  protein_coding_and_all_lncRNA.gtf protein_coding_and_all_lncRNA.info.tsv
	
	#add info
	python3  !{baseDir}/bin/RNAseq/gene_info_add.py add-gffcompre        all_lncRNA.info        gffcmp.all_lncRNA.gtf.tmap all_lncRNA.info.gffcmp
	python3  !{baseDir}/bin/RNAseq/gene_info_add.py add-bedtool-closest  all_lncRNA.info.gffcmp cloest.txt                 all_lncRNA_info.tsv
	python3  !{baseDir}/bin/RNAseq/gene_info_add.py add-geneName all_lncRNA_info.tsv gencode_protein_coding.info.tsv all_lncRNA_info.tsv2
	mv all_lncRNA_info.tsv2 all_lncRNA_info.tsv
	cat all_lncRNA_info.tsv |awk 'NR==1 || $1 ~/^TU/' >novel_lncRNA_info.tsv
	'''
}


process classify_lncRNA {

    publishDir "${params.outdir_sub}/novel_lncRNA/classify", mode: 'link'
   
    input:
    path "all_lncRNA.gtf"
    // actually, I only need "gencode_protein_coding.gtf"
    tuple path("gencode_protein_coding.gtf"), path("known_lncRNA.gtf"), path("known_lncRNA.list") 
    
    output:
    tuple path("lncRNA_class_closest_PCG.tsv"), path("novel_lncRNA_info.tsv")
    
	shell:
    '''
	set +u; source activate pipeOne_py3; set -u
    awk '$3=="exon" && $1 ~/^chr|^\\d/' all_lncRNA.gtf  > all_lncRNA.exon.gtf
    awk '$3=="exon"' gencode_protein_coding.gtf > gencode_protein_coding.exon.gtf
    
    gtf2bed < all_lncRNA.exon.gtf |cut -f1-6 > novel_final.bed.a
    gtf2bed < gencode_protein_coding.exon.gtf |cut -f1-6 >annot.bed.b
    sort -k1,1 -k2,2n novel_final.bed.a >novel_final.bed
    sort -k1,1 -k2,2n annot.bed.b >annot.bed

    python3  !{baseDir}/bin/RNAseq/gtf.py  to-gene-bed gencode_protein_coding.exon.gtf annot_gene.bed.a
    python3  !{baseDir}/bin/RNAseq/gtf.py  to-gene-bed all_lncRNA.exon.gtf  novel_final_gene.bed.b

    sort -k1,1 -k2,2n annot_gene.bed.a > annot_gene.bed
    sort -k1,1 -k2,2n novel_final_gene.bed.b >novel_final_gene.bed

    bedtools intersect -a novel_final_gene.bed -b annot_gene.bed -wa -wb -v  >non_overlaps.txt
    bedtools intersect -a novel_final_gene.bed -b annot_gene.bed -wa -wb     >overlaps.txt

    awk '$6==$12 {print}' FS="\t" OFS="\t" overlaps.txt >sense_or_intron.txt
    awk '$6!=$12 {print}' FS="\t" OFS="\t" overlaps.txt >antisense.txt

    python3  !{baseDir}/bin/RNAseq/venn.py venn novel_final_gene.bed sense_or_intron.txt novel_final_sense_intron.bed tmp  3 3

    bedtools intersect -a novel_final_sense_intron.bed -b annot.bed -wa -wb -v  >intron.txt
    bedtools intersect -a novel_final_sense_intron.bed -b annot.bed -wa -wb     >sense.txt

    # -D b: Report distance with respect to the orientation of the interval in B. 
    # That is, when B is on the - strand, “upstream” means A has higher start/stop coordinates. 
    # When B is on the + strand, “upstream” means A has lower start/stop coordinates.
    bedtools closest -a novel_final_gene.bed -b annot.bed -D b > closest.txt

    # -s Require same strandedness. By default, overlaps are reported without respect to strand.
    bedtools closest -a novel_final_gene.bed -b annot.bed -s -D b  >closest_sense.txt

    awk '{print $4"\tAntisense" }' antisense.txt|sort -u >lncRNA_classification.txt
    awk '{print $4"\tIntron" }' 	intron.txt |sort -u>>lncRNA_classification.txt
    awk '{print $4"\tSense" }' 	sense.txt |sort -u>>lncRNA_classification.txt
    awk '{print $4"\tIntergenic" }' non_overlaps.txt |sort -u>>lncRNA_classification.txt

    ## get info
    python3  !{baseDir}/bin/RNAseq/gtf.py to-info  all_lncRNA.exon.gtf all_lncRNA.info
	python3  !{baseDir}/bin/RNAseq/gtf.py to-info  gencode_protein_coding.gtf gencode_protein_coding.info.tsv
   
    ## output "lncRNA_class_closest_PCG.tsv"
    python3  !{baseDir}/bin/RNAseq/lncRNA_classify.py final lncRNA_classification.txt \
    closest.txt closest_sense.txt all_lncRNA.info gencode_protein_coding.info.tsv
    
    cat lncRNA_class_closest_PCG.tsv | awk 'NR==1 || $1 ~/^TU/' > novel_lncRNA_info.tsv
    '''
}


process epxr_gene_summary{

		publishDir "${params.outdir_sub}/novel_lncRNA/expressed", mode: 'link',
		saveAs: { filename-> 
			if( filename =~ /novel_lncRNA/) filename
			}

		input:
		path "gene_id.keep.txt" 
		path "protein_coding_and_all_lncRNA.gtf" 
		path "novel_lncRNA.gtf" 

		output:
		tuple path("protein_coding_and_all_lncRNA.info.expressed.tsv"), path("novel_lncRNA.info.expressed.tsv")
		tuple path("protein_coding_and_all_lncRNA.expressed.gtf"), path("novel_lncRNA.expressed.gtf")
		path "protein_coding_and_all_lncRNA.expressed.gtf", emit: gtf

		"""
		set +u; source activate pipeOne_py3; set -u
		python3 ${baseDir}/bin/RNAseq/gtf.py  get_by_gene_id --gtf novel_lncRNA.gtf --info gene_id.keep.txt --out novel_lncRNA.expressed.gtf
		python3 ${baseDir}/bin/RNAseq/gtf.py  get_by_gene_id --gtf protein_coding_and_all_lncRNA.gtf --info gene_id.keep.txt --out protein_coding_and_all_lncRNA.expressed.gtf
		python3 ${baseDir}/bin/RNAseq/gtf.py  to-info  --gtf protein_coding_and_all_lncRNA.expressed.gtf --out_tsv protein_coding_and_all_lncRNA.info.expressed.tsv
		python3 ${baseDir}/bin/RNAseq/gtf.py  to-info  --gtf novel_lncRNA.expressed.gtf --out_tsv novel_lncRNA.info.expressed.tsv
		"""

	}


process  sep_lnc_mRNA  {
    input:
    path tpm
    path lncinfo

    output:
    path "lncR_gene.tpm.tsv", emit: lnc_TPM
	path "prot_gene.tpm.tsv", emit: prot_TPM

    script:
    """
    python3 ${baseDir}/bin/RNAseq/summary_table.py sep_lnc_mRNA \\
        $tpm \\
        $lncinfo
    """

}