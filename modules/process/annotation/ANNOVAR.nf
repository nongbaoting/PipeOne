
process annovar_sprint{
    publishDir "${params.outdir_sub}/sprint/ANNOVAR", mode: 'link'

    input:
    tuple path("SPRINT_all.tsv"), path("SPRINT_A2I.tsv"), path("SPrint_A2I_table.tsv")
    path "annovar_db/*"

    output:
    tuple path("SPRINT_all.annovar.tsv"), 
          path("SPRINT_A2I.annovar.tsv")
		  
    path  "SPrint_A2I_table.annovar.tsv", emit: table

    """
    set +u; source activate pipeOne_py3; set -u
    # annovar annotate rnaEditing
	cut -f1-4 SPRINT_all.tsv |awk 'NR !=1 {split(\$4, a, ""); print \$1,\$3,\$3, a[1],a[2] }' OFS="\t" > rnaEditing_pos.txt
	perl ${params.annovar_BinDir}/annotate_variation.pl -geneanno -dbtype refGene -buildver hg38 rnaEditing_pos.txt annovar_db
	python3 ${baseDir}/bin/RNAseq/rnaEditing.py add_tab_anno rnaEditing_pos.txt.variant_function SPrint_A2I_table.tsv SPrint_A2I_table.annovar.tsv
	python3 ${baseDir}/bin/RNAseq/rnaEditing.py add_res_anno rnaEditing_pos.txt.variant_function SPRINT_all.tsv SPRINT_all.annovar.tsv
	python3 ${baseDir}/bin/RNAseq/rnaEditing.py add_res_anno rnaEditing_pos.txt.variant_function SPRINT_A2I.tsv SPRINT_A2I.annovar.tsv
    """


}


def gene_based='refGene'
//filter_based='cytoBand,exac03,avsnp147,dbnsfp30a'
def operation='g'
//-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f 
process Annovar_vcf {
	errorStrategy 'ignore'
	label 'bigCPU'
	
	publishDir "${params.outdir_sub}/variantAnnotateAnnovar/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }

    input:
    tuple val(id), path("${id}.output.vcf"), path("${id}.output.vcf.idx" )
	path "annovar_db/*" 
	
    output:
    path "${id}.${params.genome_build}_multianno.txt" 
    
    """
	set +u; source activate pipeOne_gatk3.8; set -u
	perl ${params.annovar_BinDir}/convert2annovar.pl -format vcf4 -filter PASS --allsample ${id}.output.vcf --outfile avinput 1> conv.log 2>&1
	perl ${params.annovar_BinDir}/table_annovar.pl avinput.${id}.avinput  annovar_db \
		--buildver ${params.genome_build} --remove --protocol $gene_based \
		--operation $operation -nastring . \
		-out ${id} --thread ${task.cpus} 1> SnpIndel_annovar.log 2>&1
    """
}

process Annovar_vcf_2_genebaseTable {
	
	publishDir "${params.outdir_sub}/annovar_table/", mode: 'link'

	input:
	path "annovar_res/*" 

	output:
	path "snp.geneBase.tsv"

	"""
	set +u; source activate pipeOne_py3; set -u
	python3 ${baseDir}/bin/RNAseq/SNP.py gene_base_table snp.geneBase.tsv annovar_res
	"""
}

workflow ANNOVAR_VCF {
	take:
	vcf
	annovar_data

	main:
	Annovar_vcf(vcf, annovar_data ) 
	Annovar_vcf_2_genebaseTable( Annovar_vcf.out.collect() )

	emit:
	snp_geneBase = Annovar_vcf_2_genebaseTable.out

}
