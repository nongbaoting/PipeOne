
process picard_RG_MD {
	//picard_AddOrReplaceReadGroups_MarkDuplicates
	tag {id}
	label 'bigMEM'

	input:
	tuple val(id), path(star2pass_sam)
	
	output:
	tuple val(id), path("${id}.dedupped.bam") 
	
	"""
	set +u; source activate pipeOne_gatk3.8; set -u
	samtools index ${star2pass_sam}
	picard -Xmx50g AddOrReplaceReadGroups \
		I=${star2pass_sam} \
		O=${id}.rg_added_sorted.bam \
		SO=coordinate RGID=${id} RGLB=mRNA RGPL=illumina RGPU=HiSeq RGSM=${id}
	samtools index ${id}.rg_added_sorted.bam
	picard -Xmx50g MarkDuplicates \
		I=${id}.rg_added_sorted.bam \
		O=${id}.dedupped.bam  \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT \
		M=${id}.metrics 

	"""
}


process SplitNCigarReads  {
	tag {id}
	label 'bigMEM'

	input:
	tuple val(id), path("dedupped.bam")
	path fasta
	path fasta_fai
	path fasta_dict
	
	output:
	tuple val(id),  path("${id}.split.bam")
	
	"""
	conda_base=`conda info --base`
	set +u; source activate pipeOne_gatk3.8; set -u
	samtools index dedupped.bam
	java -Xmx50g -jar \${conda_base}/envs/pipeOne_gatk3.8/opt/gatk-3.8/GenomeAnalysisTK.jar -T SplitNCigarReads \
		-R ${fasta} -I dedupped.bam -o ${id}.split.bam \
		-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	"""
}

process  Variant_calling {
	label 'bigMEM'
	tag {id}

	input:
	tuple val(id), path(split_bam)
	path fasta
	path fasta_fai
	path fasta_dict
	
	output:
	tuple val(id),  path("${id}.vcf"),  path("${id}.vcf.idx" )
	
	"""
	conda_base=`conda info --base`
	set +u; source activate pipeOne_gatk3.8; set -u
	samtools index ${split_bam}
	java -Xmx50g -jar \${conda_base}/envs/pipeOne_gatk3.8/opt/gatk-3.8/GenomeAnalysisTK.jar -T HaplotypeCaller \
		-R ${fasta} -I ${split_bam} \
		-dontUseSoftClippedBases -stand_call_conf 20.0 \
		-o ${id}.vcf
	"""
}

params.saveVCF = true
process  Variant_filtering {
	label 'bigMEM'
	tag {id}

	publishDir "${params.outdir_sub}/Variant_filtering/", mode: 'link', 
    saveAs: { filename ->params.saveVCF ? "$filename" : null }
	
	input:
	tuple val(id), path(input_vcf), path(vcf_idx)
	path fasta
	path fasta_fai
	path fasta_dict 
	
	output:
	tuple val(id),  path("${id}.output.vcf"),  path("${id}.output.vcf.idx")
	
	"""
	conda_base=`conda info --base`
	set +u; source activate pipeOne_gatk3.8; set -u
	java -Xmx50g -jar \${conda_base}/envs/pipeOne_gatk3.8/opt/gatk-3.8/GenomeAnalysisTK.jar -T VariantFiltration \
		-R ${fasta} \
		-V ${input_vcf} \
		-window 35 -cluster 3 \
		-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" \
		-o ${id}.output.vcf
	"""
}


workflow GATK {
    take: 
	bam
	fasta
	fasta_fai
	fasta_dict
	
    main:
    picard_RG_MD( bam )
	
	SplitNCigarReads(
		picard_RG_MD.out,
		fasta,  
		fasta_fai,
		fasta_dict )

	Variant_calling(
		SplitNCigarReads.out,
		fasta,  
		fasta_fai,
		fasta_dict )

	Variant_filtering(
		Variant_calling.out,
		fasta,  
		fasta_fai,
		fasta_dict )
	
	emit:
	id_vcf_idx = Variant_filtering.out
}