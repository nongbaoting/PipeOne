#!/usr/bin/env nextflow

params.bam  = ""
params.sra  = ""
params.reads = ""
params.cleaned  = false
params.threads = 8
params.single = ""
params.genome =""
params.saveIntermediateVariants = false
params.saveIntermediateFiles = false
params.update_GTF = false 

params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false
params.genome_build  = params.genome ? params.genomes[ params.genome ].genome_build ?:false :false
params.annovar_data_dir = params.genome ? params.genomes[ params.genome ].annovar_data_dir ?:false :false


def read_base = params.reads.split('/')[-1]
def ifPaired = read_base =~/\{1,2\}/

scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2; scripts_3; scripts_4}

if ( params.fasta ){
    fasta = file(params.fasta)
	fasta_fai = file("${params.fasta}.fai")
	fasta_dict_path = params.fasta.replaceAll(~/\.fa.*?$/, '.dict')
	fasta_dict = file(fasta_dict_path)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}


if ( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}

if(params.annovar_data_dir){
	Channel
		.fromPath("${params.annovar_data_dir}/*")
		.ifEmpty{exit 1, "annovar data not found: ${params.annovar_data_dir}"}
		.set{annovar_db }
}

if( params.star_index ){
    star_index = Channel
        .fromPath("${params.star_index}/*")
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    star_index.into{star_index; star_index_1; star_index_2}
}


if(params.reads ){
	
	// defined single or paired reads !
	if(ifPaired ){
		ch = Channel
		.fromFilePairs(params.reads)
		.ifEmpty{error "Can not find any reads matching: ${params.reads}"}
		println "Paire-end mode"
	}else{
		Channel
			.fromPath(params.reads)
			.ifEmpty{error "Can not find any reads matching: ${params.reads}"}
			.map{file ->
				def key = file.name.toString().tokenize('.').get(0)
				return tuple(key,file)
				}
			.groupTuple()
			.set{ch}
		println "sing-end mode"
	}
	ch.into{ch;ch_print;}
	//ch_print.println()
	
	if (! params.cleaned ){
		process fastp{
			tag {id}
            maxForks 2
			
			publishDir "${params.outdir}/fastp/", mode: 'link',
				saveAs: {filename -> 
					if(filename =~ /fastp.fq.gz/) "clean/${filename}"
					else "report/${id}.${filename}" 
					}
			
			input:
			set id, file(reads) from ch
			
			output:
			set id, "${id}*.fastp.fq.gz" into read_ch
			set  "fastp.html", "fastp.json"
			
			script:
			if(!ifPaired ){
				"""
				if [ "${reads}" = "${id}.fastp.fq.gz" ];
				then
					echo "input raw reads name and output clean name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads}  -o ${id}.fastp.fq.gz  -q 20 --thread  1
				"""
			}else
			{
				"""
				if [ "${reads[0]}" = "${id}.R1.fastp.fq.gz" ];
				then
					echo "input raw read name and output clean read name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  1
				"""
			}
				
		}
		
		//fastp_out=fastp_out.map{it ->a=[it[0],it[1 ..2 ] ] }
		

	}else{
		read_ch = ch
	}
	
		
}else if(params.sra && !params.reads){
	
	Channel
	.fromPath(params.sra)
	.ifEmpty{error "Can not find any SRA matching: ${params.sra}"}
	.map{file ->
		def key = file.name.toString().tokenize('.').get(0)
		return tuple(key,file)
		}
	.groupTuple()
	.set{ sra_ch}

	sra_ch.into{sra_ch; print_ch}
	//print_ch.println()
	
	process fastq_dump_and_fastp{
		tag {id}
		
		publishDir "${params.outdir}/fastp/", mode: 'link',
			saveAs: {filename -> 
				if(filename =~ /fastp.fq.gz/) "clean/${filename}"
				else "report/${id}.${filename}" 
				}
		
		input:
		set id, file(sra) from sra_ch
		
		output:
		set id, "${id}*fastp.fq.gz" into read_ch
		set "fastp.html", "fastp.json"
		
		script:
		if(!ifPaired ){
			"""
			fastq-dump --split-3  $sra
			fastp -i ${id}_1.fastq  -o ${id}.fastp.fq.gz -O -q 20 --thread  1
			"""
		}else{
			"""
			fastq-dump --split-3  $sra
			fastp -i ${id}_1.fastq -I ${id}_2.fastq -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  1
			"""
		}
		
	}
	
	//read_ch=fastp_out.map{it ->a=[it[0],it[1 .. 2]] }
    }else if(params.bam ){
	
	Channel
		.fromPath(params.bam)
		.ifEmpty{error "Can not find any bam file matching: ${params.bam}"}
		.map{file ->
			def key = file.name.toString().tokenize('.').get(0)
			return tuple(key,file)
			}
		.groupTuple()
		.set{ star_2pass_out }
			
	
	}else{
		
	exit 1, "you should input --reads or --sra or --bam"
}




if(params.reads || params.sra ){
	
	read_ch.into{ reads_star; read_star_2pass; reads_star_fusion; reads_print}
	
	
	process star_1_pass{
		tag {id}
		
		input:
		file "STARIndex/*" from  star_index.collect()
		set  id, file(reads) from reads_star
		file gtf
		file fasta
		
		output:
		file "${id}.SJ.out.tab" into star_1_out
		
		"""
		set +u; source activate pipeOne_py3; set -u
		STAR --runThreadN ${params.threads} --genomeDir STARIndex \
			--readFilesIn ${reads}  --readFilesCommand zcat \
			--outFileNamePrefix	${id}.
			rm ${id}.Aligned.out.sam 
		"""
	}


	process starIndex_for_2pass{
		
		input:
		file "star_junction/*" from star_1_out.collect()
		file gtf
		file fasta
		
		output:
		
		file "STARIndex_2pass/*" into STARIndex_2pass
		
		"""
		set +u; source activate pipeOne_py3; set -u
		find star_junction -type f |xargs -i cat {} > SJ.out.tab
		
		genomeDir=STARIndex_2pass
		mkdir \$genomeDir
		STAR --runMode genomeGenerate --genomeDir \$genomeDir --genomeFastaFiles ${fasta} \
			 --sjdbFileChrStartEnd SJ.out.tab  --runThreadN  ${params.threads}
		"""
	}


	process star_2_pass{
		
		input:
		file "STARIndex_2pass/*" from  STARIndex_2pass.collect()
		set  id, file(reads) from read_star_2pass
		file gtf
		file fasta
		
		output:
		set id, file("${id}.Aligned.out.sam") into star_2pass_out
		
		"""
		set +u; source activate pipeOne_py3; set -u
		STAR --genomeDir STARIndex_2pass --runThreadN  ${params.threads} \
			--readFilesIn ${reads} --readFilesCommand zcat \
			--outFileNamePrefix	${id}.
		"""
	}


}


process picard_AddOrReplaceReadGroups_MarkDuplicates {
	
	input:
	set id, file(star2pass_sam) from star_2pass_out
	
	output:
	set id, "${id}.dedupped.bam" into picard_out
	
	"""
	set +u; source activate pipeOne_gatk3.8; set -u
	picard AddOrReplaceReadGroups \
		I=${star2pass_sam} \
		O=${id}.rg_added_sorted.bam \
		SO=coordinate RGID=${id} RGLB=mRNA RGPL=illumina RGPU=HiSeq RGSM=${id}
	
	picard MarkDuplicates \
		I=${id}.rg_added_sorted.bam \
		O=${id}.dedupped.bam  \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT \
		M=${id}.metrics 

	"""
}


process SplitNCigarReads  {
	
	input:
	set id, "dedupped.bam" from picard_out
	file fasta
	file fasta_fai
	file fasta_dict
	
	
	output:
	set id, "${id}.split.bam" into splitN_out
	
	"""
	conda_base=`conda info --base`
	set +u; source activate pipeOne_gatk3.8; set -u
	java -jar \${conda_base}/envs/pipeOne_gatk3.8/opt/gatk-3.8/GenomeAnalysisTK.jar -T SplitNCigarReads \
		-R ${fasta} -I dedupped.bam -o ${id}.split.bam \
		-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	"""
}

process  Variant_calling {
	
	input:
	set id, file(split_bam) from splitN_out
	file fasta
	file fasta_fai
	file fasta_dict
	
	output:
	set id, "${id}.vcf", "${id}.vcf.idx" into vcf_out
	
	"""
	conda_base=`conda info --base`
	set +u; source activate pipeOne_gatk3.8; set -u
	java -jar \${conda_base}/envs/pipeOne_gatk3.8/opt/gatk-3.8/GenomeAnalysisTK.jar -T HaplotypeCaller \
		-R ${fasta} -I ${split_bam} \
		-dontUseSoftClippedBases -stand_call_conf 20.0 \
		-o ${id}.vcf
	"""
}

process  Variant_filtering {
	publishDir "${params.outdir}/Variant_filtering/", mode: 'link', 
    saveAs: { filename -> params.saveIntermediateVariants ? "$filename" : null }
	
	input:
	set id, file(input_vcf), file(vcf_idx) from vcf_out
	file fasta
	file fasta_fai
	file fasta_dict 
	
	output:
	tuple id, "${id}.output.vcf", "${id}.output.vcf.idx" into vcf_filter_out_annovar
	
	"""
	conda_base=`conda info --base`
	set +u; source activate pipeOne_gatk3.8; set -u
	java -jar \${conda_base}/envs/pipeOne_gatk3.8/opt/gatk-3.8/GenomeAnalysisTK.jar -T VariantFiltration \
		-R ${fasta} \
		-V ${input_vcf} \
		-window 35 -cluster 3 \
		-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" \
		-o ${id}.output.vcf
	"""
}



gene_based='refGene'
//filter_based='cytoBand,exac03,avsnp147,dbnsfp30a'
operation='g'
//-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f 

process variant_AnnotateAnnovar {
	
	publishDir "${params.outdir}/variantAnnotateAnnovar/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateVariants ? "$filename" : null }

    input:
    tuple id, "${id}.output.vcf", "${id}.output.vcf.idx" from  vcf_filter_out_annovar
	file "annovar_db/*" from annovar_db.collect()
	
    output:
    file "${id}.${params.genome_build}_multianno.txt" into annovar_out
    
    
    """
	set +u; source activate pipeOne_gatk3.8; set -u
	perl ${params.annovar_BinDir}/convert2annovar.pl -format vcf4 -filter PASS --allsample ${id}.output.vcf --outfile avinput 1> conv.log 2>&1
	perl ${params.annovar_BinDir}/table_annovar.pl avinput.${id}.avinput  annovar_db \
		--buildver ${params.genome_build} --remove --protocol $gene_based \
		--operation $operation -nastring . \
		-out ${id} --thread ${params.threads} 1> SnpIndel_annovar.log 2>&1
    """
}


process To_gene_base_table {
	publishDir "${params.outdir}/annovar_table/", mode: 'link'

	input:
	file "annovar_res/*" from annovar_out.collect()

	output:
	file "snp.geneBase.tsv"

	"""
	set +u; source activate pipeOne_py3; set -u
	python3 ${baseDir}/bin/SNP.py gene_base_table snp.geneBase.tsv annovar_res
	"""
}