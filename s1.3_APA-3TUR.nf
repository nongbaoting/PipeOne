#!/usr/bin/env nextflow

params.genome = ""
params.read  = ""
params.sra    = ""
params.single = false
params.cleaned  = false
params.datainfo = false
params.contrast = false
params.saveIntermediateFiles = false

params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.utr_gtf =  params.genome ? params.genomes[ params.genome ].utr_gtf  ?: false : false

params.salmon_index_3UTR = params.genome ? params.genomes[ params.genome ].salmon_index_3UTR  ?: false : false
params.replace_SalmonIndex_ID = params.genome ? params.genomes[ params.genome ].replace_SalmonIndex_ID  ?: false : false

params.apa3utr =  params.genome ? params.genomes[ params.genome ].apa3utr ?: false : false


params.threads = 12
threads = params.threads


if ( params.fasta &&  params.utr_gtf ){
    fasta = file(params.fasta)
	fasta_fai = file("${params.fasta}.fai")
    utr_gtf = file(params.utr_gtf)
	
	if(! fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
	if(! utr_gtf.exists()) exit 1, "utr_gtf file not found: ${params.utr_gtf}"
}


params.layout = 'paired'
def ifPaired = true
if (params.layout =~ /single/ ){
	ifPaired = false
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
			if(! ifPaired ){
				println("fastp single-end!")
				"""
				if [ "${reads}" == "${id}.fastp.fq.gz" ];
				then
					echo "input raw reads name and output clean name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads}  -o ${id}.fastp.fq.gz  -q 20 
				"""
			}else
			{
				println("fastp paired!")
				"""
				if [ "${reads[0]}" == "${id}.R1.fastp.fq.gz" ];
				then
					echo "input raw reads name and output clean reads name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 
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
		errorStrategy 'ignore'
		maxForks 1
		
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
			fasterq-dump --split-3 -o ${id}  $sra
			fastp -i ${id}_1.fastq  -o ${id}.fastp.fq.gz -O -q 20 --thread  2
			"""
		}else{
			"""
			fasterq-dump --split-3 -o ${id}  $sra
			fastp -i ${id}_1.fastq -I ${id}_2.fastq -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  2
			"""
		}
		
	}
	
	//read_ch=fastp_out.map{it ->a=[it[0],it[1 .. 2]] }
    }else{
	exit 1, "you should input --reads or --sra or --bam"
}

read_ch.into{ read_salmon ;read_print}



if(  params.apa3utr && (! params.salmon_index_3UTR ) ){
	apa3utr_bed = file( params.apa3utr)
	
	process utr_to_fasta{
		
		
		publishDir "${params.outdir}/ref_transcripts/", mode: 'link'
		
		input:
		file "output_utrs.bed" from apa3utr_bed
		file fasta
		
		output:
		file 'output_utrs.shortedID.fa' into cal_expr_fa
		file "fa_id_relate.txt" into utr_id_file
		
		"""
		set +u; source activate apa; set -u
		qapa fasta -f ${fasta} output_utrs.bed output_utrs.fa
		python3 ${baseDir}/bin/apa_3utr.py shorten_id output_utrs.fa output_utrs.shortedID.fa
		"""
	}

	process salmon_index_3UTR{
		
		publishDir "${params.outdir}/ref_salmon_index_3UTR/", mode: 'link'
		
		input:
		file 'transcripts.fa' from cal_expr_fa
		
		output:
		file "salmon_index_3UTR/*" into salmon_index_ch
		
		"""
		salmon index -t transcripts.fa -i salmon_index_3UTR -p 4
		"""
	}

}else if( params.salmon_index_3UTR && params.replace_SalmonIndex_ID  ){
	
	salmon_index_ch = Channel
			.fromPath("${params.salmon_index_3UTR }/*")
			.ifEmpty { exit 1, "params.salmon_index_3UTR  not found: ${params.salmon_index_3UTR }" }
	
	if( params.replace_SalmonIndex_ID ){
		utr_id_file = file(params.replace_SalmonIndex_ID)
	}

			
}else{
	exit 1, "must provide --salmon_index_3UTR ,replace_SalmonIndex_ID or --apa3utr "
}


process salmon_APA {
	
	tag {id } 
	input:
	set id, file(read) from read_salmon
	file "transcripts_index/*" from salmon_index_ch.collect()
	
	output:
	file "${id}"  into salmon_slueth, salmon_out
	
	script:
	if(params.single){
					
		"""
		salmon quant -i transcripts_index -l A -r $read -p 8 -o ${id}
		"""
	}else{
					
		"""
		salmon quant -i transcripts_index -l A -1 ${read[0] } -2 ${read[1] } -p 8 -o ${id}
		"""
	}
	
}


process qapa {
	
	publishDir "${params.outdir}/apa_3utr/", mode: 'link'
	
	input:
	file utr_gtf
	file "salmon/*" from salmon_out.collect()
	file utr_id_file
	
	output:
	set "pau_results.txt" ,"salmon" into qapa_out
	
	
	"""
	set +u; source activate apa; set -u
	python3 ${baseDir}/bin/gtf.py to_gencode2biomark5 ${utr_gtf}  db_identifiers
	python3 ${baseDir}/bin/apa_3utr.py replace_SalmonIndex_ID ${utr_id_file} salmon
	qapa quant --db db_identifiers salmon/*/quant_replace.sf > pau_results.txt
	"""
}

