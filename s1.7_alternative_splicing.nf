#!/usr/bin/env nextflow
params.bam  = ""
params.sra  = ""
params.reads = ""
params.update_GTF = 0 
params.genome = ""
params.cleaned  = false
params.saveIntermediateFiles = false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false

params.threads = 16

scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2; scripts_3; scripts_4}

ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";

def text_red = {  str ->  ANSI_RED + str + ANSI_RESET }
def text_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def text_green = {  str ->  ANSI_GREEN + str + ANSI_RESET }
def text_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def text_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def text_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def text_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def text_white = {  str -> ANSI_WHITE + str + ANSI_RESET }


if ( params.fasta ){
    fasta = file(params.fasta)
	fasta_fai = file("${params.fasta}.fai")
	fasta_dict_path = params.fasta.replaceAll(~/\.fa.*?$/, '.dict')
	fasta_dict = file(fasta_dict_path)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}


if ( params.gtf && params.update_GTF == 0 ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}else if(params.update_GTF == 1){
	gtf = file("../s1.1_lncRNA/results/annotations_and_fasta/protein_coding_and_all_lncRNA.gtf")
	if( !gtf.exists() ) exit 1, "file: ../s1.1_lncRNA/results/annotations_and_fasta/protein_coding_and_all_lncRNA.gtf does not found\
	\nPlease check step s1.1_lncRNA has complete?"
}else{
	exit 1, "GTF file not found!"
}


if( params.star_index ){
    star_index = Channel
        .fromPath("${params.star_index}/*")
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    star_index.into{star_index; star_index_1; star_index_2}
}


params.layout = 'paired'
def ifPaired = true
if (params.layout =~ /single/ ){
	ifPaired = false
}

if (params.bam){
	Channel
	.fromPath(params.bam)
	.ifEmpty{error "Can not find any bamfiles matching: ${params.bam}"}
	.map{file ->
		def key = file.name.toString().tokenize('.').get(0)
		return tuple(key,file)
		}
	.groupTuple()
	.set{ch}

	bam_infile = ch

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


if (params.reads || params.sra ){
	read_ch.into{ reads_star; read_star_2pass;  reads_print }

	process star_1_pass {
		tag {id}
		input:
		set id, file(reads) from reads_star
		file "STARIndex/*" from  star_index.collect()
		
		output:
		
		file "${id}.SJ.out.tab" into star_1_out
			
		"""
		set +u; source activate pipeOne_py3; set -u
		STAR 	--runThreadN ${params.threads } \\
				--genomeDir STARIndex \\
				--outSAMtype None \\
				--readFilesIn ${reads} \\
				--outFileNamePrefix	${id}. \\
				--readFilesCommand zcat \\
				--outFilterMultimapNmax 50 \\
				--outFilterMultimapScoreRange 3 \\
				--outFilterScoreMinOverLread 0.7 \\
				--outFilterMatchNminOverLread 0.7 \\
				--outFilterMismatchNmax 10 \\
				--alignIntronMax 500000 \\
				--alignMatesGapMax 1000000 \\
				--sjdbScore 2
		
		"""
		
	}

	// Aligned.sortedByCoord.out.bam
	// Aligned.out.bam
	// Aligned.out.sam


	process star_genomeGenerate_for_2pass{
		
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
		STAR --runMode genomeGenerate --genomeDir \$genomeDir \\
			--genomeFastaFiles ${fasta} \\
			--sjdbGTFfile ${gtf} \\
			 --sjdbFileChrStartEnd SJ.out.tab  --runThreadN  ${params.threads}
		
		"""
	}


	process star_2_pass{
		tag {id}
		
		errorStrategy 'ignore'

		publishDir "${params.outdir}/star2pass", mode: 'link'
		
		input:
		file "STARIndex/*" from  STARIndex_2pass.collect()
		set  id, file(reads) from read_star_2pass
		file gtf
		file fasta
		
		output:
		set id, "${id}.bam", "${id}.bam.bai" into star_2pass_out
		set "${id}.bam", "${id}.bam.bai" into bams_bai
		
		"""
		set +u; source activate pipeOne_py3; set -u
		STAR --runThreadN ${params.threads } \\
			--genomeDir STARIndex \\
			--readFilesIn ${reads} \\
			--outFileNamePrefix ${id}. \\
			--readFilesCommand zcat \\
			--outSAMtype BAM SortedByCoordinate \\
			--outSAMstrandField intronMotif \\
			--outSAMattributes NH HI NM MD AS XS \\
			--outSAMheaderHD @HD VN:1.4 \\
			--outFilterMultimapNmax 50 \\
			--outFilterMultimapScoreRange 3 \\
			--outFilterScoreMinOverLread 0.7 \\
			--outFilterMatchNminOverLread 0.7 \\
			--outFilterMismatchNmax 10 \\
			--alignIntronMax 500000 \\
			--alignMatesGapMax 1000000 \\
			--sjdbScore 2 --limitBAMsortRAM 1443648494946
			
		mv ${id}.Aligned.sortedByCoord.out.bam ${id}.bam
		samtools index ${id}.bam
		
		"""
	}


	bam_files = star_2pass_out


}else if(params.bam ){
	Channel
		.fromPath(params.bam)
		.ifEmpty{error "Can not find any reads matching: ${params.read}"}
		.map{file ->
			def key = file.name.toString().tokenize('.').get(0)
			return tuple(key,file)
			}
		.groupTuple()
		.set{bam_infile}
			
			
	process bamIndex{
		tag {id}
		input:
		set id, file(bam) from bam_infile
		
		output:
		set id, "${bam}", "${bam}.bai" into bam_files
		set "${id}.bam", "${id}.bam.bai" into bams_bai
		
		"""
		set +u; source activate pipeOne_py3; set -u
		samtools index ${bam}
		"""
	}
	
}


bam_files.into{ bam_files; bam_files_f; bam_files_quant_first;}
bam_files_first =  bam_files_f.first()

bams_bai.into{ bams_bai;  bam_bai_quant }


test_bam = file("${baseDir}/dat/test.bam")
test_bai = file("${baseDir}/dat/test.bam.bai")
process Single_graphs_for_generate_gtf_pickle {
	
	input:
	file gtf
	file test_bam
	file test_bai

	output:
	file "${gtf}.pickle" into gtf_pickle
	
	"""
	set +u; source activate pipeOne_AS; set -u
	spladder build -o spladder_out \\
                   -a ${gtf} \\
                   -b ${test_bam} \\
                   --merge-strat single \\
                   --no-extract-ase \\
				   --parallel 8 \\
				   -n 100 \\
				   -c 3
	"""
}
	
//def a =graph_one.toList()
//a.println()

// bam_files
	// .filter{ it[0].toString() !="SRR7062359"  }
	// .println{ text_blue(it[0]) }

// bam_files.println()
// println	 text_blue('hello ') + graph_one.collect()

process Single_graphs {
	label 'moreParralle'
	tag {id}
	
	input:
	set id, file(bam), file(bai) from bam_files
	file gtf
	file "${gtf}.pickle" from gtf_pickle
	
	output:
	file "spladder_out/spladder/genes_graph_conf3.${id}.pickle" into graphs_all_merge, graphs_all_quant
	set id, "spladder_out/spladder/genes_graph_conf3.${id}.pickle" into graphs_ch 
	"""
	set +u; source activate pipeOne_AS; set -u
	spladder build -o spladder_out \\
                   -a ${gtf} \\
                   -b ${bam} \\
                   --merge-strat single \\
                   --no-extract-ase \\
				   --parallel ${params.threads } \\
				   -n 100 \\
				   -c 3
	"""
}



process Merged_graph {
	
	input:
	file "spladder_out/spladder/*" from graphs_all_merge.collect()
	file "*"  from bams_bai.collect()
	file gtf
	file "${gtf}.pickle" from gtf_pickle
	
	output:
	file "spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle" into graph_merge
	//file "*.hdf5" into mergeGraph_hdf5_ch, hdf5_all
	
	"""
	set +u; source activate pipeOne_AS; set -u
	ls |grep "bam\$" > alignments.txt
	mkdir -p  spladder_out/tmp
	
	spladder build -o spladder_out \\
               -a ${gtf} \\
               -b alignments.txt \\
               --merge-strat merge_graphs \\
               --no-extract-ase  \\
			   --parallel 48 \\
			   -n 100 \\
				-c 3
	"""
}


// transform into [id, [bam, bai] ]

bam_files_quant_first
	.map{ it -> a = [it[0], it[1 .. 2]] }
	.set{ bam_files_quant_first }

graphs_ch.into{ graphs_ch; graphs_ch2 }

/*
mergeGraph_hdf5_ch
	.flatten()
	.map{file ->
		def key = file.name.toString().tokenize('.').get(0)
		return tuple(key, file)
		}
	.groupTuple()
	.combine(bam_files_quant_first, by: 0)
	.combine(graphs_ch, by: 0)
	.set{ quant_ch }
*/

bam_files_quant_first 
	.combine(graphs_ch, by: 0)
	.set{ quant_ch }
	
quant_ch.into{ quant_ch; quant_ch_merge; m_print }
//m_print.println()



process Quant_each {
	
	label 'moreParralle'
	input:
	
	file "spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle" from graph_merge
	//set id, file(bam), file(hdf5),"spladder_out/spladder/genes_graph_conf3.${id}.pickle" from quant_ch
	set id, file(bam), "spladder_out/spladder/genes_graph_conf3.${id}.pickle" from quant_ch
	
	file gtf
	file gtf_pickle
	
	output:
	file "spladder_out/spladder/genes_graph_conf3.merge_graphs.${id}.count.hdf5" into count_all
	
	
	"""
	set +u; source activate pipeOne_AS; set -u
	spladder build -o spladder_out -a ${gtf} -b ${id}.bam \\
                   --merge-strat merge_graphs \\
                   --no-extract-ase \\
                   --quantify-graph \\
                   --qmode single  --parallel 10 \\
				   -n 100 \\
				   -c 3
	"""
}




process QuantAll_and_EventCalling{
	
	publishDir "${params.outdir}/spladder_out/", mode: 'link'
	
	input:
	file "*" from bam_bai_quant.collect()
	//file "*" from hdf5_all.collect()
	file "spladder_out/spladder/*" from graphs_all_quant.collect()
	file "spladder_out/spladder/*" from count_all.collect()
	file "spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle" from graph_merge
	file gtf
	file gtf_pickle
	
	output:
	file "spladder_out/merge_graphs*" into spladder_out
	
	"""
	set +u; source activate pipeOne_AS; set -u
	ls |grep "bam\$" > alignments.txt
	mkdir -p  spladder_out/tmp
	
	echo "merge counts ...."
	spladder build -o spladder_out \\
               -a ${gtf} \\
               -b alignments.txt \\
               --merge-strat merge_graphs \\
               --no-extract-ase \\
               --quantify-graph \\
               --qmode collect --parallel 48 \\
			   -n 100 \\
			   -c 3
			   
	echo "Event Calling ...."		   
	spladder build -o spladder_out \\
		-a ${gtf} \\
		-b alignments.txt \\
		--parallel ${params.threads } \\
		--event-types  exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip \\
		--output-txt 
			
	"""
}


process PSI_table{

	publishDir "${params.outdir}/spladder_out_table/", mode: 'link'

	input:
	file "spladder_out/*" from spladder_out.collect()
	
	output:
	file "*psi.txt.gz" into spladder_PSI_out
	
	shell:
	''' 
	set +u; source activate pipeOne_py3; set -u
	mkdir -p  spladder_PSI/
	for i in `ls spladder_out/*.confirmed.txt.gz`
	do
		python3 !{baseDir}/bin/alternative_splicing.py spladder_PSI_tab ${i}
	done
	'''
}