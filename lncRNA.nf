#!/usr/bin/env nextflow

params.genome = "hg19"
params.species = 'human'
params.genecode_gtf = params.genome ? params.genomes[ params.genome ].genecode_gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.lncpedia_gtf  = params.genome ? params.genomes[ params.genome ].lncpedia_gtf  ?: false : false
params.genecode_lncRNA_gtf  = params.genome ? params.genomes[ params.genome ].genecode_lncRNA_gtf   ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2_index ?:false :false
params.clusterOptions=""
params.bams = ""
params.reads  = "/public1/pub/guohh/other/human_for_xiong/rawdata/*_{1,2}.fq.gz"
params.single = false
params.cleaned  = false

genecode_lncRNA_gtf = file(params.genecode_lncRNA_gtf)

if ( params.fasta ){
    fasta = file(params.fasta)
	fasta_fai = file("${params.fasta}.fai")
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}else{exit 1, "No reference genome file specified!"}

if ( params.genecode_gtf ){
    gtf = file(params.genecode_gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.genecode_gtf}"
}else{exit 1, "No reference gtf file specified!"}

if ( params.lncpedia_gtf ){
    lncpedia_gtf = file(params.lncpedia_gtf)
    if( !lncpedia_gtf.exists() ) exit 1, "lncpedia_gtf file not found: ${params.lncpedia_gtf}"
}else{exit 1, "No lncpedia gtf file specified!"}


scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2}
params.threads = 12
threads = params.threads


params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded =false

if (params.bams){
	Channel
	.fromPath(params.bams)
	.ifEmpty{error "Can not find any reads matching: ${params.reads}"}
	.map{file ->
		def key = file.name.toString().tokenize('.').get(0)
		return tuple(key,file)
		}
	.groupTuple()
	.set{ch}

	bams_ch = ch

}

if(params.reads ){
	
	// defined single or paired reads !
	if(!params.single){
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

	if (! params.cleaned ){
	
		process fastp{
			tag {id}
			
			publishDir "${params.outdir}/fastp/", mode: 'copy',
				saveAs: {filename -> 
					if(filename =~ /fastq.gz/) "clean/${filename}"
					else "report/${id}.${filename}" 
					}
			
			input:
			set id, file(reads) from ch
			
			output:
			set id, "${id}.R1.fastq.gz","${id}.R2.fastq.gz" into fastp_out
			set  "fastp.html", "fastp.json"
			
			"""
			fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastq.gz -O ${id}.R2.fastq.gz -q 20 --thread  4
			"""
		
		}
		
		read_ch=fastp_out.map{it ->a=[it[0],it[1 .. 2]] }

		}else{
			read_ch =ch
		}
		
	read_ch.into{ hisat2_ch; reads_kallisto }
		
	
	if ( params.hisat2_index && !params.bams ){
		
		hisat2_indices = Channel
			.fromPath("${params.hisat2_index}*")
			.ifEmpty { exit 1, "hisat2_index not found: ${params.hisat2_index}" }
			hisat2_base= params.hisat2_index.split('/')[-1]
	
		process hisat2 {
			tag {id}
			publishDir "${params.outdir}/hisat2/", mode: 'copy',
				saveAs: {filename -> 
					if(filename =~ /bam/) filename
					else if (filename =~/log/) "logs/${filename}"
					}
			
			input:
			set id, file(reads)   from hisat2_ch
			file "hisat2_index/*" from hisat2_indices.collect()
			
			output:
			set id, "${id}.hisat2.sortbycoordinate.bam" into hisat2_out
			file "${id}.hisat2.log" into hisat2_log
			
			
			
			script:
			def rnastrandness = ''
			if (params.forward_stranded && !params.unstranded){
				rnastrandness = params.single ? '--rna-strandness F' : '--rna-strandness FR'
			} else if (params.reverse_stranded && !params.unstranded){
				rnastrandness = params.single ? '--rna-strandness R' : '--rna-strandness RF'
			}
			
			if(! params.single){
				"""
				zcat ${reads[0]} > read1.fq
				zcat ${reads[1]} > read2.fq
				hisat2 -p $threads --dta  $rnastrandness  -x hisat2_index/$hisat2_base \
				-1 read1.fq  -2 read2.fq   2> ${id}.hisat2.log | samtools sort -@ 8 - -o ${id}.hisat2.sortbycoordinate.bam 
				samtools index ${id}.hisat2.sortbycoordinate.bam
				rm read*fq
				"""
			}else{
				
				"""
				zcat ${reads[0]} > read1.fq
				hisat2 -p $threads --dta  $rnastrandness -x  hisat2_index/$hisat2_base \
				-U  read1.fq   2> ${id}.hisat2.log  | samtools sort -@ 8 - -o ${id}.hisat2.sortbycoordinate.bam
				samtools index ${id}.hisat2.sortbycoordinate.bam 
				rm read*fq 
				"""
			}
	 
		}

		bams_ch = hisat2_out
	
	}
	
}else if( !params.bams ){
	exit 1, "you should input --bams or --reads "

}

bams_ch.into{ bam_stringite; bam_stringtieFPKM; bam_featurecounts;}

process stringtie {
    tag { id }
    publishDir "${params.outdir}/stringtie", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
            else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
            else "$filename"
        }

    input:
    set id, file(bam) from bam_stringite
    file gtf
	
    output:
    file "${id}_transcripts.gtf" into stringtie_gtf
   

    script:
    def st_direction = ''
    if (params.forward_stranded && !params.unstranded){
        st_direction = "--fr"
    } else if (params.reverse_stranded && !params.unstranded){
        st_direction = "--rf"
    }
	
    """
    stringtie $bam $st_direction \\
			-o ${id}_transcripts.gtf \\
			-p $threads -G $gtf
    """
}


process taco{
	publishDir "${params.outdir}/taco/", mode: 'copy'
	input:
	file gtf
	file "gtf/*" from stringtie_gtf.collect()
	
	output:
	file "taco_out/assembly.gtf" into taco_out
	
	
	"""
	ls gtf | awk '{print "gtf/"\$0}' > gtf_list
	taco_run -p ${threads} -o taco_out gtf_list
	"""
}


process gffcompare{
	publishDir "${params.outdir}/gffcompare/", mode: 'copy'
	input:
	file fasta
	file fasta_fai
	file gtf
	file lncpedia_gtf
	file "assembly.gtf" from taco_out.collect()
	
	
	output:
	set  "novel_lncRNA_candidate.gtf","novel_lncRNA_candidate.fa" into novel_lncRNA_candidate
	file "novel_lncRNA_candidate.fa" into novel_lncRNA_candidate_fa
	file "*"
	file "gffcmp.assembly.gtf.tmap" into gffcompare_tmap
	
	"""
	#!/bin/bash
	cat ${gtf} ${lncpedia_gtf} >all.gtf
	gffcompare -r all.gtf assembly.gtf 
	awk '\$3 =="x"||\$3=="u"||\$3=="i"{print}' gffcmp.assembly.gtf.tmap > novel.gtf.tmap    
    awk '\$10 >=200{print}' novel.gtf.tmap > novel.lncRNA.gtf.tmap
	cut -f5 novel.lncRNA.gtf.tmap  >novel_lncRNA_candidate.list
	gtf.py get_by_trans_id assembly.gtf novel_lncRNA_candidate.list  novel_lncRNA_candidate.gtf
	gffread novel_lncRNA_candidate.gtf -g ${fasta} -w novel_lncRNA_candidate.fa -W
	"""
}

novel_lncRNA_candidate_fa.into{lncRNA_candidate_fa_1; lncRNA_candidate_fa_2; }
cpatpath = Channel.fromPath("${params.cpatpath}/dat/*")

process cpat {
	
	publishDir "${params.outdir}/CPAT/", mode: 'copy'
	
	input:
	file "novel_lncRNA_candidate.fa" from lncRNA_candidate_fa_1.collect()
	file  "cpatpath/dat/*" from cpatpath.collect()
	
	output:
	file "CPAT.out" into cpat_out
	
	script:
	if(params.species=="human"){
        '''
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x cpatpath/dat/Human_Hexamer.tsv \
                                       -d cpatpath/dat/Human_logitModel.RData \
                                       -o CPAT.out
        '''
    }else if (params.species=="mouse"){
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x cpatpath/dat/Mouse_Hexamer.tsv \
                                       -d cpatpath/dat/Mouse_logitModel.RData \
                                       -o CPAT.out
        '''

    }else if (params.species=="zebrafish"){
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x cpatpath/dat/zebrafish_Hexamer.tsv \
                                       -d cpatpath/dat/zebrafish_logitModel.RData \
                                       -o CPAT.out
        '''
    }else {
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x cpatpath/dat/fly_Hexamer.tsv \
                                       -d cpatpath/dat/fly_logitModel.RData \
                                       -o CPAT.out
        '''
    }

}

process PLEK {
    validExitStatus 0,1,2
	publishDir "${params.outdir}/PLEK/", mode: 'copy'
	
	input:
	file "novel_lncRNA_candidate.fa" from lncRNA_candidate_fa_2.collect()
	
	output:
	file "PLEK.out" into plekout
    
	"""
    #!/bin/bash
	PLEK.py -fasta novel_lncRNA_candidate.fa -out PLEK -thread $threads
	ln -s PLEK PLEK.out
	"""
}


process prepare_reference_gtf{
    publishDir "${params.outdir}/reference_gtf_info/", mode: 'copy'
    label 'local'
    
	input:
	file "gencode_annotation.gtf" from gtf
	file "lncpedia.gtf" from lncpedia_gtf
	file genecode_lncRNA_gtf
    
    output:
    set "gencode_protein_coding.gtf", "known_lncRNA.gtf","known_lncRNA.list" into known_gtf_list_ch
    file "known_lncRNA.gtf" into known_lncRNA_ch
    set "known_lncRNA.list", "nonOverlap_lncpedia.list","gencode_protein_coding_tran_id.list"
    
    """
    #!/bin/bash
	gffcompare -r ${genecode_lncRNA_gtf} lncpedia.gtf 
	awk '\$3 =="x"||\$3=="u"||\$3=="i"{print}' gffcmp.lncpedia.gtf.tmap > nonOverlap.gtf.tmap    
	cut -f5 nonOverlap.gtf.tmap    > nonOverlap_lncpedia.list
	gtf.py get_by_trans_id lncpedia.gtf nonOverlap_lncpedia.list  nonOverlap_lncpedia.gtf
    cat gencode_annotation.gtf  |grep "protein_coding" > gencode_protein_coding.gtf
    cat ${genecode_lncRNA_gtf} nonOverlap_lncpedia.gtf > known_lncRNA.gtf
    gtf.py get_id gencode_protein_coding.gtf transcript_id gencode_protein_coding_tran_id.list
    gtf.py get_id  known_lncRNA.gtf transcript_id  known_lncRNA.list
    """
	
}

known_gtf_list_ch.into{ known_gtf_list_ch; known_gtf_list_ch_1}

process filter_coding_potentail {
    publishDir "${params.outdir}/filter_coding_potentail/", mode: 'copy',
		saveAs: { filename-> 
			if( filename =~ /(list|fa|gtf)$/)filename
			}
	publishDir "${params.outdir}/novel_lncRNA/", mode: 'copy',
		saveAs: { filename-> 
			if( filename =~ /novel_lncRNA/)filename
			}
	
	publishDir "${params.outdir}/novel_lncRNA/tmp/", mode: 'copy'
			
    scratch false
	
    
    input:
    file gtf
	file fasta
	file fasta_fai
    set "gencode_protein_coding.gtf", "known_lncRNA.gtf","known_lncRNA.list" from known_gtf_list_ch.collect()
	set "novel_lncRNA_candidate.gtf","novel_lncRNA_candidate.fa" from novel_lncRNA_candidate.collect()
    file "CPAT.out" from cpat_out.collect()
    file "PLEK.out" from plekout.collect()
    
    
    output:
    file "protein_coding_and_all_lncRNA.gtf" into cal_expr_gtf
    file "protein_coding_and_all_lncRNA.fa"  into cal_expr_fa
    set "all_lncRNA.gtf","all_lncRNA.list","novel_lncRNA.fa","novel_lncRNA.gtf", "novel_lncRNA.list" into novel_list
	file "TUCP*"
	file "all_lncRNA.list" into cal_deg_ch
	file "*"
	
    """
	echo "filter coding potential!"
    awk '\$1 == "Coding"{print \$3}' PLEK.out |cut -f2 -d '>' >PLEK.coding
	cpat_cutoff.pl CPAT.out ${params.species}
	cat  CPAT.coding  PLEK.coding |sort -u > all_coding_potential.list
	cat  novel_lncRNA_candidate.fa |grep '>'|cut -d ' ' -f1 |cut -d '>' -f2 >all.list
	nong_get_venn_results.pl  -a all.list -i all_coding_potential.list -o protein_coding.tmp -r  novel_lncRNA.list.allExons
	#all novel lncRNA
	get_seq_by_id.py novel_lncRNA_candidate.fa novel_lncRNA.list.allExons novel_lncRNA.fa.allExons
	gtf.py get_by_trans_id novel_lncRNA_candidate.gtf novel_lncRNA.list.allExons  novel_lncRNA.gtf.allExons
	
	#TUCP
	get_seq_by_id.py novel_lncRNA_candidate.fa all_coding_potential.list TUCP.fa
	gtf.py get_by_trans_id novel_lncRNA_candidate.gtf all_coding_potential.list  TUCP.gtf
	
	##get multi-exons lncRNAs
    cat novel_lncRNA.gtf.allExons known_lncRNA.gtf > all_lncRNA.tmp.gtf
	gtf.py to_info all_lncRNA.tmp.gtf all_lncRNA.tmp.info
	awk '\$8>1' all_lncRNA.tmp.info |cut -f1 > all_lncRNA.multi.exons.list
	gtf.py get_by_trans_id all_lncRNA.tmp.gtf all_lncRNA.multi.exons.list  all_lncRNA.gtf
	##get multi-exons novel_lncRNA
	gtf.py get_by_trans_id novel_lncRNA_candidate.gtf all_lncRNA.multi.exons.list  novel_lncRNA.gtf
	gtf.py get_id  novel_lncRNA.gtf transcript_id  novel_lncRNA.list 
	gffread novel_lncRNA.gtf -g ${fasta} -w novel_lncRNA.fa -W
	
	#for expression calculate
	cat novel_lncRNA.list known_lncRNA.list >all_lncRNA.list
    cat gencode_protein_coding.gtf all_lncRNA.gtf > protein_coding_and_all_lncRNA.gtf
	gffread all_lncRNA.gtf -g ${fasta} -w all_lncRNA.fa -W
	gffread protein_coding_and_all_lncRNA.gtf -g ${fasta} -w protein_coding_and_all_lncRNA.fa -W
    """
}
cal_expr_gtf.into{cal_expr_gtf; cal_expr_gtf_1; cal_expr_gtf_2;cal_expr_gtf_3;}

process featureCounts{
    tag { id }
    publishDir "${params.outdir}/featureCounts/samples", mode: 'copy'
	input:
	set id, file(bam) from bam_featurecounts
	file "ref_gtf" from cal_expr_gtf.collect()
	
	output:
    file "${id}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
  
    script:
    def featureCounts_direction = 0
    if (params.forward_stranded && !params.unstranded) {
        featureCounts_direction = 1
    } else if (params.reverse_stranded && !params.unstranded){
        featureCounts_direction = 2
    }
    """
    featureCounts -a ref_gtf -g gene_id -o ${id}_gene.featureCounts.txt -p -s $featureCounts_direction $bam
   
    """
}

process merge_featureCounts {
    
    publishDir "${params.outdir}/featureCounts/", mode: 'copy'

    input:
    file input_files from featureCounts_to_merge.collect()
    
    
    output:
    file 'merged_gene_counts.txt' into featureCounts_res

    script:
    """
    merge_featurecounts.py -o merged_gene_counts.txt -i $input_files
	sed -i -e '1s/ENSEMBL_ID\t//' -e 's/_gene.featureCounts.txt//g' merged_gene_counts.txt
    """
}

process format_lncRNA_info {
	
	publishDir "${params.outdir}/novel_lncRNA/", mode: 'copy'
	label 'local'
	
	input:
	file "scripts/*" from scripts_1.collect()
	set "all_lncRNA.gtf","all_lncRNA.list","novel_lncRNA.fa","novel_lncRNA.gtf", "novel_lncRNA.list" from novel_list.collect()
	set "gencode_protein_coding.gtf", "known_lncRNA.gtf","known_lncRNA.list" from known_gtf_list_ch_1.collect()
	
	output:
	set "all_lncRNA_info.tsv", "novel_lncRNA_info.tsv" into lncRNA_res
	
	'''
    gffcompare -r gencode_protein_coding.gtf all_lncRNA.gtf
	#convert to bed and closest to bed
	awk '$3=="exon"{print}'  gencode_protein_coding.gtf >gencode_protein_coding.exon.gtf
	awk '$3=="exon"{print}'  all_lncRNA.gtf             >all_lncRNA.exon.gtf
    
	gtf2bed < gencode_protein_coding.exon.gtf |cut -f1-6 >gencode_protein_coding.exon.gtf.bed
	gtf2bed < all_lncRNA.exon.gtf |cut -f1-6 >all_lncRNA.exon.gtf.bed
    sort -k1,1 -k2,2n gencode_protein_coding.exon.gtf.bed > gencode_protein_coding.exon.gtf.bed.sorted
    sort -k1,1 -k2,2n all_lncRNA.exon.gtf.bed > all_lncRNA.exon.gtf.bed.sorted
	bedtools closest -a all_lncRNA.exon.gtf.bed.sorted -b gencode_protein_coding.exon.gtf.bed.sorted  -d |sort -u > cloest.exon.txt
	
    # convert to gene bed and cloest to bed
    python3 scripts/gtf2.py to-gene-bed gencode_protein_coding.exon.gtf gencode_protein_coding.exon.gtf.genebed
    python3 scripts/gtf2.py to-gene-bed all_lncRNA.exon.gtf all_lncRNA.exon.gtf.genebed
    sort -k1,1 -k2,2n gencode_protein_coding.exon.gtf.genebed > gencode_protein_coding.exon.gtf.genebed.sorted
    sort -k1,1 -k2,2n all_lncRNA.exon.gtf.genebed > all_lncRNA.exon.gtf.genebed.sorted
    bedtools closest -a all_lncRNA.exon.gtf.genebed.sorted -b gencode_protein_coding.exon.gtf.genebed.sorted  -d |sort -u > cloest.gene.txt
	
    #gtf to info
	python3 scripts/gtf2.py to-info  all_lncRNA.gtf all_lncRNA.info
	python3 scripts/gtf2.py to-info  gencode_protein_coding.gtf gencode_protein_coding.info.tsv
	
	#add info
	python3 scripts/gene_info_add.py add-gffcompre  all_lncRNA.info gffcmp.all_lncRNA.gtf.tmap all_lncRNA.info.gffcmp
    python3 scripts/gene_info_add.py  add-2-bedtool-closest  all_lncRNA.info.gffcmp cloest.gene.txt cloest.exon.txt all_lncRNA_info.tsv.tmp1 >other.txt
    python3 scripts/gene_info_add.py add-geneName all_lncRNA_info.tsv.tmp1 gencode_protein_coding.info.tsv all_lncRNA_info.tmp3
	mv all_lncRNA_info.tmp3 all_lncRNA_info.tsv
	cat all_lncRNA_info.tsv |awk 'NR==1 || $1 ~/^TU/' >novel_lncRNA_info.tsv
	'''
}


params.datainfo = ""
cal_deg_ch.into{cal_deg_ch;cal_deg_ch_1;cal_deg_ch_2}
if(params.datainfo){
	datainfo = file(params.datainfo)
	process DEG_by_DESeq2_with_featureCounts {
		publishDir "${params.outdir}/deg_DESeq2/featureCounts/", mode: 'copy'
		label 'local'
		input:
		file 'merged_gene_counts.txt' from featureCounts_res.collect()
		file "all_lncRNA.list" from cal_deg_ch.collect()
		file "ref_gtf" from cal_expr_gtf_1.collect()
		file "datainfo.txt" from datainfo
		
		output:
		file '*tsv'
		file '*tiff'
		file "files_cache__"
		
		script:
		"""
		gtf.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
		deseq2.R merged_gene_counts.txt
		
		mkdir -p files_cache__
		ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
		"""
	}

}


if (params.reads){
	
	process kallisto_index{
		
		input:
		file "protein_coding_and_all_lncRNA.fa" from cal_expr_fa
		
		output:
		file "kallisto_index" into kallisto_index_ch
		
		"""
		kallisto index --index=kallisto_index protein_coding_and_all_lncRNA.fa
		"""
	}
	
	process kallisto_abundance_estimates{
		publishDir "${params.outdir}/kallisto/", mode: 'copy'
		
		tag {id } 
		input:
		set id, file(reads) from reads_kallisto
		file "kallisto/*" from kallisto_index_ch.collect()
		
		output:
		file "${id}"  into kallisto_slueth,kallisto_res
		 
		
		script:
		def st_direction = ''
		if (params.forward_stranded && !params.unstranded) {
			st_direction = '--fr-stranded'
		} else if (params.reverse_stranded && !params.unstranded){
			st_direction = '--rf-stranded'
		}
		"""
		kallisto quant -b 100 $st_direction --index=kallisto/kallisto_index --output-dir=${id} --threads=${threads}  ${reads}
		
		"""
	}
	
	process kallisto_merge {
		label 'local'
		publishDir "${params.outdir}/kallisto/", mode: 'copy'
		
		input:
		file "ref_gtf" from cal_expr_gtf_1.collect()
		file "samples/*" from kallisto_res.collect()
		
		output:
		file "kallisto_gene_est_counts.tsv" into kallisto_tab
		file "files_cache__"
		file "*tsv"
		
		"""
		gtf.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
		tximport.R
		mkdir -p files_cache__
		ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
		"""
	}
	
	if(params.datainfo){
		process differential_anlysis_by_slueth{
			publishDir "${params.outdir}/deg_slueth/", mode: 'copy'
			
			input:
			file "samples/*" from kallisto_slueth.collect()
			file "datainfo.txt" from datainfo
			output:
			file "*"
			
			"""
			echo test
			"""
		
		}
		
		process DEG_by_DESeq2_with_kallisto {
			publishDir "${params.outdir}/deg_DESeq2/kallisto/", mode: 'copy'
			label 'local'
			input:
			
			file "all_lncRNA.list" from cal_deg_ch_1.collect()
			file "ref_gtf" from cal_expr_gtf_2.collect()
			file "datainfo.txt" from datainfo
			file "kallisto_gene_est_counts.tsv" from  kallisto_tab.collect()
			
			output:
			file '*tsv'
			file '*tiff'
			file "files_cache__"
			
			script:
			"""
			gtf.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
			deseq2.R kallisto_gene_est_counts.tsv
			
			mkdir -p files_cache__
			ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
			"""
	}
	
	}
	
}
