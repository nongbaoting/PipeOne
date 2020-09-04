#!/usr/bin/env nextflow

params.genome = ""
params.species = params.genome ? params.genomes[ params.genome ].species ?: false : false
params.genecode_gtf = params.genome ? params.genomes[ params.genome ].genecode_gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.lncpedia_gtf  = params.genome ? params.genomes[ params.genome ].lncpedia_gtf  ?: false : false
params.genecode_lncRNA_gtf  = params.genome ? params.genomes[ params.genome ].genecode_lncRNA_gtf   ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2_index ?:false :false

params.bam = ""
params.reads  = ""
params.single = false
params.cleaned  = false
params.featureCounts = false
params.saveIntermediateFiles = false
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


println("genome: " + params.genome + "\n" + "species: "+ params.species )

scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2}
params.threads = 16
threads = params.threads

params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded =false

params.layout = 'paired'
def ifPaired = true
if (params.layout =~ /(?i)single(?-i)/ ){
	ifPaired = false
}

if (params.bam){
	Channel
	.fromPath(params.bam)
	.ifEmpty{error "Can not find any bam files matching: ${params.bam}"}
	.map{file ->
		def key = file.name.toString().tokenize('.').get(0)
		return tuple(key,file)
		}
	.groupTuple()
	.set{ch}

	bam_ch = ch

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


read_ch.into{ hisat2_ch; reads_kallisto;reads_salmon ;reads_print}
//reads_print.println()


if ( params.hisat2_index && !params.bam ){
		
    hisat2_indices = Channel
        .fromPath("${params.hisat2_index}*")
        .ifEmpty { exit 1, "hisat2_index not found: ${params.hisat2_index}" }
        hisat2_base= params.hisat2_index.split('/')[-1]

    process hisat2 {
			
        tag {id}
        publishDir "${params.outdir}/hisat2/", mode: 'link',
            saveAs: {filename -> 
                if(filename =~ /bam/   )  "bam/$filename"
                else if (filename =~/log/) "logs/${filename}"
				else null
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
            rnastrandness = ifPaired ? '--rna-strandness FR' : '--rna-strandness F'
        } else if (params.reverse_stranded && !params.unstranded){
            rnastrandness = ifPaired ? '--rna-strandness RF' : '--rna-strandness R' 
        }
        
        if(ifPaired){
            """
			set +u; source activate lncRNA; set -u
            hisat2 -p $threads --dta  $rnastrandness  -x hisat2_index/$hisat2_base \
            -1 ${reads[0]} -2 ${reads[1]}   2> ${id}.hisat2.log | samtools sort -@ 8 - -o ${id}.hisat2.sortbycoordinate.bam 
            samtools index ${id}.hisat2.sortbycoordinate.bam
			conda deactivate
            """
        }else{
            
            """
			set +u; source activate lncRNA; set -u
            hisat2 -p $threads --dta  $rnastrandness -x  hisat2_index/$hisat2_base \
            -U  ${reads}   2> ${id}.hisat2.log  | samtools sort -@ 8 - -o ${id}.hisat2.sortbycoordinate.bam
            samtools index ${id}.hisat2.sortbycoordinate.bam
			conda deactivate
            """
        }
 
    }

    bam_ch = hisat2_out


}else if( !params.bam ){
    exit 1, "you should input --bam or --reads "

}

bam_ch.into{ bam_stringite; bam_stringtieFPKM; bam_featurecounts;}

process stringtie {
    tag { id }
    publishDir "${params.outdir}/stringtie", mode: 'link',
        saveAs: {filename ->
			if( params.saveIntermediateFiles ){
				if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
				else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
				else "$filename"
			}else null
            
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
	set +u; source activate lncRNA; set -u
    stringtie $bam $st_direction \\
			-o ${id}_transcripts.gtf \\
			-p $threads -G $gtf
    """
}



process taco{
	
	publishDir "${params.outdir}/taco/", mode: 'link', 
    	saveAs: { filename -> if( params.saveIntermediateFiles  && filename =~ /gtf/ ) "$filename" 
							  else null }

	input:
	file gtf
	file "stringtie_gtf/*" from stringtie_gtf.collect()
	
	output:
	file "taco_out/assembly.gtf" into taco_out
	
	
	"""
	set +u; source activate lncRNA; set -u
	ls stringtie_gtf | awk '{print "stringtie_gtf/"\$0}' > gtf_list
	taco_run -p ${threads} -o taco_out gtf_list
	"""
}

// to get known annotate genes/lncRNA and novel lncRNAs
process gffcompare{
	publishDir "${params.outdir}/gffcompare/", mode: 'link',
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
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
	set +u; source activate lncRNA; set -u
	#!/bin/bash
	cat ${gtf} ${lncpedia_gtf} >all.gtf
	gffcompare -r all.gtf assembly.gtf 

	## get multi-exons lncRNAs or single exon length >= 2000	
	python3 ${baseDir}/bin/lncRNA.py gffcompare_tmap gffcmp.assembly.gtf.tmap novel_lncRNA_candidate.list
	
	python3 ${baseDir}/bin/gtf.py get_by_trans_id assembly.gtf novel_lncRNA_candidate.list novel_lncRNA_candidate.gtf
	gffread novel_lncRNA_candidate.gtf -g ${fasta} -w novel_lncRNA_candidate.fa -W
	"""
}

novel_lncRNA_candidate_fa.into{
	lncRNA_candidate_fa_1;
	lncRNA_candidate_fa_2;
	lncRNA_candidate_CPPred; 
	}
	
cpatpath = params.cpatpath

process cpat{
	
	publishDir "${params.outdir}/coding_potential/CPAT/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	file "novel_lncRNA_candidate.fa" from lncRNA_candidate_fa_1.collect()
	
	
	output:
	file "CPAT.out" into cpat_out
	
	shell:
	
	if(params.species =~ /(?i)human(?-i)/){
        '''
		set +u; source activate lncRNA; set -u
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x !{cpatpath}/dat/Human_Hexamer.tsv \
                                       -d !{cpatpath}/dat/Human_logitModel.RData \
                                       -o CPAT.out
        '''
    }else if (params.species =~ /(?i)mouse(?-i)/ ){
		
        '''
		set +u; source activate lncRNA; set -u
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x !{cpatpath}/dat/Mouse_Hexamer.tsv \
                                       -d !{cpatpath}/dat/Mouse_logitModel.RData \
                                       -o CPAT.out
        '''

    }else if (params.species =~/(?i)zebrafish(?-i)/){
		
        '''
		set +u; source activate lncRNA; set -u
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x !{cpatpath}/dat/zebrafish_Hexamer.tsv \
                                       -d !{cpatpath}/dat/zebrafish_logitModel.RData \
                                       -o CPAT.out
        '''
    }else if (params.species =~ /(?i)Fly(?-i)/) {
        '''
		set +u; source activate lncRNA; set -u
        cpat.py -g novel_lncRNA_candidate.fa  \
                                       -x !{cpatpath}/dat/fly_Hexamer.tsv \
                                       -d !{cpatpath}/dat/Fly_logitModel.RData \
                                       -o CPAT.out
        '''
    }

}

process PLEK {
	
    validExitStatus 0,1,2
	publishDir "${params.outdir}/coding_potential/PLEK/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	file "novel_lncRNA_candidate.fa" from lncRNA_candidate_fa_2.collect()
	
	output:
	file "PLEK.out" into plekout
    
	"""
    #!/bin/bash
	set +u; source activate lncRNA; set -u;
	PLEK.py -fasta novel_lncRNA_candidate.fa -out PLEK -thread $threads
	ln -s PLEK PLEK.out
    printf "type\tPLEK_score\ttx_id\n" >PLEK.txt
    more PLEK.out |cut -f1 -d ' ' |sed 's/>//'  >> PLEK.txt
	"""
}


process CPPred {
	publishDir "${params.outdir}/coding_potential/CPPred/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	file "novel_lncRNA_candidate.fa" from lncRNA_candidate_CPPred.collect()
	
	output:
	file "CPPred.out" into CPPred_out
    
	"""
	set +u; source activate lncRNA; set -u
    python2 ${params.CPPred}/bin/CPPred.py -i novel_lncRNA_candidate.fa \
		-hex ${params.CPPred}/Hexamer/Human_Hexamer.tsv -r ${params.CPPred}/Human_Model/Human.range \
		-mol ${params.CPPred}/Human_Model/Human.model -spe Human -o CPPred.out
	"""
}

process prepare_reference_gtf {
    publishDir "${params.outdir}/reference_gtf_info/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles  ? "$filename" : null }
	
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
	set +u; source activate lncRNA; set -u
	gffcompare -r ${genecode_lncRNA_gtf} lncpedia.gtf 
	awk '\$3 =="x"||\$3=="u"||\$3=="i"{print}' gffcmp.lncpedia.gtf.tmap > nonOverlap.gtf.tmap    
	cut -f5 nonOverlap.gtf.tmap > nonOverlap_lncpedia.list
	python3 ${baseDir}/bin/gtf.py get_by_trans_id lncpedia.gtf nonOverlap_lncpedia.list  nonOverlap_lncpedia.gtf
    cat gencode_annotation.gtf  |grep "protein_coding" > gencode_protein_coding.gtf
    cat ${genecode_lncRNA_gtf} nonOverlap_lncpedia.gtf > known_lncRNA.gtf
    python3 ${baseDir}/bin/gtf.py get_id gencode_protein_coding.gtf transcript_id gencode_protein_coding_tran_id.list
    python3 ${baseDir}/bin/gtf.py get_id known_lncRNA.gtf transcript_id known_lncRNA.list
    """
}

known_gtf_list_ch.into{ known_gtf_list_ch; known_gtf_list_ch_1;  known_gtf_list_ch_2;known_gtf_list_ch_3;}

process filter_coding_potentail{

    publishDir "${params.outdir}/annotations_and_fasta/", mode: 'link',
		saveAs: { filename ->  params.saveIntermediateFiles  ? "$filename" : null }

	publishDir "${params.outdir}/novel_lncRNA/", mode: 'link',
		saveAs: { filename-> 
			if( filename =~ /novel_lncRNA/) filename
			}
	
	publishDir "${params.outdir}/novel_lncRNA/tmp/", mode: 'link',
		saveAs: { filename ->  params.saveIntermediateFiles  ? "$filename" : null }

	publishDir "${params.outdir}/coding_potential/", mode: 'link',
		saveAs: { filename ->  filename ==  "coding_potential_sum.tsv" ? "$filename" : null }

			
    input:
    file gtf
	file fasta
	file fasta_fai
    set "gencode_protein_coding.gtf", "known_lncRNA.gtf","known_lncRNA.list" from known_gtf_list_ch.collect()
	set "novel_lncRNA_candidate.gtf", "novel_lncRNA_candidate.fa" from novel_lncRNA_candidate.collect()
    file "CPAT.out" from cpat_out.collect()
    file "PLEK.out" from plekout.collect()
	file "CPPred.out" from CPPred_out.collect()
    
    
    output:
    file "protein_coding_and_all_lncRNA.gtf" into cal_expr_gtf
    file "protein_coding_and_all_lncRNA.fa" into cal_expr_fa
    set "all_lncRNA.gtf","all_lncRNA.list","novel_lncRNA.fa","novel_lncRNA.gtf", "novel_lncRNA.list" into novel_list
	file "all_lncRNA.gtf" into all_lnc_gtf
    file "TUCP*"
	file "all_lncRNA.list" into cal_deg_ch
	file "coding_potential_sum.tsv"
	
    """
	set +u; source activate lncRNA; set -u;
	echo "filter coding potential!"
	python3 ${baseDir}/bin/gtf.py to_info novel_lncRNA_candidate.gtf novel_lncRNA_candidate.info
	python3 ${baseDir}/bin/lncRNA.py combine_coding_prediction novel_lncRNA_candidate.info CPAT.out PLEK.out CPPred.out ${params.species}  novel_lncRNA.list
	
	# all non-coding RNA
	python3 ${baseDir}/bin/get_seq_by_id.py novel_lncRNA_candidate.fa novel_lncRNA.list novel_lncRNA.fa
	python3 ${baseDir}/bin/gtf.py get_by_trans_id novel_lncRNA_candidate.gtf novel_lncRNA.list  novel_lncRNA.gtf
	
	# TUCP
	python3 ${baseDir}/bin/get_seq_by_id.py novel_lncRNA_candidate.fa TUCP.list TUCP.fa
	python3 ${baseDir}/bin/gtf.py get_by_trans_id novel_lncRNA_candidate.gtf TUCP.list  TUCP.gtf
	
	
	# for expression calculate
	cat novel_lncRNA.list known_lncRNA.list >all_lncRNA.list
	cat novel_lncRNA.gtf  known_lncRNA.gtf  >all_lncRNA.gtf
    cat gencode_protein_coding.gtf novel_lncRNA.gtf known_lncRNA.gtf > protein_coding_and_all_lncRNA.tmp.gtf
	python3 ${baseDir}/bin/gtf.py add_gene_name protein_coding_and_all_lncRNA.tmp.gtf protein_coding_and_all_lncRNA.gtf
	gffread all_lncRNA.gtf -g ${fasta} -w all_lncRNA.fa -W
	gffread protein_coding_and_all_lncRNA.gtf -g ${fasta} -w protein_coding_and_all_lncRNA.fa -W
    """
}
cal_expr_gtf.into{cal_expr_gtf; cal_expr_gtf_1; cal_expr_gtf_2;cal_expr_gtf_3;}

if( params.featureCounts){
	process featureCounts{
		
		tag { id }
		publishDir "${params.outdir}/featureCounts/samples", mode: 'link'
		
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
		source activate lncRNA
		featureCounts -a ref_gtf -g gene_id -o ${id}_gene.featureCounts.txt -p -s $featureCounts_direction $bam
	   
		"""
	}

	process merge_featureCounts {
		
		publishDir "${params.outdir}/featureCounts/", mode: 'link'

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

}

process format_lncRNA_info{
	publishDir "${params.outdir}/novel_lncRNA/", mode: 'link'
	errorStrategy 'ignore'
    
	input:
	file "scripts/*" from scripts_1.collect()
	set "all_lncRNA.gtf","all_lncRNA.list","novel_lncRNA.fa","novel_lncRNA.gtf", "novel_lncRNA.list" from novel_list.collect()
	set "gencode_protein_coding.gtf", "known_lncRNA.gtf","known_lncRNA.list" from known_gtf_list_ch_1.collect()
	file "protein_coding_and_all_lncRNA.gtf" from cal_expr_gtf_2.collect()
	
	output:
	set "all_lncRNA_info.tsv", "novel_lncRNA_info.tsv" into lncRNA_res
	file "protein_coding_and_all_lncRNA.info.tsv"
	
	shell:
	'''
	set +u; source activate lncRNA; set -u
    gffcompare -r gencode_protein_coding.gtf all_lncRNA.gtf
	
	#conver to bed and colsest to bed
	awk '$3=="exon"{print}'  gencode_protein_coding.gtf >gencode_protein_coding.exon.gtf
	awk '$3=="exon"{print}'  all_lncRNA.gtf             >all_lncRNA.exon.gtf
	gtf2bed < gencode_protein_coding.exon.gtf |cut -f1-6 >gencode_protein_coding.exon.gtf.bed
	gtf2bed < all_lncRNA.exon.gtf |cut -f1-6 >all_lncRNA.exon.gtf.bed
	bedtools closest -a all_lncRNA.exon.gtf.bed -b gencode_protein_coding.exon.gtf.bed  -d |sort -u > cloest.txt
	
	#gtf to info
	python3 !{baseDir}/bin/gtf.py to-info  all_lncRNA.gtf all_lncRNA.info
	python3 !{baseDir}/bin/gtf.py to-info  gencode_protein_coding.gtf gencode_protein_coding.info.tsv
	python3 !{baseDir}/bin/gtf.py to-info  protein_coding_and_all_lncRNA.gtf protein_coding_and_all_lncRNA.info.tsv
	
	#add info
	python3 scripts/gene_info_add.py add-gffcompre        all_lncRNA.info        gffcmp.all_lncRNA.gtf.tmap all_lncRNA.info.gffcmp
	python3 scripts/gene_info_add.py add-bedtool-closest  all_lncRNA.info.gffcmp cloest.txt                 all_lncRNA_info.tsv
	python3 scripts/gene_info_add.py add-geneName all_lncRNA_info.tsv gencode_protein_coding.info.tsv all_lncRNA_info.tsv2
	mv all_lncRNA_info.tsv2 all_lncRNA_info.tsv
	cat all_lncRNA_info.tsv |awk 'NR==1 || $1 ~/^TU/' >novel_lncRNA_info.tsv
	'''
}


process classify_lncRNA {
    publishDir "${params.outdir}/novel_lncRNA/classify", mode: 'link'
    errorStrategy 'ignore'
    
    input:
    file "scripts/*" from scripts_2.collect()
    file "all_lncRNA.gtf" from all_lnc_gtf
    //file "" from protein_coding_gtf
    // actually, I only need "gencode_protein_coding.gtf"
    set "gencode_protein_coding.gtf", "known_lncRNA.gtf","known_lncRNA.list" from known_gtf_list_ch_3;
    
    output:
    set "lncRNA_class_closest_PCG.tsv", "novel_lncRNA_info.tsv"
    
	shell:
    '''
	set +u; source activate lncRNA; set -u
    awk '$3=="exon" && $1 ~/^chr|^\\d/' all_lncRNA.gtf  > all_lncRNA.exon.gtf
    awk '$3=="exon"' gencode_protein_coding.gtf > gencode_protein_coding.exon.gtf
    
    gtf2bed < all_lncRNA.exon.gtf |cut -f1-6 > novel_final.bed.a
    gtf2bed < gencode_protein_coding.exon.gtf |cut -f1-6 >annot.bed.b
    sort -k1,1 -k2,2n novel_final.bed.a >novel_final.bed
    sort -k1,1 -k2,2n annot.bed.b >annot.bed

    python3 !{baseDir}/bin/gtf.py  to-gene-bed gencode_protein_coding.exon.gtf annot_gene.bed.a
    python3 !{baseDir}/bin/gtf.py  to-gene-bed all_lncRNA.exon.gtf  novel_final_gene.bed.b

    sort -k1,1 -k2,2n annot_gene.bed.a > annot_gene.bed
    sort -k1,1 -k2,2n novel_final_gene.bed.b >novel_final_gene.bed

    bedtools intersect -a novel_final_gene.bed -b annot_gene.bed -wa -wb -v  >non_overlaps.txt
    bedtools intersect -a novel_final_gene.bed -b annot_gene.bed -wa -wb     >overlaps.txt

    awk '$6==$12 {print}' FS="\t" OFS="\t" overlaps.txt >sense_or_intron.txt
    awk '$6!=$12 {print}' FS="\t" OFS="\t" overlaps.txt >antisense.txt

    python3 !{baseDir}/bin/venn.py venn novel_final_gene.bed sense_or_intron.txt novel_final_sense_intron.bed tmp  3 3

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
    python3 !{baseDir}/bin/gtf.py to-info  all_lncRNA.exon.gtf all_lncRNA.info
	python3 !{baseDir}/bin/gtf.py to-info  gencode_protein_coding.gtf gencode_protein_coding.info.tsv
   
    ## output "lncRNA_class_closest_PCG.tsv"
    python3 scripts/lncRNA_classify.py final lncRNA_classification.txt \
    closest.txt closest_sense.txt all_lncRNA.info gencode_protein_coding.info.tsv
    
    cat lncRNA_class_closest_PCG.tsv |awk 'NR==1 || $1 ~/^TU/' > novel_lncRNA_info.tsv
    '''
}


params.datainfo = ""
cal_deg_ch.into{cal_deg_ch;cal_deg_ch_1;cal_deg_ch_2 }

if( params.datainfo && params.featureCounts ){
	datainfo = file(params.datainfo)
	process DEG_by_DESeq2_with_featureCounts {
		publishDir "${params.outdir}/deg_DESeq2/featureCounts/", mode: 'link'
		
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
		python3 ${baseDir}/bin/gtf.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
		Rscript ${baseDir}/bin/deseq2.R merged_gene_counts.txt
		
		mkdir -p files_cache__
		ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
		"""
	}

}


if (params.reads){
	
	process salmon_index{
		
		publishDir "${params.outdir}/salmon/salmon_index/", mode: 'link', 
    		saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
		input:
		file 'transcripts.fa' from cal_expr_fa.collect()
		
		output:
		file "salmon_index/*" into salmon_index_ch
		
		"""
		# set +u; source activate lncRNA; set -u
		salmon index -t transcripts.fa -i salmon_index -p 4
		"""
	}

	process salmon {
		
		publishDir "${params.outdir}/salmon/samples", mode: 'link',
    		saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
		
		tag {id } 
		input:
		set id, file(reads) from reads_salmon
		file "transcripts_index/*" from salmon_index_ch.collect()
		
		output:
		file "${id}"  into salmon_slueth,salmon_res
		 
		script:
		if(params.single){
						
			"""
			# set +u; source activate lncRNA; set -u
			salmon quant -i transcripts_index -l A -r $reads -p 8 -o ${id} --gcBias

			"""
		}else{
						
			"""
			# set +u; source activate lncRNA; set -u
			salmon quant -i transcripts_index -l A -1 ${reads[0] } -2 ${reads[1] } -p 8 -o ${id} --gcBias

			"""
		}
		
	}

	process salmon_merge {
		
		publishDir "${params.outdir}/salmon/", mode: 'link'
		
		input:
		file "ref_gtf" from cal_expr_gtf_1.collect()
		file "samples/*" from salmon_res.collect()
		
		output:
		file "salmon_gene_est_counts.tsv" into salmon_tab
		// file "files_cache__"
		file "*tsv"
		
		"""
		python3 ${baseDir}/bin/lncRNA.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
		Rscript ${baseDir}/bin/tximport_salmon.R
		# mkdir -p files_cache__
		# ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
		"""
	}

	salmon_tab.into{ salmon_tab; salmon_tab_1 }
	
	
	if(params.datainfo){
		process differential_anlysis_by_slueth{
			publishDir "${params.outdir}/deg_slueth/", mode: 'link'
			
			input:
			file "samples/*" from kallisto_slueth.collect()
			file "datainfo.txt" from datainfo
			output:
			file "*"
			
			"""
			echo test
			"""
		
		}
		
		process DEG_by_DESeq2_with_salmon {
			publishDir "${params.outdir}/deg_DESeq2/salmon/", mode: 'link'
			
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
			python3 ${baseDir}/bin/gtf.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
			Rscript  ${baseDir}/bin/deseq2.R kallisto_gene_est_counts.tsv
			
			mkdir -p files_cache__
			ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
			"""
	 }
	
	}
	
}

