#!/usr/bin/env nextflow
/*  @Date         : 2021/01/31 19:33:26
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.genome_build = params.genome ? params.genomes[ params.genome ].genome_build ?: false : false
params.outdir_sub = params.outdir + '/s1.8_SNP'
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.fasta_fai = params.genome ? params.genomes[ params.genome ].fasta_fai ?: false : false
params.fasta_dict = params.genome ? params.genomes[ params.genome ].fasta_dict ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false
params.annovar_data_dir = params.genome ? params.genomes[ params.genome ].annovar_data_dir ?:false :false

// include funtions must be placed after params !
include { Header; check_file; input_reads; 
    get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include {GATK }    from "${baseDir}/modules/process/SNP/GATK.nf"
include {STAR_AS } from "${baseDir}/modules/process/mapping/STAR.nf"
include {ANNOVAR_VCF } from "${baseDir}/modules/process/annotation/ANNOVAR.nf"
include { mark_feature } from "${baseDir}/modules/modules.nf"
//include {  } from  './modules/modules.nf'

// define files
gtf = check_file(params.gtf, '--gtf')
fasta = check_file(params.fasta, '--fasta')
fasta_fai = check_file(params.fasta_fai, '--fasta_fai')
fasta_dict = check_file(params.fasta_dict, '--fasta_dict')
star_indices = get_dir_files(params.star_index, '--star_index')
annovar_data  = get_dir_files(params.annovar_data_dir, '--annovar_data_dir')
annovar_BinDir = Channel.fromPath("${params.annovar_BinDir}/*.pl")

workflow SNPbam {
    take:
    bam
    
    main:
    GATK(bam, fasta, fasta_fai, fasta_dict)
    ANNOVAR_VCF( GATK.out.id_vcf_idx, annovar_data.collect() )
    mark_feature(ANNOVAR_VCF.out.snp_geneBase, "snp.geneBase.csv", "SNP")

}

workflow SNPfastq{
    take: 
    reads
   
    main:
    STAR_AS(reads, fasta, gtf, star_indices.collect() )
    SNPbam( STAR_AS.out.id_bam)

}