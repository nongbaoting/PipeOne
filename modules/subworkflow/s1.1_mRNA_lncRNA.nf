#!/usr/bin/env nextflow
/*  @Date         : 2021/01/30 17:05:49
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting

*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir  + "/s1.1_mRNA_lncRNA"
params.genome = ""
params.species = params.genome ? params.genomes[ params.genome ].species ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.fasta_fai = params.genome ? params.genomes[ params.genome ].fasta_fai ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.lncpedia_gtf  = params.genome ? params.genomes[ params.genome ].lncpedia_gtf  ?: false : false
params.genecode_lncRNA_gtf  = params.genome ? params.genomes[ params.genome ].genecode_lncRNA_gtf   ?: false : false
params.hisat2_idx = params.genome ? params.genomes[ params.genome ].hisat2_index ?:false :false

params.featureCounts = false
params.saveIntermediateFiles = false

// include funtions must be placed after params !
include { Header; check_file; input_reads; get_base_index } from "${baseDir}/modules/functions.nf"
include { mark_feature as mark_feature_lnc ; mark_feature as mark_feature_prot } from "${baseDir}/modules/modules.nf"
include {hisat2; } from  "${baseDir}/modules/process/mapping/hisat2.nf"
include {stringtie; } from  "${baseDir}/modules/process/assemble/stringtie.nf"
include {taco; } from  "${baseDir}/modules/process/assemble/taco.nf"
include {gffcompare; prepare_ref_gtf; filter_coding_potentail; format_lncRNA_info;
         classify_lncRNA; epxr_gene_summary; sep_lnc_mRNA } from "${baseDir}/modules/process/lncRNA/mRNA_lncRNA.nf"
include { cal_coding_P } from "${baseDir}/modules/process/coding_potential/coding_potential.nf"
include {salmon_idx_cal } from  "${baseDir}/modules/process/quant/salmon.nf"
//include {  } from  './modules/modules.nf'

// check file 
(hisat2_base, hisat2_indices) = get_base_index(params.hisat2_idx )
gtf     = check_file(params.gtf, '--gtf')
lncpedia_gtf = check_file(params.lncpedia_gtf, '--lncpeida_gtf')
genecode_lncRNA_gtf = check_file(params.genecode_lncRNA_gtf, '--genecode_lncRNA_gtf')
fasta    = check_file(params.fasta, '--fasta')
fasta_fai = check_file(params.fasta_fai, '--fasta_fai')

workflow mRNA_lncRNA {
    take:
    reads

    main:
    prepare_ref_gtf(gtf, lncpedia_gtf, genecode_lncRNA_gtf)
    hisat2(reads, hisat2_base, hisat2_indices )
    stringtie(hisat2.out.bam, gtf)
    taco(gtf, stringtie.out.collect() )
    gffcompare(fasta, fasta_fai, gtf, lncpedia_gtf, taco.out )
    cal_coding_P(gffcompare.out.nov_lnc_candidate_fa)
    filter_coding_potentail(
            gtf, fasta, fasta_fai, 
            cal_coding_P.out.cpat,
            cal_coding_P.out.PLEK, 
            cal_coding_P.out.cppred ,
            prepare_ref_gtf.out.known_gtf_lst.collect(), 
            gffcompare.out.nov_lnc_candidate.collect() )

    format_lncRNA_info(
        filter_coding_potentail.out.nov_lst.collect(),
        prepare_ref_gtf.out.known_gtf_lst.collect(),
        filter_coding_potentail.out.cal_expr_gtf )

    classify_lncRNA(
        filter_coding_potentail.out.all_lnc_gtf,
        prepare_ref_gtf.out.known_gtf_lst.collect() )

    salmon_idx_cal(
        reads,  
        fasta,
        filter_coding_potentail.out.cal_expr_gtf )

    epxr_gene_summary(
        salmon_idx_cal.out.expr_gene_id,
        filter_coding_potentail.out.cal_expr_gtf,
        filter_coding_potentail.out.nov_lnc_gtf )

    sep_lnc_mRNA( 
        salmon_idx_cal.out.keep_TPM, 
        format_lncRNA_info.out.all_lnc_info )
    
    mark_feature_lnc(sep_lnc_mRNA.out.lnc_TPM,  'lncR_gene.tpm.csv',  'lncRNA')
    mark_feature_prot(sep_lnc_mRNA.out.prot_TPM,'prot_gene.tpm.csv', 'mRNA')

    emit:
    bam = hisat2.out.bam
    tpm = salmon_idx_cal.out.TPM
    gtf = epxr_gene_summary.out.gtf

}



