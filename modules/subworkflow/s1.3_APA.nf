#!/usr/bin/env nextflow
/*  @Date         : 2021/01/30 23:36:59
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params

params.outdir_sub = params.outdir + '/s1.3_APA'
params.utr_gtf =  params.genome ? params.genomes[ params.genome ].utr_gtf  ?: false : false
params.salmon_index_3UTR = params.genome ? params.genomes[ params.genome ].salmon_index_3UTR  ?: false : false
params.replace_SalmonIndex_ID = params.genome ? params.genomes[ params.genome ].replace_SalmonIndex_ID  ?: false : false
params.apa3utr =  params.genome ? params.genomes[ params.genome ].apa3utr ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
include { mark_feature } from "${baseDir}/modules/modules.nf"

// include funtions must be placed after params !
include { Header; check_file; input_reads; get_dir_files; get_base_index } from '../functions.nf'
include {salmon_APA; qapa; qapa_filter } from  "${baseDir}/modules/process/quant/qapa.nf"
include {salmon_idx_cal } from  "${baseDir}/modules/process/quant/salmon.nf"
//include {  } from  './modules/modules.nf'


salmon_indices = get_dir_files(params.salmon_index_3UTR, '--salmon_index_3UTR' )
replace_SalmonIndex_ID = check_file(params.replace_SalmonIndex_ID, '--replace_SalmonIndex_ID')
utr_gtf = check_file(params.utr_gtf, '--utr_gtf' )
apa3utr = check_file(params.apa3utr, '--apa3utr' )
fasta   = check_file(params.fasta,   '--fasta'   )

workflow APA {
    take: 
    reads
    TPM

    main:
    salmon_APA(reads, salmon_indices.collect() )
    qapa(utr_gtf, replace_SalmonIndex_ID, salmon_APA.out.collect() ) 
    qapa_filter(qapa.out, TPM,  utr_gtf )
    mark_feature(qapa_filter.out.pau, "APA_pau-distal-proximal.csv", "APA")
   

}

workflow APA_salmonGene {
    take: 
    reads

    main:
    salmon_idx_cal(
        reads,  
        fasta,
        utr_gtf
    )
    APA(reads, salmon_idx_cal.out.tpm )

}

