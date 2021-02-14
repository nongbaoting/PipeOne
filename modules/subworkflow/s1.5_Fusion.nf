#!/usr/bin/env nextflow
/*  @Date         : 2021/01/31 16:13:14
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + '/s1.5_Fusion'
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false
params.blacklisted =  params.genome ? params.genomes[ params.genome ].blacklisted ?:false :false
params.proteinDomains =  params.genome ? params.genomes[ params.genome ].proteinDomains ?:false :false
params.cytobands =  params.genome ? params.genomes[ params.genome ].cytobands ?:false :false

params.ctat_dir =  params.genome ? params.genomes[ params.genome ].ctat_dir ?:false :false


// include funtions must be placed after params !
include { Header; check_file; input_reads; get_base_index;get_dir_files } from "${baseDir}/modules/functions.nf"
include { arriba; merge_arriba } from "${baseDir}/modules/process/fusion/arriba.nf"
include { mark_feature } from "${baseDir}/modules/modules.nf"

//include {  } from  './modules/modules.nf'

star_indices = get_dir_files(params.star_index, '--star_index')
fasta = check_file(params.fasta, '--fasta')
blacklisted = check_file(params.blacklisted, '--blacklisted')
proteinDomains = check_file(params.proteinDomains, '--proteinDomains')
cytobands = check_file(params.cytobands, '--cytobands')


workflow Fusion {
    take: 
    reads
    gtf 

    main:
    arriba(reads, star_indices.collect(), gtf, fasta, blacklisted, proteinDomains, cytobands )
    merge_arriba(arriba.out.tsv.collect() )
    mark_feature(merge_arriba.out, "fusion_arriba_out.csv", "Fusion")
}