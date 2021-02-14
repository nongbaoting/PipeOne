#!/usr/bin/env nextflow
/*  @Date         : 2021/01/31 16:38:30
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + '/s1.6_RNAediting'
params.sprint_index = params.genome ? params.genomes[ params.genome ].sprint_index ?:false :false
params.sprint_repeat = params.genome ? params.genomes[ params.genome ].sprint_repeat ?:false :false

params.snp = params.genome ? params.genomes[ params.genome ].snp ?:false :false
params.annovar_data_dir = params.genome ? params.genomes[ params.genome ].annovar_data_dir ?:false :false


params.miniEditing_reads = 10
// include funtions must be placed after params !
include { Header; check_file; input_reads; get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include{ sprint; merge_sprint_A2I} from "${baseDir}/modules/process/RNAediting/SPRINT.nf"
include {annovar_sprint } from "${baseDir}/modules/process/annotation/ANNOVAR.nf"
include { mark_feature } from "${baseDir}/modules/modules.nf"
//include {  } from  './modules/modules.nf'


sprint_repeat = check_file(params.sprint_repeat, '--sprint_repeat')
annovar_data  = get_dir_files(params.annovar_data_dir, '--annovar_data_dir')
if ( params.sprint_index  ){
    sprint_indices = Channel
        .fromPath("${params.sprint_index}*")
        .ifEmpty { exit 1, "sprint_index not found: ${params.sprint_index}" }
    sprint_base = params.sprint_index.split('/')[-1]
	
}
annovar_BinDir = Channel.fromPath("${params.annovar_BinDir}/*.pl")

workflow RNAediting {
    take:
    reads

    main:
    sprint(reads,sprint_base, sprint_indices.collect(), sprint_repeat )
    merge_sprint_A2I(sprint.out.tsv.collect() )
    annovar_sprint(merge_sprint_A2I.out.collect(), annovar_data.collect() )
    mark_feature(annovar_sprint.out.table, "RNA-editing-rate.csv", "RNAEditing" )
}