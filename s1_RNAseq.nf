#!/usr/bin/env nextflow
/*  @Date         : 2020/09/23 21:30:53
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/

nextflow.enable.dsl=2
// defined params
params.reads =''
params.singleEnd = false
params.cleaned  = false
params.genome =''
params.saveIntermediateFiles = false
params.outdir = "./results"
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.update_GTF = false
params.prepare_ref = false

// include funtions must be placed after params !
include {Header; mkdir_tmp; check_file; input_reads; get_base_index } from './modules/functions.nf'
include {Prepare_References } from "./modules/subworkflow/s0_prepare_ref.nf"
include {QC; } from './modules/subworkflow/s1.0_QC.nf'
include {mRNA_lncRNA } from  './modules/subworkflow/s1.1_mRNA_lncRNA.nf'
include {CircRNA_CIRIquant_bam; CircRNA_CIRIquant_fastq} from "./modules/subworkflow/s1.2_circRNA.nf"
include {APA; } from  './modules/subworkflow/s1.3_APA.nf'
include {RetroTrans} from  './modules/subworkflow/s1.4_RetroTrans.nf'
include {Fusion} from "./modules/subworkflow/s1.5_Fusion.nf"
include {RNAediting} from "./modules/subworkflow/s1.6_RNAediting.nf"
include {AS} from "./modules/subworkflow/s1.7_AS.nf"
include {SNPbam; SNPfastq} from "./modules/subworkflow/s1.8_SNP.nf"

log.info Header()
mkdir_tmp('./tmp')
// or: mRNA_lncRNA,circRNA,APA,RetroTrans,Fusion,RNAediting,AS,SNP  or mix: 1,2,RetroTrans,Fusion,RNAediting,AS,8
params.run_s1 = '1,2,3,4,5,6,7,8' 

def runList = []
def gtf
def ch_reads

workflow {

    // define prepare or run
    if( params.prepare_ref ){
        println("Preparing reference data!")
        Prepare_References()
        
    }else{
        println("RNA-seq data processing!")
        gtf     = check_file(params.gtf, '--gtf') 
        ch_reads = input_reads(params.reads, params.singleEnd )
        runList = params.run_s1.toString().split(',')
        QC(ch_reads)
        reads = QC.out.reads
    

        if( runList.contains('1') || runList.contains('mRNA_lncRNA') ){  
            mRNA_lncRNA(reads)
            if(params.update_GTF){
                gtf = mRNA_lncRNA.out.gtf
            }
        } 

        if( runList.contains('2') || runList.contains('circRNA' )    ){

            if( runList.contains('1') || runList.contains('mRNA_lncRNA') ){
                CircRNA_CIRIquant_bam( reads, mRNA_lncRNA.out.bam, gtf)
            }else{
                CircRNA_CIRIquant_fastq(reads, gtf )
            }
        } 

        if( runList.contains('3') || runList.contains('APA' )         ){ 
            if( runList.contains('1') || runList.contains('mRNA_lncRNA') ){
                APA(reads, mRNA_lncRNA.out.tpm)
            }else{
                APA_salmonGene(reads )
            }
        } 

        if( runList.contains('4') || runList.contains('RetroTrans')  ){ RetroTrans(reads) }

        if( runList.contains('5') || runList.contains('Fusion')      ){ Fusion(reads, gtf) }

        if( runList.contains('6') || runList.contains('RNAediting')  ){ RNAediting(reads) }

        if( runList.contains('7') || runList.contains('AS')          ){ AS(reads, gtf) }

    
        if( runList.contains('8') || runList.contains('SNP')         ){ 

            if( runList.contains('7') || runList.contains('AS') ){
                SNPbam(AS.out.id_bam)
            }else{
                SNPfastq(reads)
            }
        }

    }

}