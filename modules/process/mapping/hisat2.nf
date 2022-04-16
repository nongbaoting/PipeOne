params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded =false

process hisat2 {
		label 'bigCPU'
        tag {id}
        publishDir "${params.outdir_sub}/hisat2/", mode: 'copy',
            saveAs: {filename -> 
                if(filename =~ /bam/ && params.saveIntermediateFiles  )  "bam/$filename"
                else if (filename =~/log/) "logs/${filename}"
				else null
                }
        
        input:
        tuple val(id), path(reads)
        val hisat2_base
        path "hisat2_index/*"
        
        output:
        tuple val(id), path("${id}.hisat2.sortbycoordinate.bam") , emit: bam
        path "${id}.hisat2.log", emit: log
        
        script:
        def rnastrandness = ''
        if (params.forward_stranded && !params.unstranded){
            rnastrandness =params.singleEnd ?  '--rna-strandness F' :  '--rna-strandness FR'
        } else if (params.reverse_stranded && !params.unstranded){
            rnastrandness = params.singleEnd ? '--rna-strandness R'  :  '--rna-strandness RF'
        }
        
        if(! params.singleEnd){
            """
			set +u; source activate pipeOne_lncRNA; set -u
            hisat2 -p $task.cpus --dta  $rnastrandness  -x hisat2_index/$hisat2_base \
            -1 ${reads[0]} -2 ${reads[1]}   2> ${id}.hisat2.log | samtools sort -@ 2 - -o ${id}.hisat2.sortbycoordinate.bam 
            samtools index ${id}.hisat2.sortbycoordinate.bam
			conda deactivate
            """
        }else{
            
            """
			set +u; source activate pipeOne_lncRNA; set -u
            hisat2 -p $task.cpus --dta  $rnastrandness -x  hisat2_index/$hisat2_base \
            -U  ${reads}   2> ${id}.hisat2.log  | samtools sort -@ 2 - -o ${id}.hisat2.sortbycoordinate.bam
            samtools index ${id}.hisat2.sortbycoordinate.bam
			conda deactivate
            """
        }
 
    }