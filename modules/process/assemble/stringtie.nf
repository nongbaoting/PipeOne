params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded =false

process stringtie {
    tag { id }
    publishDir "${params.outdir_sub}/stringtie", mode: 'link',
        saveAs: {filename ->
			if( params.saveIntermediateFiles ){
				if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
				else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
				else "$filename"
			}else null
            
        }

    input:
    tuple val(id), path(bam) 
    path gtf
	
    output:
    path "${id}_transcripts.gtf" 
   
    script:
    def st_direction = ''
    if (params.forward_stranded && !params.unstranded){
        st_direction = "--fr"
    } else if (params.reverse_stranded && !params.unstranded){
        st_direction = "--rf"
    }
	
    """
	set +u; source activate pipeOne_lncRNA; set -u
    stringtie $bam $st_direction \\
			-o ${id}_transcripts.gtf \\
			-p $task.cpus -G $gtf
    """
}

