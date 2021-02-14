

process fastp{
    tag {id}

    publishDir "${params.outdir_sub}/fastp/", mode: 'link',
        saveAs: {filename -> 
            if(filename =~ /fastp.fq.gz/) "clean/${filename}"
            else "report/${id}.${filename}" 
            }

    input:
    tuple val(id), path(reads) 

    output:
    tuple val(id), path("${id}*.fastp.fq.gz"), emit: reads
    tuple  path("fastp.html"), path("fastp.json"), emit: log

    script:
    if( params.singleEnd ){
       
        """
        set +u; source activate pipeOne_lncRNA; set -u
        if [ "${reads}" == "${id}.fastp.fq.gz" ];
        then
            echo "input raw reads name and output clean name are identical, please check your input files!"
            exit 1
        fi
        
        fastp -i ${reads}  -o ${id}.fastp.fq.gz  -q 20 
        """
    }else
    {
        
        """
        set +u; source activate pipeOne_lncRNA; set -u
        if [ "${reads[0]}" == "${id}.R1.fastp.fq.gz" ];
        then
            echo "input raw reads name and output clean reads name are identical, please check your input files!"
            exit 1
        fi
        
        fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 
        """
    }
				
}