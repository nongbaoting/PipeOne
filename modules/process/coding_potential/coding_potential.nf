

	
cpatpath = params.cpat_dir

process cpat{
	
	publishDir "${params.outdir_sub}/coding_potential/CPAT/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	path "novel_lncRNA_candidate.fa" 
	
	output:
	path "CPAT.out" 
	
	shell:
	
	if(params.species =~ /(?i)human(?-i)/){
        '''
		set +u; source activate pipeOne_py2; set -u
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x !{cpatpath}/dat/Human_Hexamer.tsv \
                                       -d !{cpatpath}/dat/Human_logitModel.RData \
                                       -o CPAT.out
        '''
    }else if (params.species =~ /(?i)mouse(?-i)/ ){
		
        '''
		set +u; source activate pipeOne_py2; set -u
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x !{cpatpath}/dat/Mouse_Hexamer.tsv \
                                       -d !{cpatpath}/dat/Mouse_logitModel.RData \
                                       -o CPAT.out
        '''

    }else if (params.species =~/(?i)zebrafish(?-i)/){
		
        '''
		set +u; source activate pipeOne_py2; set -u
        cpat.py -g novel_lncRNA_candidate.fa \
                                       -x !{cpatpath}/dat/zebrafish_Hexamer.tsv \
                                       -d !{cpatpath}/dat/zebrafish_logitModel.RData \
                                       -o CPAT.out
        '''
    }else if (params.species =~ /(?i)Fly(?-i)/) {
        '''
		set +u; source activate pipeOne_py2; set -u
        cpat.py -g novel_lncRNA_candidate.fa  \
                                       -x !{cpatpath}/dat/fly_Hexamer.tsv \
                                       -d !{cpatpath}/dat/Fly_logitModel.RData \
                                       -o CPAT.out
        '''
    }

}

process PLEK {
	
	publishDir "${params.outdir_sub}/coding_potential/PLEK/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	path "novel_lncRNA_candidate.fa"
	
	output:
	path "PLEK.out" 
    
	"""
    #!/bin/bash
	set +u; source activate pipeOne_py2; set -u;
	PLEK.py -fasta novel_lncRNA_candidate.fa -out PLEK -thread $task.cpus
	ln -s PLEK PLEK.out
    printf "type\tPLEK_score\ttx_id\n" >PLEK.txt
    more PLEK.out |cut -f1 -d ' ' |sed 's/>//'  >> PLEK.txt
	"""
}


process cppred {
	publishDir "${params.outdir_sub}/coding_potential/CPPred/", mode: 'link', 
    	saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
	
	input:
	path "novel_lncRNA_candidate.fa" 
	
	output:
	path "CPPred.out" 
    
	"""
	set +u; source activate pipeOne_py2; set -u
    python2 ${params.CPPred_dir}/bin/CPPred.py -i novel_lncRNA_candidate.fa \
		-hex ${params.CPPred_dir}/Hexamer/Human_Hexamer.tsv -r ${params.CPPred_dir}/Human_Model/Human.range \
		-mol ${params.CPPred_dir}/Human_Model/Human.model -spe Human -o CPPred.out
	"""
}

workflow cal_coding_P{
	take:
	fasta

	main:
	cpat(fasta)
	PLEK(fasta)
	cppred(fasta)

	emit:
	cpat = cpat.out
	PLEK = PLEK.out
	cppred = cppred.out

}