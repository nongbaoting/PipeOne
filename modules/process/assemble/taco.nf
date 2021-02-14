process taco{
	
	publishDir "${params.outdir_sub}/taco/", mode: 'link', 
    	saveAs: { filename -> if( params.saveIntermediateFiles  && filename =~ /gtf/ ) "$filename" 
							  else null }

	input:
	path gtf
	path "stringtie_gtf/*" 
	
	output:
	path "taco_out/assembly.gtf" 
	
	
	"""
	set +u; source activate pipeOne_lncRNA; set -u
	ls stringtie_gtf | awk '{print "stringtie_gtf/"\$0}' > gtf_list
	taco_run -p $task.cpus -o taco_out gtf_list
	"""
}
