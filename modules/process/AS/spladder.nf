
placeholder_bam = file("${baseDir}/dat/placeholder.bam")
placeholder_bai = file("${baseDir}/dat/placeholder.bam.bai")
process generate_gtf_pickle {
	
	input:
	path gtf
	path test_bam
	path test_bai

	output:
	path "${gtf}.pickle" 
	
	"""
	set +u; source activate pipeOne_AS; set -u
	spladder build -o spladder_out \\
                   -a ${gtf} \\
                   -b ${test_bam} \\
                   --merge-strat single \\
                   --no-extract-ase \\
				   --parallel ${task.cpus} \\
				   -n 100 \\
				   -c 3
	"""
}

process Single_graphs {
	label 'moreParralle'
	tag {id}
	
	input:
	tuple val(id), path(bam), path(bai) 
	path gtf
	path "${gtf}.pickle" 
	
	output:
	path "spladder_out/spladder/genes_graph_conf3.${id}.pickle", emit: graph
	tuple val(id), path("spladder_out/spladder/genes_graph_conf3.${id}.pickle"), emit: id_graph
	"""
	set +u; source activate pipeOne_AS; set -u
	spladder build -o spladder_out \\
                   -a ${gtf} \\
                   -b ${bam} \\
                   --merge-strat single \\
                   --no-extract-ase \\
				   --parallel ${task.cpus } \\
				   -n 100 \\
				   -c 3
	"""
}

process Merged_graph {
	
	input:
	path "spladder_out/spladder/*" 
	path "*"  
	path gtf
	path "${gtf}.pickle" 
	
	output:
	path "spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle"
	//path "*.hdf5" into mergeGraph_hdf5_ch, hdf5_all
	
	"""
	set +u; source activate pipeOne_AS; set -u
	ls |grep "bam\$" > alignments.txt
	mkdir -p  spladder_out/tmp
	
	spladder build -o spladder_out \\
               -a ${gtf} \\
               -b alignments.txt \\
               --merge-strat merge_graphs \\
               --no-extract-ase  \\
			   --parallel 48 \\
			   -n 100 \\
				-c 3
	"""
}


process Quant_each {
	tag {id}
	label 'moreParralle'
	input:
	
	path "spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle" 
	//set id, path(bam), path(hdf5),"spladder_out/spladder/genes_graph_conf3.${id}.pickle" from quant_ch
	tuple val(id), path(bam), path("spladder_out/spladder/genes_graph_conf3.${id}.pickle")
	path gtf
	path gtf_pickle
	
	output:
	path "spladder_out/spladder/genes_graph_conf3.merge_graphs.${id}.count.hdf5" 
	
	
	"""
	set +u; source activate pipeOne_AS; set -u
	spladder build -o spladder_out -a ${gtf} -b ${id}.bam \\
                   --merge-strat merge_graphs \\
                   --no-extract-ase \\
                   --quantify-graph \\
                   --qmode single  --parallel 10 \\
				   -n 100 \\
				   -c 3
	"""
}


process QuantAll_and_EventCalling{
	
	publishDir "${params.outdir_sub}/spladder_out/", mode: 'link'
	
	input:
	path "*" //bam_bai
	path "spladder_out/spladder/*" 
	path "spladder_out/spladder/*" 
	path "spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle" 
	path gtf
	path gtf_pickle
	
	output:
	path "spladder_out/merge_graphs*" 
	
	"""
	set +u; source activate pipeOne_AS; set -u
	ls |grep "bam\$" > alignments.txt
	mkdir -p  spladder_out/tmp
	
	echo "merge counts ...."
	spladder build -o spladder_out \\
               -a ${gtf} \\
               -b alignments.txt \\
               --merge-strat merge_graphs \\
               --no-extract-ase \\
               --quantify-graph \\
               --qmode collect --parallel 48 \\
			   -n 100 \\
			   -c 3
			   
	echo "Event Calling ...."		   
	spladder build -o spladder_out \\
		-a ${gtf} \\
		-b alignments.txt \\
		--parallel ${task.cpus } \\
		--event-types  exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip \\
		--output-txt 
			
	"""
}

process PSI_table{
	
	publishDir "${params.outdir_sub}/spladder_out_table/", mode: 'link'

	input:
	path "spladder_out/*" 
	
	output:
	path "*psi.txt.gz" 
	
	shell:
	''' 
	set +u; source activate pipeOne_py3; set -u
	mkdir -p  spladder_PSI/
	for i in `ls spladder_out/*.confirmed.txt.gz`
	do
		python3 !{baseDir}/bin/RNAseq/alternative_splicing.py spladder_PSI_tab ${i}
	done
	'''
}

process  mark_feature_spladder {
    publishDir "${params.outdir}/tables", mode: 'link'

    input:
    path "mydir/*"

    output:
    path "*.confirmed.psi.csv"

    
    """
    python3 ${baseDir}/bin/RNAseq/summary_table.py mark_feature_spladder mydir/ 
    """

}


workflow spladder {
    take:
    id_bam_bai
    gtf

    main:
    generate_gtf_pickle(gtf, placeholder_bam, placeholder_bai)
    Single_graphs(id_bam_bai, gtf, generate_gtf_pickle.out )
    bam_bai = id_bam_bai.map{ it -> a =  it[1 .. 2 ] }.collect()
    id_BamBai = id_bam_bai.map{ it -> a = [it[0], it[1 .. 2 ]] }
    id_BamBai_graph = id_BamBai.combine( Single_graphs.out.id_graph, by: 0 )
    Merged_graph(
        Single_graphs.out.graph.collect(),
        bam_bai, 
        gtf,
        generate_gtf_pickle.out 
        )
    
    Quant_each(
        Merged_graph.out,
        id_BamBai_graph,
        gtf,
        generate_gtf_pickle.out
        )
    QuantAll_and_EventCalling(
        bam_bai,
        Single_graphs.out.graph.collect(),
        Quant_each.out.collect(),
        Merged_graph.out,
        gtf,
        generate_gtf_pickle.out
    )

    PSI_table(QuantAll_and_EventCalling.out.collect() )

	mark_feature_spladder(PSI_table.out.collect() )
}