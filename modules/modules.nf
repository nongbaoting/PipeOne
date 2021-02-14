

process  mark_feature {
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path infile
    val outfile_name
    val marks

    when: infile.size() > 1
    
    output:
    path "${outfile_name}" optional true

    
    """
    python3 ${baseDir}/bin/RNAseq/summary_table.py mark_feature \\
        $infile  \\
        $outfile_name \\
        $marks
    """

}

