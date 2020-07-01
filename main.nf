params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------

    idplot: compare similar sequences to a reference
    ================================================

    Output is written to <outdir>/idplot.html and can
    be opened with your internet browser.

    required
    --------
    --reference  Fasta sequence that will act as the root sequence
    --fasta      One ('my.fasta') or multiple ('*.fasta') query sequences
                 to compare to `--reference`.

    options
    -------
    --outdir     Base results directory for output.
                 Default: '/.results'
    --window     The sliding window size across the reference genome upon
                 which to calculate similarity.
                 Default: 1000

    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// required arguments
params.reference = false
if( !params.reference ) { exit 1, "--reference is not defined" }
file reference = file(params.reference)
if( !reference.exists() ) { exit 1, "Reference [${reference}] does not exist." }
params.fasta = false
if( !params.fasta ) { exit 1, "--fasta is not defined" }

Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set { minimap_ch }

process minimap2 {
    tag "${reference.baseName}"

    input:
    file(reference)
    file(query) from minimap_ch.collect()

    output:
    file("${reference.baseName}.aln") into report_ch

    script:
    """
    minimap2 --cs=long ${reference} ${query} | paftools.js view -l 100000000000 -f aln - > ${reference.baseName}.aln
    """
}

process process_aln {
    publishDir path: "${params.outdir}/"

    input:
    file(aln) from report_ch
    file(reference)

    output:
    file("idplot.html")

    script:
    template 'idplot.py'
}
