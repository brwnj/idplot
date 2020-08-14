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
    --gard       Run GARD for breakpoint detection in addition to 3seq.
                 Default: false
    --cpus       Threads for multi-threaded processes.
                 Default: 1
    --mpi        Breakpoint detection only takes advantage of more than 2
                 CPUs when MPI is available (Docker, local installs, but
                 not Singularity).
                 Default: false

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
    .set { mafft_ch }


process mafft {
    cpus params.cpus.toInteger()

    input:
    file(reference)
    file(query) from mafft_ch.collect()

    output:
    file("${reference.baseName}.msa.fasta") into (msa_gard_ch, msa_threeseq_ch, msa_report_ch)

    script:
    """
    cat ${reference} > mafft_input.fasta
    cat ${query} >> mafft_input.fasta
    mafft --auto --thread ${task.cpus} --quiet --maxiterate 1000 --globalpair mafft_input.fasta > ${reference.baseName}.msa.fasta
    """
}


process threeseq {
    input:
    file(msa) from msa_threeseq_ch

    output:
    file("${msa.baseName}.3s.rec") into threeseq_report_ch

    script:
    """
    echo Y | 3seq -g pvaltable500 500
    echo Y | 3seq -f ${msa} -ptable pvaltable500 -L500 -id ${msa.baseName}
    """
}


process gard {
    tag "${msa}"
    cpus params.cpus.toInteger()

    input:
    file(msa) from msa_gard_ch

    output:
    file("${msa.baseName}.json") into gard_output_ch

    when:
    params.gard

    script:
    cmd = params.mpi ? "mpirun -np ${task.cpus} --allow-run-as-root --oversubscribe HYPHYMPI" : "hyphy"
    // rv: None, GDD, Gamma
    """
    $cmd CPU=${task.cpus} gard --type nucleotide --code Universal \
        --alignment ${msa} --rv GDD --model JTT --rate-classes 4 --output ${msa.baseName}.json
    """
}

gard_report_ch = (params.gard ? gard_output_ch : [""])

process idplot {
    publishDir path: "${params.outdir}/"

    input:
    file(msa) from msa_report_ch
    file(json) from gard_report_ch
    file(rec) from threeseq_report_ch

    output:
    file("idplot.html")

    script:
    template 'idplot.py'
}
