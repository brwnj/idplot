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
                 Default: 500
    --gard       Run GARD for breakpoint detection in addition to 3seq.
                 Default: false
    --gff        Reference annotation file in .gff or .gff3 format.
                 Coordinates will be adjusted to include gaps introduced
                 during sequence alignment.
                 Default: false
    --cpus       Threads for multi-threaded processes.
                 Default: 1
    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}
// required arguments
params.reference = false
if( !params.reference ) { exit 1, "--reference is not defined" }
reference = file(params.reference)
if( !reference.exists() ) { exit 1, "Reference [${reference}] does not exist." }
params.fasta = false
if( !params.fasta ) { exit 1, "--fasta is not defined" }

gff = params.gff ? file(params.gff) : false

Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set { mafft_ch }


process mafft {
    publishDir path: "${params.outdir}/", mode: "copy"
    cpus params.cpus.toInteger()

    input:
    path(reference)
    path(query) from mafft_ch.collect()

    output:
    path("${reference.baseName}.msa.fasta") into (msa_gard_ch, msa_json_ch, msa_report_ch)

    script:
    """
    cat ${reference} > mafft_input.fasta
    cat ${query} >> mafft_input.fasta
    mafft --auto --thread ${task.cpus} --maxiterate 1000 --globalpair mafft_input.fasta > ${reference.baseName}.msa.fasta
    """
}


process gard {
    publishDir path: "${params.outdir}/gard/", mode: "copy"
    tag "${msa}"
    cpus params.cpus.toInteger()

    input:
    path(msa) from msa_gard_ch

    output:
    path("${msa.baseName}.json") into (gard_output_ch, gard_output_secondary_ch)

    when:
    params.gard

    script:
    cmd = params.nompi ? "hyphy" : "mpirun -np ${task.cpus} --allow-run-as-root --oversubscribe HYPHYMPI"
    // rv: None, GDD, Gamma
    """
    $cmd CPU=${task.cpus} gard --type nucleotide --code Universal \
        --alignment ${msa} --rv GDD --model JTT --rate-classes 4 --output ${msa.baseName}.json
    """
}

gard_json_ch = (params.gard ? gard_output_ch : [""])
gard_report_ch = (params.gard ? gard_output_secondary_ch : [""])

process jsontofasta {
    input:
    path(msa) from msa_json_ch
    file(json) from gard_json_ch

    output:
    path("*.fa") into selection_ch

    when:
    params.gard

    script:
    template "jsontofasta.py"
}


process fasttree {
    input:
    path(fasta) from selection_ch.flatten()

    output:
    path("${fasta.simpleName}.tree") into tree_output_ch

    script:
    """
    fasttree -nt -gamma -spr 4 -quiet ${fasta} > ${fasta.simpleName}.tree
    """
}

tree_report_ch = (params.gard ? tree_output_ch : [""])

process idplot {
    publishDir path: "${params.outdir}/", mode: "copy"

    input:
    path(msa) from msa_report_ch
    file(json) from gard_report_ch
    file(trees) from tree_report_ch.collect()
    file(gff)

    output:
    path("idplot.html")

    script:
    template "idplot.py"
}
