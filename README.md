# idplot

Compare similar sequences to a reference.

# Setup

Nextflow is used to run the pipeline. Its installation instructions
can be found at https://www.nextflow.io/ or installed via conda by
way of bioconda:

```
conda install -c conda-forge -c bioconda nextflow
```

# Usage

Executing the workflow using nextflow:

```
nextflow run brwnj/idplot -latest \
    --reference data/MN996532.fasta \
    --fasta 'data/query_seqs/*.fasta'
```

The reference sequence (`--reference`) should be a fasta with only one
sequence in it. Query sequences (`--fasta`) may either be single sequence
files or multi-sequence fasta files and you can specify more than one
using wildcards ('*').

Example sequences are found in `data/query_sequences`.

Output is written to <outdir>/idplot.html and can
be opened with your internet browser.

An example report is available at: https://brwnj.github.io/idplot/

```
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
```
