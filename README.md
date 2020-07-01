# idplot

Compare similar sequences to a reference.

Sequences are aligned to the reference using `minimap2` and are assumed to
be reasonably similar. The generated plots include a multiple sequence
alignment, aggregated gaps added to the reference by the query sequences,
and percent ID per query relative to the reference sequence.

See the example plot at: https://brwnj.github.io/idplot/

### Multiple sequence alignment

The reference sequence is fully colored in. Hovering along the reference
shows the base for a given color. Query sequences are colored at mismatches
and gaps (gray).

### Gaps

The second row represents gaps added into the reference and is the sum of
all gaps across the query sequences. Zooming is allowed on both the x and
y planes to facilitate actually seeing data aside from very large gaps.

### Percent ID (ANI)

Percent ID is calculated across the window (default 1000 bp) with the value
being plotted at the center point. A 1000 bp window will have 500 bp dead
spots at the beginning and end of the reference length.

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
