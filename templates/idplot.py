#!/usr/bin/env python

import json
from collections import deque
from itertools import groupby, islice

import click
import numpy as np

TEMPLATE = """
<!DOCTYPE html>
<html>

<head>
    <meta charset="utf-8" />
    <meta name="author" content="Joe Brown" />
    <title>idplot</title>
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>
    <!-- datatables, bootstrap 4 with select, responsive, scroller, searchpanes and html5 buttons -->
    <script type="text/javascript"
        src="https://cdn.datatables.net/v/bs4/dt-1.10.20/b-1.6.1/b-html5-1.6.1/r-2.2.3/sc-2.0.1/sp-1.0.1/sl-1.3.1/datatables.min.js"></script>
    <script type="text/javascript"
        src="https://cdnjs.cloudflare.com/ajax/libs/selectize.js/0.12.6/js/standalone/selectize.min.js"></script>
    <script type="text/javascript"
        src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js"></script>
    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>

    <link rel="stylesheet" type="text/css"
        href="https://cdn.datatables.net/v/bs4/dt-1.10.20/b-1.6.1/b-html5-1.6.1/r-2.2.3/sc-2.0.1/sp-1.0.1/sl-1.3.1/datatables.min.css" />
    <link rel="stylesheet" type="text/css"
        href="https://cdnjs.cloudflare.com/ajax/libs/selectize.js/0.12.6/css/selectize.bootstrap3.min.css">
    <link rel="stylesheet" type="text/css"
        href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">

    <style type="text/css">
        .disabled_div {
            pointer-events: none;
            opacity: 0.4
        }
    </style>
</head>

<body>
    <div class="container-fluid w-90 pt-3">
        <div class="tab-content">
            <div class="tab-pane fade show active" id="chrom_coverage" role="tabpanel" aria-labelledby="scaled_tab">
                <div class="row">
                    <div class="col-12">
                        <h3 class="d-none d-xl-block">
                            Reference: {{target}}
                        </h3>
                        <h3 class="d-xl-none">
                            Ref: {{target}}
                        </h3>
                    </div>
                </div>
                <div class="tab-content">
                    <div class="tab-pane fade show active mb-3" id="scaled" role="tabpanel"
                        aria-labelledby="scaled_tab">
                        <div style="height:600px" id="scaled_plot_placeholder">
                            <div
                                class="mt-3 d-flex justify-content-center align-items-center bg-light text-muted h-100">
                                <div class="d-flex flex-column">
                                    <div>
                                        Loading...
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div class="row pt-3 mb-3" id="scaled_plot" hidden></div>
                    </div>
                    <div class="tab-pane fade mb-3" id="cov" role="tabpanel" aria-labelledby="cov_tab">
                        <div style="height:600px" id="cov_plot_placeholder">
                            <div
                                class="mt-3 d-flex justify-content-center align-items-center bg-light text-muted h-100">
                                <div class="d-flex flex-column">
                                    <div>
                                        Loading...
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div class="row pt-3 mb-3" id="cov_plot" hidden></div>
                    </div>
                </div>
            </div>
        </div>
        <div class="container-fluid mb-5">
            <div class="row pt-3 pb-3">
                <div class="col-12">
                    <div class="table-responsive">
                        <table id="ped_table" class="table table-hover table-striped table-sm display nowrap"
                            width="100%"></table>
                    </div>
                </div>
            </div>
        </div>
    </div>
</body>

<script>
    const data = {{data}}
    const cov_color = 'rgba(108,117,125,0.2)'

    const plot_layout = {
        title: "",
        margin: { t: 10, b: 40 },
        height: 600,
        xaxis: { title: "Position", autorange: true, showgrid: false, showlines: false, zeroline: false },
        yaxis3: { title: "ANI", fixedrange: true, range: [0,1], showgrid: true, showticklabels: true, tickmode: 'array', tick0: 0, dtick: 0.2, zeroline: true, domain: [0, 0.50] },
        yaxis2: { title: "Gaps", fixedrange: false, showgrid: true, showticklabels: true, domain: [0.55, 0.75] },
        yaxis: { fixedrange: true, showgrid: true, domain: [0.75, 1]},
        hovermode: "closest",
        showlegend: false,
        grid: {rows: 3, columns: 1, subplots:[['xy'], ['xy2'], ['xy3']], pattern: 'independent',},

    }
    var ped = false
    var cov_traces = []
    var scaled_traces = []
    var scatter_point_color = 'rgba(31,120,180,0.5)'
    var scatter_point_color_light = 'rgba(255,255,255,0.1)'
    var gene_search_obj
    var ped_table

    const build_plots = (arr) => {
        // hide the placeholder
        \$('#scaled_plot_placeholder').prop('hidden', true)
        // show the plot
        \$('#scaled_plot').prop('hidden', false)
        plot_layout.xaxis.autorange = true
        scaled_traces = arr

        let scaled_plot = document.getElementById("scaled_plot")
        Plotly.react(scaled_plot, scaled_traces, plot_layout)
        scaled_plot.removeAllListeners("plotly_click")
        scaled_plot.removeAllListeners("plotly_doubleclick")
        scaled_plot.on("plotly_click", handle_plot_click)
        scaled_plot.on("plotly_doubleclick", handle_plot_doubleclick)
        \$("#scaled_plot").removeClass("disabled_div")
    }

    const handle_plot_click = (click_data) => {
        if (click_data.points[0].data.tracktype == 'gff') {
            let genes = click_data.points[0].text.split(";")
            for (var i = 0; i < genes.length; i++) {
                window.open("https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + genes[i], "_blank")
            }
        } else if (click_data.points[0].data.tracktype == 'vcf') {
            // nothing yet
            return
        } else if (click_data.points[0].data.tracktype == 'bed') {
            // nothing yet
            return
        } else {
            let sample_id = click_data.points[0].data.text
            if (sample_id) {
                highlight_plot_traces(sample_id)
                if (ped) {
                    search_datatable(sample_id)
                }
            }
        }
    }

    const handle_scatter_selection = (event) => {
        // plot click event
        if (event === undefined) {
            return
        }
        var sample_ids = event.points.map(point => point.text)
        var colors = []
        for (var i = 0; i < data.depth.bins.samples.length; i++) {
            colors.push(scatter_point_color_light)
        }
        event.points.forEach((pt) => {
            colors[pt.pointNumber] = scatter_point_color
        })

        Plotly.restyle('inferred_sex', 'marker.color', [colors])
        Plotly.restyle('bin_counts', 'marker.color', [colors])
        if ("pca_1" in data.depth.pca) {
            Plotly.restyle('pca_1', 'marker.color', [colors])
            Plotly.restyle('pca_2', 'marker.color', [colors])
        }
    }

    const reset_line_plots = () => {
        plot_layout.xaxis.autorange = true
        // de-select scaled plot
        Plotly.react("scaled_plot", scaled_traces, plot_layout)
        // cov_layout.xaxis.range = [0, 1.5]
        // cov_layout.yaxis.range = [0, 1.]
        // de-select cov plot
        // Plotly.react("cov_plot", cov_traces, cov_layout)
    }

    const handle_plot_doubleclick = () => {
        reset_line_plots()
        // de-select in table
        if ("ped" in data) {
            reset_datatable()
        }
    }

    const highlight_plot_traces = (sample_id) => {
        let s_traces = []
        let c_traces = []
        let k_traces = []
        let highlight_color;
        for (var i = 0; i < scaled_traces.length; i++) {
            // let trace = scaled_traces[i]
            let trace = \$.extend(true, {}, scaled_traces[i])
            // limit to significant sample traces
            if (trace.name == "significant") {
                // de-prioritize; gray
                if (trace.text != sample_id) {
                    trace.marker.color = cov_color
                }
                else {
                    highlight_color = scaled_traces[i].marker.color
                }
            }
            s_traces.push(trace)
        }
        for (var i = 0; i < cov_traces.length; i++) {
            let trace = \$.extend(true, {}, cov_traces[i])
            if (trace.text != sample_id) {
                trace.marker.color = cov_color
            } else {
                trace.marker.color = highlight_color
            }
            c_traces.push(trace)
        }
        // Plotly.react("cov_plot", c_traces, cov_layout)
        Plotly.react("scaled_plot", s_traces, plot_layout)
    }

    \$(document).ready(function () {
        build_plots(data)
        // build_table()
        // build_global_qc()
    })

</script>

</html>
"""
COLORS = [
    "#1f77b4",  # plotly and d3 - blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # dark yellow
    "#17becf",  # teal
    "#7CB5EC",  # highcharts
    "#434348",
    "#90ED7D",
    "#F7A35C",
    "#8085E9",
    "#F15C80",
    "#E4D354",
    "#2B908F",
    "#F45B5B",
    "#91E8E1",
    "#4E79A7",  # tableau
    "#F28E2C",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC949",
    "#AF7AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AB",
]
COLORSCALE = [
    [0, "#2ca02c"],
    [0.2, "#2ca02c"],
    [0.2, "#1f77b4"],
    [0.4, "#1f77b4"],
    [0.4, "#ff7f0e"],
    [0.6, "#ff7f0e"],
    [0.6, "#d62728"],
    [0.8, "#d62728"],
    [0.8, "#7f7f7f"],
    [1, "#7f7f7f"],
]


def sliding_window(iterable, size=2, step=1, fillvalue=None):
    # iterable = list(iterable)
    # iterable.extend([fillvalue] * (size - step))
    if size < 0 or step < 1:
        raise ValueError
    it = iter(iterable)
    q = deque(islice(it, size), maxlen=size)
    if not q:
        return
    # initial padding if iterable shorter than size
    q.extend(fillvalue for _ in range(size - len(q)))
    while True:
        yield iter(q)  # iter() to avoid accidental outside modifications
        try:
            q.append(next(it))
        except StopIteration:  # Python 3.5 pep 479 support
            return
        q.extend(next(it, fillvalue) for _ in range(step - 1))


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--alignments", type=click.File("r"), default="${aln}")
@click.option("--reference", type=click.File("r"), default="${reference}")
@click.option("--window", default=${params.window}, type=int, help="window size")
@click.option("--output", type=click.File("w"), default="idplot.html", help="output file name")
def main(alignments, reference, window, output):
    """minimap2 alignments file, reference fasta, and window size"""
    reference_length = False
    reference_name = False
    mismatches = False
    query_gaps = False
    query_seq = False
    target = False

    match_library = dict()
    samples = []
    # A - green, C - blue, G - yellow, T - red
    nuc_map = {
        "a": 0,
        "c": 1,
        "g": 2,
        "t": 3,
        "u": 3,
        "A": "",
        "C": "",
        "G": "",
        "T": "",
        "-": 4,
    }

    traces = list()
    z_values = list()
    y_values = list()
    text_values = list()
    traces = list()
    query_color_index = 0

    for header, group in groupby(alignments, key=lambda x: x[0] == ">"):
        if header:
            line = next(group)
            toks = line.strip().split("\t")
            query = toks[0].strip(">")
            target = toks[5]
            if not reference_length:
                reference_length = int(toks[6])
                gaps = np.empty((reference_length,), dtype=int)
            reference_name = toks[5]
            plot_line_color = COLORS[query_color_index % len(COLORS)]
            query_color_index += 1
        else:
            click.echo(query)
            mismatches = np.ones((reference_length,), dtype=int)
            query_msa = np.empty((reference_length,), dtype=object)
            ref_seq = ""
            aln_idx = 0
            for aln in group:
                if aln.startswith("Ref+"):
                    toks = aln.split()
                    aln_idx = int(toks[1]) - 1
                    ref_seq = toks[2].strip()
                if aln.startswith("Qry+"):
                    toks = aln.split()
                    for rbase, qbase in zip(ref_seq, toks[2].strip()):
                        if qbase == "-":
                            query_msa[aln_idx] = 4
                        if rbase == "-":
                            mismatches[aln_idx] += 1
                            gaps[aln_idx] += 1
                        else:
                            if rbase.isupper():
                                mismatches[aln_idx] -= 1
                            query_msa[aln_idx] = nuc_map[qbase]
                            aln_idx += 1

            identities = list()
            for mismatch_window in sliding_window(mismatches, window):
                mw = list(mismatch_window)
                pid = max((window - sum(mw)) / window, 0)
                identities.append(pid)

            match_library[query] = identities
            trace = dict(
                x=np.arange(window / 2, len(identities) + (window / 2)).tolist(),
                y=identities,
                text=query,
                xaxis="x",
                yaxis="y3",
                connectgaps=False,
                hoverinfo="text",
                type="scatter",
                mode="lines",
                name="significant",
                # include color as primary colors occupied by area traces
                marker={"width": 1, "color": plot_line_color},
            )
            traces.append(trace)
            z_values.append(query_msa.tolist())
            y_values.append(query)
            text_values.append([])

    # reference gaps
    trace = dict(
        x=np.arange(0, len(identities)).tolist(),
        y=gaps.tolist(),
        xaxis="x",
        yaxis="y2",
        type="bar",
        # include color as primary colors occupied by area traces
        marker={"color": COLORS[7]},
    )
    traces.append(trace)

    ref_bases = []
    letters = []
    for line in reference:
        if line.startswith(">"):
            continue
        for base in line.strip():
            ref_bases.append(nuc_map[base.lower()])
            if base == "U" or base == "u":
                base = "T"
            letters.append(base.upper())
    z_values.append(ref_bases)
    y_values.append(reference_name)
    text_values.append(letters)

    traces.append(
        dict(
            x=np.arange(0, len(z_values[0])).tolist(),
            y=y_values,
            z=z_values,
            hovertext=text_values,
            xaxis="x",
            yaxis="y",
            type="heatmap",
            colorscale=COLORSCALE,
            showscale=False,
        )
    )

    data_json = json.dumps(traces).encode("utf-8", "ignore").decode("utf-8")
    data_json = data_json.replace("NaN", "null")

    template = TEMPLATE.replace("{{target}}", reference_name)
    template = template.replace("{{data}}", data_json)
    print(template, file=output)


if __name__ == "__main__":
    main()
