#!/usr/bin/env python

import gzip
import json
from collections import defaultdict, deque
from itertools import groupby, islice


TEMPLATE = """<!DOCTYPE html>
<html>

<head>
    <meta charset="utf-8" />
    <meta name="author" content="Joe Brown" />
    <title>idplot</title>
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>
    <script type="text/javascript" src="https://d3js.org/d3.v3.min.js"></script>

    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <link href="https://fonts.googleapis.com/css2?family=Righteous&display=swap" rel="stylesheet">

    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta1/dist/css/bootstrap.min.css" rel="stylesheet"
        integrity="sha384-giJF6kkoqNQ00vy+HMDP7azOuL0xtbfIcaT9wjKHr8RbDVddVHyTfAAsrekwKmP1" crossorigin="anonymous">
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta1/dist/js/bootstrap.bundle.min.js"
        integrity="sha384-ygbV9kiqUc6oa4msXn9868pTtWMgiQaeYH7/t7LECLbyPA2x65Kgf80OJFdroafW"
        crossorigin="anonymous"></script>


    <style type="text/css">
        .disabled_div {
            pointer-events: none;
            opacity: 0.4
        }

        .nav {
            background-color: #000000 !important;
        }

        .brand {
            font-size: 1.2rem;
            font-family: 'Righteous', cursive;
            line-height: 1.7;
        }

        .chart-row {
            overflow-x: auto;
            height: 235px;
            max-width: 100%;
        }

        .tree-view {
            width: 430px;
        }

        .tree-view>svg {
            cursor: pointer;
            pointer-events: all;
        }

        .tree-container {
            display: flex;
        }

        .dropdown-header {
            padding: .25rem .5rem !important;
            font-weight: bold;
        }

        .meta-value {
            overflow-x: scroll;
            /* white-space: nowrap !important; */
            font-size: .8rem;
            padding-right: .5rem;
        }

        *[id] {
            scroll-margin-top: 79px;
        }

        .dropdown-details {
            width: 400px;
        }

        .small {
            font-size: 80%;
            font-weight: 400;
        }

        .btn-vsm {
            vertical-align: inherit !important;
            padding: 0rem 0rem !important;
            font-size: inherit !important;
            line-height: 1 !important;
        }

        .code {
            color: #6b6b6b;
            word-break: break-word;
        }

        .form-control::placeholder {
            color: #6c757d;
            opacity: 1;
            text-overflow: ellipsis;
        }

        .plot-color {
            font-weight: 900;
            font-size: 1.2rem;
        }

        .type-select {
            width: 230px;
        }

        .btn-group-sm>.btn,
        .btn-sm {
            padding: .27rem .5rem !important;
        }
    </style>
</head>

<body>
    <div class="row bg-dark mx-0 p-1 sticky-top nav">
        <div class="col-2"><a class="brand text-white text-decoration-none"
                href="https://github.com/brwnj/idplot">idplot</a></div>
        <div class="col-10 d-flex align-items-center justify-content-end" id="meta-header">
            <div class="input-group input-group-sm type-select pe-2 d-none" id="annotation-select">
                <span class="input-group-text">Annotation</span>
                <select class="form-select" id="annotation-type">
                </select>
            </div>
            <div class="dropdown">
                <button class="btn btn-sm btn-primary dropdown-toggle" type="button" id="details"
                    data-bs-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    Run details
                </button>
                <div class="dropdown-menu dropdown-menu-right dropdown-details" aria-labelledby="details">
                    <h6 class="dropdown-header">Reference</h6>
                    <div class="container meta-value" id="meta-reference"><code></code></div>
                    <h6 class="dropdown-header">Alignment length</h6>
                    <div class="container meta-value" id="meta-length"></div>
                    <h6 class="dropdown-header">Nextflow command</h6>
                    <div class="container meta-value" id="meta-cli"></div>
                    <h6 class="dropdown-header">Launch directory</h6>
                    <div class="container meta-value" id="meta-dir"></div>
                    <h6 class="dropdown-header">Workflow container</h6>
                    <div class="container meta-value" id="meta-container"></div>
                </div>
            </div>
        </div>
    </div>
    <div class="container-fluid w-90">
        <div class="row p-2 bg-light d-none" id="dendrograms-row-wrapper">
            <div class="col-4">
                <h5>GARD refinements</h5>
            </div>
            <div class="col-8">
                <h5>GARD breakpoint trees (iteration <span id="iteration-number">00</span>)</h5>
            </div>
            <div class="col-4 ps-0" id="gard-plot"></div>
            <div class="col-8">
                <div class="row chart-row ms-0" id="dendrograms-row">
                    <div class="col-9 tree-container px-0" id="dendrograms"></div>
                </div>
            </div>
        </div>
        <div class="row">
            <div class="col-12">
                <div class="row pt-2 mb-3" id="grid-plot"></div>
            </div>
        </div>
        <div class="row p-2 d-flex">
            <div class="col-6">
                <h5>Sequences</h5>
            </div>
            <div class="col-6 d-flex align-items-center justify-content-end">
                <button class="btn btn-primary" id="export-button" type="button" data-toggle="tooltip"
                    title="Export all sequences as .fasta file">Export all</button>
            </div>
            <div class="col-6">
                <div class="text-muted small">Sequence selection is based on plot zoom level</div>
            </div>
            <div class="col-6 d-flex py-2 align-items-center justify-content-end">
                <div class="form-check form-check-inline me-0">
                    <input class="form-check-input" type="checkbox" value="" id="remove-gaps">
                    <label class="form-check-label" for="remove-gaps" title="Remove gaps (-) from sequence exports"
                        data-toggle="tooltip">
                        Remove gaps
                    </label>
                </div>
            </div>

            <div class="col-12 pt-2 px-4 mb-4" id="sequence-selection">
            </div>
        </div>
    </div>
</body>

<script>
    /*
  d3.phylogram.js

  Copyright (c) 2013, Ken-ichi Ueda

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer. Redistributions in binary
  form must reproduce the above copyright notice, this list of conditions and
  the following disclaimer in the documentation and/or other materials
  provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/

    if (!d3) { throw "d3 wasn't included!" };
    (function () {
        d3.phylogram = {}
        d3.phylogram.rightAngleDiagonal = function () {
            var projection = function (d) { return [d.y, d.x]; }

            var path = function (pathData) {
                return "M" + pathData[0] + ' ' + pathData[1] + " " + pathData[2];
            }

            function diagonal(diagonalPath, i) {
                var source = diagonalPath.source,
                    target = diagonalPath.target,
                    midpointX = (source.x + target.x) / 2,
                    midpointY = (source.y + target.y) / 2,
                    pathData = [source, { x: target.x, y: source.y }, target];
                pathData = pathData.map(projection);
                return path(pathData)
            }

            diagonal.projection = function (x) {
                if (!arguments.length) return projection;
                projection = x;
                return diagonal;
            };

            diagonal.path = function (x) {
                if (!arguments.length) return path;
                path = x;
                return diagonal;
            };

            return diagonal;
        }

        d3.phylogram.radialRightAngleDiagonal = function () {
            return d3.phylogram.rightAngleDiagonal()
                .path(function (pathData) {
                    var src = pathData[0],
                        mid = pathData[1],
                        dst = pathData[2],
                        radius = Math.sqrt(src[0] * src[0] + src[1] * src[1]),
                        srcAngle = d3.phylogram.coordinateToAngle(src, radius),
                        midAngle = d3.phylogram.coordinateToAngle(mid, radius),
                        clockwise = Math.abs(midAngle - srcAngle) > Math.PI ? midAngle <= srcAngle : midAngle > srcAngle,
                        rotation = 0,
                        largeArc = 0,
                        sweep = clockwise ? 0 : 1;
                    return 'M' + src + ' ' +
                        "A" + [radius, radius] + ' ' + rotation + ' ' + largeArc + ',' + sweep + ' ' + mid +
                        'L' + dst;
                })
                .projection(function (d) {
                    var r = d.y, a = (d.x - 90) / 180 * Math.PI;
                    return [r * Math.cos(a), r * Math.sin(a)];
                })
        }

        // Convert XY and radius to angle of a circle centered at 0,0
        d3.phylogram.coordinateToAngle = function (coord, radius) {
            var wholeAngle = 2 * Math.PI,
                quarterAngle = wholeAngle / 4

            var coordQuad = coord[0] >= 0 ? (coord[1] >= 0 ? 1 : 2) : (coord[1] >= 0 ? 4 : 3),
                coordBaseAngle = Math.abs(Math.asin(coord[1] / radius))

            // Since this is just based on the angle of the right triangle formed
            // by the coordinate and the origin, each quad will have different
            // offsets
            switch (coordQuad) {
                case 1:
                    coordAngle = quarterAngle - coordBaseAngle
                    break
                case 2:
                    coordAngle = quarterAngle + coordBaseAngle
                    break
                case 3:
                    coordAngle = 2 * quarterAngle + quarterAngle - coordBaseAngle
                    break
                case 4:
                    coordAngle = 3 * quarterAngle + coordBaseAngle
            }
            return coordAngle
        }

        d3.phylogram.styleTreeNodes = function (vis) {
            vis.selectAll('g.leaf.node')
                .append("svg:circle")
                .attr("r", 4.5)
                .attr('stroke', 'black')
                .attr('stroke-width', '1px')
                .attr('fill', function (d) { return strain_colors(d.name) || 'white' });

            vis.selectAll('g.root.node')
                .append('svg:circle')
                .attr("r", 4.5)
                .attr('fill', 'black')
                .attr('stroke', 'black')
                .attr('stroke-width', '1px');
        }

        function scaleBranchLengths(nodes, w) {
            // Visit all nodes and adjust y pos width distance metric
            var visitPreOrder = function (root, callback) {
                callback(root)
                if (root.children) {
                    for (var i = root.children.length - 1; i >= 0; i--) {
                        visitPreOrder(root.children[i], callback)
                    };
                }
            }
            visitPreOrder(nodes[0], function (node) {
                node.rootDist = (node.parent ? node.parent.rootDist : 0) + (node.length || 0)
            })
            var rootDists = nodes.map(function (n) { return n.rootDist; });
            var yscale = d3.scale.linear()
                .domain([0, d3.max(rootDists)])
                .range([0, w]);
            visitPreOrder(nodes[0], function (node) {
                node.y = yscale(node.rootDist)
            })
            return yscale
        }

        d3.phylogram.build = function (selector, nodes, options) {
            options = options || {}
            var w = options.width || d3.select(selector).style('width') || d3.select(selector).attr('width'),
                h = options.height || d3.select(selector).style('height') || d3.select(selector).attr('height'),
                w = parseInt(w),
                h = parseInt(h);
            var tree = options.tree || d3.layout.cluster()
                .size([h, w])
                .sort(function (node) { return node.children ? node.children.length : -1; })
                .children(options.children || function (node) {
                    return node.branchset
                });
            var diagonal = options.diagonal || d3.phylogram.rightAngleDiagonal();
            var vis = options.vis || d3.select(selector).append("svg:svg")
                .attr("width", w + 100)
                .attr("height", h + 10)
                //.attr("width", w + 300)
                //.attr("height", h + 30)
                .append("svg:g")
                .attr("transform", "translate(20, 20)");
            var nodes = tree(nodes);

            if (options.skipBranchLengthScaling) {
                var yscale = d3.scale.linear()
                    .domain([0, w])
                    .range([0, w]);
            } else {
                var yscale = scaleBranchLengths(nodes, w)
            }

            if (!options.skipTicks) {
                vis.selectAll('line')
                    .data(yscale.ticks(10))
                    .enter().append('svg:line')
                    .attr('y1', 0)
                    .attr('y2', h)
                    .attr('x1', yscale)
                    .attr('x2', yscale)
                    .attr("stroke", "#ddd");

                vis.selectAll("text.rule")
                    .data(yscale.ticks(10))
                    .enter().append("svg:text")
                    .attr("class", "rule")
                    .attr("x", yscale)
                    .attr("y", 0)
                    .attr("dy", -3)
                    .attr("text-anchor", "middle")
                    .attr('font-size', '8px')
                    .attr('fill', '#000')
                    .text(function (d) { return Math.round(d * 100) / 100; });
            }

            var link = vis.selectAll("path.link")
                .data(tree.links(nodes))
                .enter().append("svg:path")
                .attr("class", "link")
                .attr("d", diagonal)
                .attr("fill", "none")
                .attr("stroke", "#000")
                .attr("stroke-width", "1px");

            var node = vis.selectAll("g.node")
                .data(nodes)
                .enter().append("svg:g")
                .attr("class", function (n) {
                    if (n.children) {
                        if (n.depth == 0) {
                            return "root node"
                        } else {
                            return "inner node"
                        }
                    } else {
                        return "leaf node"
                    }
                })
                .attr("transform", function (d) { return "translate(" + d.y + "," + d.x + ")"; })

            d3.phylogram.styleTreeNodes(vis)

            if (!options.skipLabels) {
                //vis.selectAll('g.inner.node')
                //  .append("svg:text")
                //    .attr("dx", -6)
                //    .attr("dy", -6)
                //    .attr("text-anchor", 'end')
                //    .attr('fill', '#ccc')
                //    .text(function(d) { return d.length; });

                vis.selectAll('g.leaf.node').append("svg:text")
                    .attr("dx", 8)
                    .attr("dy", 3)
                    .attr("text-anchor", "start")
                    .attr('font-family', 'Helvetica Neue, Helvetica, sans-serif')
                    .attr('font-size', '10px')
                    .attr('fill', 'black')
                    .text(function (d) { return d.name; });
            }

            return { tree: tree, vis: vis }
        }

    }());

    /**
     * Newick format parser in JavaScript.
     *
     * Copyright (c) Jason Davies 2010.
     *
     * Permission is hereby granted, free of charge, to any person obtaining a copy
     * of this software and associated documentation files (the "Software"), to deal
     * in the Software without restriction, including without limitation the rights
     * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
     * copies of the Software, and to permit persons to whom the Software is
     * furnished to do so, subject to the following conditions:
     *
     * The above copyright notice and this permission notice shall be included in
     * all copies or substantial portions of the Software.
     */
    (function (exports) {
        exports.parse = function (s) {
            var ancestors = [];
            var tree = {};
            var tokens = s.split(/\\s*(;|\\(|\\)|,|:)\\s*/);
            for (var i = 0; i < tokens.length; i++) {
                var token = tokens[i];
                switch (token) {
                    case '(': // new branchset
                        var subtree = {};
                        tree.branchset = [subtree];
                        ancestors.push(tree);
                        tree = subtree;
                        break;
                    case ',': // another branch
                        var subtree = {};
                        ancestors[ancestors.length - 1].branchset.push(subtree);
                        tree = subtree;
                        break;
                    case ')': // optional name next
                        tree = ancestors.pop();
                        break;
                    case ':': // optional length next
                        break;
                    default:
                        var x = tokens[i - 1];
                        if (x == ')' || x == '(' || x == ',') {
                            tree.name = token;
                        } else if (x == ':') {
                            tree.length = parseFloat(token);
                        }
                }
            }
            return tree;
        };
    })(
        // exports will be set in any commonjs platform; use it if it's available
        typeof exports !== "undefined" ?
            exports :
            // otherwise construct a name space.  outside the anonymous function,
            // "this" will always be "window" in a browser, even in strict mode.
            this.Newick = {}
    );

    const data = {{data}}
    const cov_color = 'rgba(108,117,125,0.2)'
    const colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#7CB5EC", "#434348",
        "#90ED7D", "#F7A35C", "#8085E9", "#F15C80", "#E4D354", "#2B908F",
        "#F45B5B", "#91E8E1", "#4E79A7", "#F28E2C", "#E15759", "#76B7B2",
        "#59A14F", "#EDC949", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AB",
    ]
    const colorscale = [
        [0, "#2ca02c"], [0.2, "#2ca02c"], [0.2, "#1f77b4"], [0.4, "#1f77b4"], [0.4, "#ff7f0e"],
        [0.6, "#ff7f0e"], [0.6, "#d62728"], [0.8, "#d62728"], [0.8, "#7f7f7f"], [1, "#7f7f7f"],
    ]
    let grid_traces = []
    let selected_trees = false
    let selected_index

    jQuery('.dropdown-menu').on("click.bs.dropdown", (e) => {
        e.stopPropagation()
        e.preventDefault()
    })

    const strain_colors = (strain_id) => {
        for (const [i, strain] of Object.keys(data.queries).entries()) {
            if (strain_id == strain) {
                return colors[i % colors.length]
            }
        }
        return "#000000"
    }

    const plot_layout = () => {
        let layout = {
            title: "",
            margin: { t: 10, b: 40, r: 40 },
            height: 550 + (Object.keys(data.queries).length * 10),
            xaxis: { title: "Position", autorange: true, showgrid: false, showlines: false, zeroline: false, rangeslider: {} },
            yaxis: { title: "", fixedrange: true, showgrid: false, showspikes: false, domain: [0.65, 1], automargin: true },
            yaxis2: { title: "ANI", showgrid: true, showticklabels: true, tickmode: 'array', tickvals: [0, 0.2, 0.4, 0.6, 0.8, 1], zeroline: true, domain: [0, 0.60] },
            yaxis3: {},
            yaxis4: {},
            hovermode: "closest",
            showlegend: false,
            grid: { rows: 2, columns: 1, subplots: [['xy2', 'xy']], pattern: 'independent' },
            shapes: [],
        }
        if (data.gard) {
            layout.yaxis2.domain = [0, 0.45]
            layout.yaxis.domain = [0.60, 1]
            layout.yaxis4 = {
                title: "GARD", fixedrange: true, range: [-2, 3], showticklabels: false, showgrid: false, zeroline: false, domain: [0.45, 0.60]
            }
            layout.grid.rows += 1
            layout.grid.subplots[0].push('xy4')
        }
        return layout
    }

    const plot_config = {
        displaylogo: false,
        modeBarButtonsToRemove: ["select2d", "lasso2d"],
    }

    const get_ani_traces = (q, window) => {
        let traces = []
        for (const [strain_id, strain_data] of Object.entries(q)) {
            let trace = {
                x: Array.from({ length: strain_data.identity.length }, (v, k) => (window / 2) + k),
                y: strain_data.identity,
                text: strain_id,
                xaxis: "x",
                yaxis: "y2",
                connectgaps: false,
                hoverinfo: "text+x+y",
                type: "scatter",
                mode: "lines",
                name: "significant",
                marker: {
                    width: 1,
                    color: strain_colors(strain_id)
                }
            }
            traces.push(trace)
        }
        return traces
    }

    const get_msa_traces = (queries, reference) => {
        let z = []
        let y = []
        let text = []
        for (const [strain_id, strain_data] of Object.entries(queries)) {
            z.push(strain_data.z)
            y.push(strain_id)
            text.push([])
        }
        z.push(reference.msa)
        y.push(reference.name)
        text.push(Array.from(reference.seq))

        return {
            x: Array.from({ length: reference.seq.length }, (v, k) => k),
            y: y,
            z: z,
            text: text,
            hoverinfo: "text+x+y",
            hoverongaps: false,
            xaxis: "x",
            yaxis: "y",
            type: "heatmap",
            tracktype: "msa",
            colorscale: colorscale,
            showscale: false,
        }
    }

    const get_gard_trace = () => {
        let x = []
        let y = []
        let text = []
        if (!(selected_trees)) {
            return []
        }
        for (const [idx, arr] of selected_trees.entries()) {
            let start = arr[0]
            let end = arr[1]
            // `range` replacement
            let a = Array.from({ length: end - start + 1 }, (v, k) => k + start)
            for (let i of a) {
                x.push(i)
                y.push(idx % 2)
                text.push(`\${start}-\${end}`)
            }
            x.push("")
            y.push("")
            text.push("")
        }
        return {
            x: x,
            y: y,
            text: text,
            xaxis: "x",
            yaxis: "y4",
            hoverinfo: "text",
            type: "scatter",
            name: "breakpoints",
            tracktype: "breakpoints",
            connectgaps: false,
            showlegend: false,
            line: {
                width: 10, color: colors[7]
            },
        }
    }

    const get_gapped_sequence = (start, end, known_gaps = 0) => {
        let seq = data.reference.seq.slice(start, end + known_gaps)
        // count the selected number of gaps
        let gaps = [...seq].filter(l => l === '-').length
        // restart if we added new gaps
        if (gaps > known_gaps) {
            return get_gapped_sequence(start, end, known_gaps = gaps)
        }
        return [seq, gaps]
    }

    const get_annotation_trace = (traces) => {
        if (!(data.gff)) {
            return []
        }

        if (Object.keys(data.gff).length == 0) {
            return []
        }

        let offset = 1.1
        let start = 0
        let end = 0
        let gaps = 0
        let x = []
        let y = []
        let text = []

        select = document.getElementById("annotation-type")
        regions = data.gff[select.value]

        for (let i = 0; i < regions.length; i++) {

            // account for any gaps introduced prior to the first region
            if (i == 0) {
                let s = get_gapped_sequence(0, regions[i][0])
                gaps += s[1]
            }

            region = regions[i]
            let t = region[2].replaceAll(";", "<br>")
            offset = i % 2 == 0 ? 1.1 : 1.05

            // add existing gaps
            start = region[0] + gaps
            end = region[1] + gaps

            // grab the gapped reference sequence for this region
            let seqdata = get_gapped_sequence(start, end)
            let seq = seqdata[0]
            let g = seqdata[1]
            // update gaps tracker
            gaps += g

            // add current gaps
            end += g

            x.push(start)
            y.push(offset)
            text.push(t)
            x.push(end)
            y.push(offset)
            text.push(t)
            x.push("")
            y.push("")
            text.push("")
        }
        traces.push({
            x: x,
            y: y,
            text: text,
            xaxis: "x",
            yaxis: "y2",
            type: "scattergl",
            name: "annotation",
            connectgaps: false,
            showlegend: false,
            line: { width: 2, color: "black" },
            mode: "lines+markers",
            hoverinfo: "text+x+name",
            hoverlabel: { namelength: -1 },
            marker: {
                size: 6,
                symbol: "square",
                color: "black",
                line: { width: 1, color: "white" },
            },
        })
        return traces
    }

    const build_grid_plots = () => {
        let ani_traces = get_ani_traces(data.queries, data.window)
        let msa_trace = get_msa_traces(data.queries, data.reference)
        let gard_trace = get_gard_trace()
        let annotation_trace = get_annotation_trace(ani_traces)
        // global var
        grid_traces = [msa_trace, gard_trace, annotation_trace, ...ani_traces]

        // let layout = plot_layout()
        let grid_plot = document.getElementById("grid-plot")
        let p_obj = Plotly.react(grid_plot, grid_traces, plot_layout(), plot_config)
        grid_plot.removeAllListeners("plotly_click")
        grid_plot.removeAllListeners("plotly_doubleclick")
        grid_plot.on("plotly_click", handle_plot_click)
        grid_plot.on("plotly_doubleclick", handle_plot_doubleclick)
        grid_plot.on("plotly_relayout", draw_sequences)
        grid_plot.on("plotly_afterplot", () => {
            let yticks = jQuery("#grid-plot .yaxislayer-above > .ytick > text")
            for (const [key, value] of Object.entries(yticks)) {
                tick = jQuery(value)
                if (Object.keys(data.queries).includes(tick.text())) {
                    tick.css({
                        fill: strain_colors(tick.text()),
                        "font-weight": 600,
                    })
                }
            }
        })

        grid_plot.classList.remove("disabled_div")
        return p_obj
    }

    const update_grid_plot = () => {
        grid_traces[1] = get_gard_trace()
        let grid_plot = document.getElementById("grid-plot")
        Plotly.react(grid_plot, grid_traces, plot_layout(), plot_config)
    }

    const handle_dendrogram_click = (start, end, scroll=true) => {
        let treeplot = document.getElementById(`\${start}-\${end}-dendrogram`)
        if (treeplot == null) {
            return
        }
        jQuery(".tree-view").removeClass("border-primary bg-white")
        jQuery(".tree-view").addClass("border-light bg-light")

        treeplot.classList.add("border-primary")
        treeplot.classList.add("bg-white")
        treeplot.classList.remove("border-light")

        if (scroll) {
            treeplot.scrollIntoView({ behavior: "smooth", block: "start", inline: "center" })
        }

        // highlight this region in grid-plot
        let shape = [{
            type: "rect",
            x0: start,
            y0: 0.47,
            x1: end,
            y1: 0.57,
            xref: "x",
            yref: "paper",
            line: {
                width: 1,
                color: "#007bff",
            }
        }]
        Plotly.relayout("grid-plot", { shapes: shape })
    }

    const handle_plot_click = (click_data) => {
        if (click_data.points[0].data.tracktype == 'breakpoints') {
            let bp = click_data.points[0].text

            let coords = bp.split("-")
            let start = coords[0]
            let end = coords[1]
            // selected_coords = [start, end]
            handle_dendrogram_click(start, end)
        } else {
            let sample_id = click_data.points[0].data.text
            if (sample_id) {
                highlight_plot_traces(sample_id)
            }
        }
    }

    const handle_plot_doubleclick = () => {
        jQuery(".tree-view").removeClass("border-primary bg-white")
        jQuery(".tree-view").addClass("border-light bg-light")
        build_grid_plots()
    }

    const highlight_plot_traces = (sample_id) => {
        let s_traces = []
        let c_traces = []
        let k_traces = []
        let highlight_color;
        for (var i = 0; i < grid_traces.length; i++) {
            // let trace = grid_traces[i]
            let trace = jQuery.extend(true, {}, grid_traces[i])
            // limit to significant sample traces
            if (trace.name == "significant") {
                if (sample_id.includes(trace.text)) {
                    highlight_color = grid_traces[i].marker.color
                }
                else {
                    trace.marker.color = cov_color
                }
            }
            s_traces.push(trace)
        }
        Plotly.react("grid-plot", s_traces, plot_layout(), plot_config)
    }

    const build_newick = (str, div_id) => {
        let newick = Newick.parse(str)
        let nodes = []
        function build_newick_nodes(node, callback) {
            nodes.push(node)
            if (node.branchset) {
                for (let i; i < node.branchset.length; i++) {
                    build_newick_nodes(node.branchset[i])
                }
            }
        }
        build_newick_nodes(newick)
        let d = document.getElementById(div_id)
        d3.phylogram.build(d, newick, {
            width: 330,
            height: 180
        })
        d.addEventListener("mouseover", (e) => {
            let coords = e.target.parentElement.id.split("-")
            handle_dendrogram_click(coords[0], coords[1], scroll=false)
        })
    }

    const zoom_grid_plot = (id) => {
        let coords = id.split('-')
        let start = coords[0]
        let end = coords[1]

        let layout = plot_layout()
        layout.xaxis.range = [start, end]
        layout.xaxis.autorange = false
        Plotly.relayout("grid-plot", layout)
    }

    const build_dendrograms = () => {
        return new Promise((resolve, reject) => {
            document.getElementById("dendrograms-row-wrapper").classList.remove("d-none")
            let d = document.getElementById("dendrograms")
            d.innerHTML = ""
            for (const arr of selected_trees) {
                let start = arr[0]
                let end = arr[1]
                let newick = arr[2]
                let field_id = `\${start}-\${end}`
                d.insertAdjacentHTML("beforeend", `
                    <div class="tree-view border border-light rounded pt-1 pl-1" id="\${field_id}-dendrogram">
                        <div class="small">Region: <button type="button" class="btn btn-vsm btn-link tree-zoom" id="\${field_id}-zoom" onClick="zoom_grid_plot(this.id)">\${field_id}</div>
                    </div>
                `)
                build_newick(newick, `\${field_id}-dendrogram`)
            }
            resolve()
        })
    }

    const update_details = (name, length, cli, dir, container) => {
        document.getElementById("meta-reference").innerHTML = name
        document.getElementById("meta-length").innerHTML = length
        document.getElementById("meta-cli").innerHTML = `<code>\${cli}<code>`
        document.getElementById("meta-dir").innerHTML = `<code>\${dir}<code>`
        document.getElementById("meta-container").innerHTML = `<code>\${container}<code>`
    }

    const get_selected_range = () => {
        let p = document.getElementById("grid-plot")
        let start = p.layout.xaxis.range[0] < 0 ? 0 : Math.round(p.layout.xaxis.range[0])
        let end = Math.round(p.layout.xaxis.range[1])
        return [start, end]
    }

    const get_seq_by_id = (id) => {
        let [start, end] = get_selected_range()
        let seq
        if (id == data.reference.name) {
            seq = data.reference.seq.substring(start, end)
        } else {
            seq = data.queries[id].seq.substring(start, end)
        }
        if (jQuery("#remove-gaps").is(":checked")) {
            seq = seq.replaceAll("-", "")
        }
        return seq
    }

    const build_gard_plot = () => {
        let x = []
        let y = []
        let selected_x
        let selected_y
        cum_daic = 0
        for (const [imp_id, bp_obj] of Object.entries(data.gard.improvements)) {
            cum_daic += bp_obj.deltaAICc
            x.push(parseInt(imp_id))
            y.push(cum_daic)
            selected_x = parseInt(imp_id)
            selected_y = cum_daic
        }

        let layout = {
            margin: { t: 0, b: 30, r: 20 },
            height: 235,
            paper_bgcolor: '#f8f9fa',
            xaxis: {
                title: "Iteration",
                rangemode: "tozero",
            },
            yaxis: {
                title: "Cumulative &#x394;AIC",
            },
            annotations: [{
                text: "Selected: " + selected_x,
                x: selected_x,
                y: selected_y,
                showarrow: true,
                arrowhead: 6,
                ax: 0,
                ay: 40,
                font: {size: 14},
                bgcolor: 'rgba(255, 255, 255, 0.8)',
            }]
        }

        trace = [{
            x: x,
            y: y,
            mode: "lines+markers",
            marker: { size: 10, color: "rgb(211,211,211,0.5)", line: { color: 'rgb(40,40,40)', width: 1 } },
            line: {
                color: "black",
                width: 1
            }
        }]

        let p = document.getElementById("gard-plot")
        let p_obj = Plotly.react(p, trace, layout, { displayModeBar: false })
        p.removeAllListeners("plotly_click")
        p.on("plotly_click", handle_iteration_click)
        return p_obj
    }

    const get_trees = (idx) => {
        let start = 0
        let trees = []
        for (let end of data.gard.improvements[idx].breakpoints) {
            end = end[0]
            trees.push([start, end, data.trees[`\${start}_\${end}`].replace(";", "")])
            start = end + 1
        }
        trees.push([start, data.reference.seq.length, data.trees[`\${start}_\${data.reference.seq.length}`]])
        return trees
    }

    const update_tree_data = () => {
        build_dendrograms().then(() => {
            handle_dendrogram_click(0, data.gard.improvements[selected_index].breakpoints[0][0])
        })
        update_grid_plot()
    }

    const handle_iteration_click = (click_data) => {
        document.getElementById("iteration-number").innerHTML = click_data.points[0].x
        selected_index = click_data.points[0].x
        annotations = [{
            text: 'Selected: ' + click_data.points[0].x,
            x: click_data.points[0].x,
            y: parseFloat(click_data.points[0].y.toPrecision(4)),
            showarrow: true,
            arrowhead: 6,
            ax: 0,
            ay: 40,
            font: {size: 14},
            bgcolor: 'rgba(255, 255, 255, 0.8)',
        }]
        Plotly.relayout("gard-plot", { annotations: annotations })

        selected_trees = get_trees(click_data.points[0].x)
        update_tree_data()
    }

    const init_sequences = () => {
        let seq_sel = document.getElementById("sequence-selection")
        seq_sel.insertAdjacentHTML('beforeend', `
            <label for="\${data.reference.name}-seq">\${data.reference.name} (Reference)</label>
            <div class="input-group input-group-sm pb-2">
                <input class="form-control text-monospace" type="text" placeholder="\${get_seq_by_id(data.reference.name)}" id="\${data.reference.name}-seq" readonly="">
                <button class="btn btn-primary" type="button" id="\${data.reference.name}-copy-btn" title="Copy selected region" data-toggle="tooltip" onclick="copy('\${data.reference.name}')">Copy</button>
                <button class="btn btn-primary blast-btn" type="button" id="\${data.reference.name}-blast-btn" title="Send selected region to BLAST" data-toggle="tooltip" onclick="blast('\${data.reference.name}')">BLAST</button>
            </div>
            `
        )
        for (const query in data.queries) {
            seq_sel.insertAdjacentHTML('beforeend', `
                <label for="\${query}-seq"><span class="plot-color" style="color:\${strain_colors(query)}">|</span> \${query}</label>
                <div class="input-group input-group-sm pb-2">
                    <input class="form-control text-monospace" type="text" placeholder="\${get_seq_by_id(query)}" id="\${query}-seq" readonly="">
                    <button class="btn btn-primary" type="button" id="\${query}-copy-btn" title="Copy selected region" data-toggle="tooltip" onclick="copy('\${query}')">Copy</button>
                    <button class="btn btn-primary blast-btn" type="button" id="\${query}-blast-btn" title="Send selected region to BLAST" data-toggle="tooltip" onclick="blast('\${query}')">BLAST</button>
                </div>
                `
            )
        }
    }

    const toggle_blast_button = () => {
        // disable for anything longer than 8k...
        let [start, end] = get_selected_range()
        if (end - start > 8000) {
            document.querySelectorAll('.blast-btn').forEach(elem => {
                elem.disabled = true
                jQuery(`#\${elem.id}`).attr("data-original-title", "Selection too long (>8kb)").tooltip("update")
            })
        } else {
            document.querySelectorAll('.blast-btn').forEach(elem => {
                elem.disabled = false
                jQuery(`#\${elem.id}`).attr("data-original-title", "Send selected region to BLAST").tooltip("update")
            })
        }
    }

    const draw_sequences = () => {
        // reference
        document.getElementById(`\${data.reference.name}-seq`).placeholder = get_seq_by_id(data.reference.name)
        // queries
        for (const query in data.queries) {
            document.getElementById(`\${query}-seq`).placeholder = get_seq_by_id(query)
        }
        toggle_blast_button()
    }

    const copy = (id) => {
        var \$temp = \$("<input>")
        \$("body").append(\$temp)

        \$temp.val(get_seq_by_id(id)).select()
        document.execCommand("copy")
        \$temp.remove()

        \$(`#\${id}-copy-btn`).attr("title", "Copied!").tooltip("_fixTitle").tooltip("show").attr("title", "Copy selected region").tooltip("_fixTitle")
    }

    const blast = (id) => {
        let seq = get_seq_by_id(id)
        let url = new URL(`https://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY=\${seq}&DATABASE=nt&PROGRAM=blastn&CMD=put`)
        window.open(url)
    }

    // https://robkendal.co.uk/blog/2020-04-17-saving-text-to-client-side-file-using-vanilla-js
    const download_to_file = (content, filename, contentType = 'text/plain') => {
        const a = document.createElement('a');
        const file = new Blob([content], { type: contentType });

        a.href = URL.createObjectURL(file);
        a.download = filename;
        a.click();

        URL.revokeObjectURL(a.href);
    };

    jQuery("#export-button").on("click", () => {
        let fa = [">" + data.reference.name, get_seq_by_id(data.reference.name)]
        for (const query in data.queries) {
            fa.push(">" + query)
            fa.push(get_seq_by_id(query))
        }
        let [start, end] = get_selected_range()
        let filename = `\${data.reference.name}_\${start}_\${end}.fasta`
        download_to_file(fa.join("\\\\n") + "\\\\n", filename)
    })

    jQuery("#remove-gaps").on("change", () => {
        draw_sequences()
    })

    jQuery("#annotation-type").on("change", () => {
        build_grid_plots()
    })

    jQuery(document).ready(function () {
        if (data.gard) {
            document.getElementById("iteration-number").innerHTML = Object.keys(data.gard.improvements).length - 1
            selected_trees = get_trees(Object.keys(data.gard.improvements).length - 1)
        }
        if (data.gff) {
            // show the select
            document.getElementById("annotation-select").classList.remove("d-none")
            let select = document.getElementById("annotation-type")
            // populate the select
            for (const atype of Object.keys(data.gff)) {
                let opt = document.createElement("option")
                opt.value = atype
                opt.innerHTML = atype
                if (atype == "gene") {
                    opt.selected = true
                }
                select.appendChild(opt)
            }
        }
        update_details(data.reference.name, data.reference.seq.length, data.meta.cli, data.meta.dir, data.meta.container)
        build_grid_plots().then(() => {
            if (data.gard) {
                build_dendrograms().then(() => {
                    build_gard_plot().then(() => {
                        handle_dendrogram_click(0, data.gard.improvements[Object.keys(data.gard.improvements).length - 1].breakpoints[0][0], scroll = false)
                    })
                })
            }
        })
        init_sequences()
        jQuery('[data-toggle="tooltip"]').tooltip()
        toggle_blast_button()
    })
</script>
</html>
"""

# A - green, C - blue, G - yellow, T - red
nuc_map = {
    "a": 0,
    "c": 1,
    "g": 2,
    "t": 3,
    "u": 3,
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "-": 4,
}


def gzopen(f):
    if f.endswith(".gz"):
        return gzip.open(f, "rt")
    else:
        return open(f)


def read_fasta(fh):
    for header, group in groupby(fh, lambda line: line[0] == ">"):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = "".join(line.strip() for line in group)
            yield name, seq


def parse_alignments(fasta):
    reference = False
    queries = []
    with open(alignments) as fh:
        for (name, seq) in read_fasta(fh):
            if reference:
                queries.append(dict(name=name, seq=seq))
            else:
                reference = dict(name=name, seq=seq)
                msa = list()
                for b in seq:
                    try:
                        msa.append(nuc_map[b])
                    except KeyError:
                        msa.append("")
                reference["msa"] = msa
    return reference, queries


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


def process_queries(refseq, query_seqs):
    query_vals = dict()
    for query in query_seqs:
        mismatches = list()
        query_msa = list()
        # mismatches, msa of the query
        for rbase, qbase in zip(refseq, query["seq"]):
            if rbase == qbase:
                mismatches.append(0)
                query_msa.append("")
            else:
                mismatches.append(1)
                try:
                    query_msa.append(nuc_map[qbase])
                except KeyError:
                    # handle ambiguous bases
                    query_msa.append("")

        assert len(mismatches) == len(reference["seq"])

        identities = list()
        for mismatch_window in sliding_window(mismatches, window):
            mw = list(mismatch_window)
            pid = max((window - sum(mw)) / window, 0)
            identities.append(pid)

        query_vals[query["name"]] = dict(identity=identities, z=query_msa, seq=query["seq"])
    return query_vals


def parse_gard(filepath):
    if not filepath:
        return False

    with open(filepath) as fh:
        jdata = json.load(fh)
        # small amount of bs
        jdata.pop('analysis', None)
        # large amount of bs
        jdata.pop('siteBreakPointSupport', None)
        return jdata


def parse_trees(filepaths):
    if not filepaths:
        return False

    t = dict()
    for f in filepaths.split(" "):
        _, start, end = f.rpartition(".")[0].rsplit("_", 2)
        with open(f) as fh:
            newick = fh.readline().strip()
            t[f"{start}_{end}"] = newick
    return t


def parse_gff(filepath):
    if not filepath:
        return False
    gff_data = defaultdict(list)
    with gzopen(filepath) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            toks = line.strip().split("\\t")
            # make sure we ignore fasta lines when present
            if len(toks) < 9:
                continue
            gff_data[toks[2]].append([int(toks[3]), int(toks[4]), toks[8]])
    return gff_data


alignments = "$msa"
json_input = "$json" if "$json" != "input.2" else False
trees = "$trees" if "$trees" != "input.3" else False
gff = "$gff" if "$gff" != "input.4" else False
window = $params.window
output = "idplot.html"
nextflow_command = "$workflow.commandLine"
launch_directory = "$workflow.launchDir"
workflow_container = "$workflow.container"

reference, queries = parse_alignments(alignments)
queries = process_queries(reference["seq"], queries)
gard_results = parse_gard(json_input)
tree_results = parse_trees(trees)
gff_results = parse_gff(gff)
data = {
        "reference": reference,
        "queries": queries,
        "gard": gard_results,
        "trees": tree_results,
        "gff": gff_results,
        "window": window,
        "meta": {
            "cli": nextflow_command,
            "dir": launch_directory,
            "container": workflow_container,
        },
    }

with open(output, "w") as fh:
    data_json = json.dumps(data).encode("utf-8", "ignore").decode("utf-8")
    data_json = data_json.replace("NaN", "null")
    template = TEMPLATE.replace("{{data}}", data_json)
    print(template, file=fh)
