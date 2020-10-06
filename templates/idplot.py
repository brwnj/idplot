#!/usr/bin/env python

import json
from collections import defaultdict, deque
from itertools import groupby, islice


TEMPLATE = """
<!DOCTYPE html>
<html>

<head>
    <meta charset="utf-8" />
    <meta name="author" content="Joe Brown" />
    <title>idplot</title>
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>
    <script type="text/javascript" src="https://d3js.org/d3.v3.min.js"></script>
    <script type="text/javascript"
        src="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.bundle.min.js"></script>
    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <link href="https://fonts.googleapis.com/css2?family=Righteous&display=swap" rel="stylesheet">
    <link rel="stylesheet" type="text/css"
        href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">

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
            height: 220px;
            max-width: 100%;
        }

        .tree-view {
            width: 430px;
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

        .code {
            color: #6b6b6b;
            word-break: break-word;
        }

        .highlight-select {
            background-color: #ccc;
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
    </style>
</head>

<body>
    <div class="row bg-dark mx-0 p-1 sticky-top nav">
        <div class="col-2"><a class="brand text-white text-decoration-none"
                href="https://github.com/brwnj/idplot">idplot</a></div>
        <div class="col-10 d-flex align-items-center justify-content-end" id="meta-header">
            <div class="dropdown">
                <button class="btn btn-sm btn-primary dropdown-toggle" type="button" id="details"
                    data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
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
        <div class="row p-2 d-none" id="dendrograms-row-wrapper">
            <div class="col-12">
                <h5>Breakpoints</h5>
            </div>
            <div class="row chart-row" id="dendrograms-row">
                <div class="col-12 tree-container" id="dendrograms"></div>
            </div>
        </div>
        <div class="row">
            <div class="col-12">
                <div class="row pt-2 mb-3" id="grid_plot"></div>
            </div>
        </div>
        <div class="row p-2">
            <div class="col-12">
                <h5>Sequences</h5>
                <div class="text-muted small">Sequence selection is based on plot zoom level</div>
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

    \$('.dropdown-menu').on("click.bs.dropdown", (e) => {
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
        let subplots = []
        let layout = {
            title: "",
            margin: { t: 10, b: 40, r: 40 },
            height: 550,
            xaxis: { title: "Position", autorange: true, showgrid: false, showlines: false, zeroline: false, rangeslider: {} },
            yaxis: { title: "", fixedrange: true, showgrid: false, showspikes: false, domain: [0.65, 1] },
            yaxis2: { title: "ANI", fixedrange: true, range: [0, 1], showgrid: true, showticklabels: true, tickmode: 'array', tick0: 0, dtick: 0.2, zeroline: true, domain: [0, 0.60] },
            yaxis3: {},
            yaxis4: {},
            hovermode: "closest",
            showlegend: false,
            grid: { rows: 2, columns: 1, subplots: [['xy2', 'xy']], pattern: 'independent' },
        }
        if (data.gard) {
            layout.yaxis2.domain = [0, 0.40]
            layout.yaxis.domain = [0.75, 1]
            layout.yaxis4 = {
                title: "GARD", fixedrange: true, range: [-2, 3], showticklabels: false, showgrid: false, zeroline: false, domain: [0.45, 0.70]
            }
            layout.grid.rows += 1
            layout.grid.subplots[0].push('xy4')
        }
        return layout
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
        text.push(reference.seq)

        return {
            x: Array.from({ length: reference.seq.length }, (v, k) => k),
            y: y,
            z: z,
            hovertext: text,
            xaxis: "x",
            yaxis: "y",
            type: "heatmap",
            tracktype: "msa",
            colorscale: colorscale,
            showscale: false,
        }
    }

    const get_gard_trace = (breakpoints) => {
        if (!breakpoints) {
            return []
        }
        let x = []
        let y = []
        let text = []

        for (const breakpoint in breakpoints) {
            let idx = breakpoint
            let coords = breakpoints[breakpoint].bps[0]
            let start = coords[0]
            let end = coords[1]
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
            type: "scattergl",
            name: "breakpoints",
            tracktype: "breakpoints",
            connectgaps: false,
            showlegend: false,
            line: {
                width: 10, color: colors[7]
            },
        }
    }

    const build_plots = () => {
        let ani_traces = get_ani_traces(data.queries, data.window)
        let msa_trace = get_msa_traces(data.queries, data.reference)
        let gard_trace = get_gard_trace(data.gard.breakpointData)
        // global var
        grid_traces = [msa_trace, gard_trace, ...ani_traces]

        let grid_plot = document.getElementById("grid_plot")
        Plotly.react(grid_plot, grid_traces, plot_layout())
        grid_plot.removeAllListeners("plotly_click")
        grid_plot.removeAllListeners("plotly_doubleclick")
        grid_plot.on("plotly_click", handle_plot_click)
        grid_plot.on("plotly_doubleclick", handle_plot_doubleclick)
        grid_plot.on("plotly_relayout", draw_sequences)

        grid_plot.classList.remove("disabled_div")
    }

    const handle_plot_click = (click_data) => {
        if (click_data.points[0].data.tracktype == 'breakpoints') {
            let bp = click_data.points[0].text
            let treeplot = document.getElementById(`\${bp}-dendrogram`)
            treeplot.scrollIntoView({ behavior: "smooth", block: "start", inline: "center" })
            // change the background color
            treeplot.classList.add("highlight-select");
            // sleep; change the background back to white
            setTimeout(() => { treeplot.classList.remove("highlight-select") }, 800)
        } else {
            let sample_id = click_data.points[0].data.text
            if (sample_id) {
                highlight_plot_traces(sample_id)
            }
        }
    }

    const handle_plot_doubleclick = () => {
        build_plots()
    }

    const highlight_plot_traces = (sample_id) => {
        let s_traces = []
        let c_traces = []
        let k_traces = []
        let highlight_color;
        for (var i = 0; i < grid_traces.length; i++) {
            // let trace = grid_traces[i]
            let trace = \$.extend(true, {}, grid_traces[i])
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
        Plotly.react("grid_plot", s_traces, plot_layout())
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
    }

    const build_dendrograms = (tree_data) => {
        if (tree_data) {
            document.getElementById("dendrograms-row-wrapper").classList.remove("d-none")
            let d = document.getElementById("dendrograms")
            for (const idx in tree_data) {
                let tree_obj = tree_data[idx]
                let field_id = `\${tree_obj.bps[0][0]}-\${tree_obj.bps[0][1]}`
                let newick_str = tree_obj.tree
                d.insertAdjacentHTML("beforeend", `
                <div class="tree-view" id="\${field_id}-dendrogram">
                    <div class="small pl-3">Region: \${field_id}</div>
                </div>
                `)
                build_newick(newick_str, `\${field_id}-dendrogram`)
            }
        }
    }

    const update_details = (name, length, cli, dir, container) => {
        document.getElementById("meta-reference").innerHTML = name
        document.getElementById("meta-length").innerHTML = length
        document.getElementById("meta-cli").innerHTML = `<code>\${cli}<code>`
        document.getElementById("meta-dir").innerHTML = `<code>\${dir}<code>`
        document.getElementById("meta-container").innerHTML = `<code>\${container}<code>`
    }

    const get_selected_range = () => {
        let p = document.getElementById("grid_plot")
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
        return seq
    }

    const init_sequences = () => {
        let seq_sel = document.getElementById("sequence-selection")
        seq_sel.insertAdjacentHTML('beforeend', `
            <label for="\${data.reference.name}-seq">\${data.reference.name} (Reference)</label>
            <div class="input-group input-group-sm pb-2">
                <input class="form-control" type="text" placeholder="\${get_seq_by_id(data.reference.name)}" id="\${data.reference.name}-seq" readonly="">
                <div class="input-group-append">
                    <button class="btn btn-primary" type="button" id="\${data.reference.name}-copy-btn" title="Copy selected region" data-toggle="tooltip" onclick="copy('\${data.reference.name}')">Copy</button>
                    <button class="btn btn-primary blast-btn" type="button" id="\${data.reference.name}-blast-btn" title="Send selected region to BLAST" data-toggle="tooltip" onclick="blast('\${data.reference.name}')">BLAST</button>
                </div>
            </div>
            `
        )
        for (const query in data.queries) {
            seq_sel.insertAdjacentHTML('beforeend', `
                <label for="\${query}-seq"><span class="plot-color" style="color:\${strain_colors(query)}">|</span> \${query}</label>
                <div class="input-group input-group-sm pb-2">
                    <input class="form-control" type="text" placeholder="\${get_seq_by_id(query)}" id="\${query}-seq" readonly="">
                    <div class="input-group-append">
                        <button class="btn btn-primary" type="button" id="\${query}-copy-btn" title="Copy selected region" data-toggle="tooltip" onclick="copy('\${query}')">Copy</button>
                        <button class="btn btn-primary blast-btn" type="button" id="\${query}-blast-btn" title="Send selected region to BLAST" data-toggle="tooltip" onclick="blast('\${query}')">BLAST</button>
                    </div>
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
                \$(`#\${elem.id}`).attr("data-original-title", "Selection too long (>8kb)").tooltip("update")
            })
        } else {
            document.querySelectorAll('.blast-btn').forEach(elem => {
                elem.disabled = false
                \$(`#\${elem.id}`).attr("data-original-title", "Send selected region to BLAST").tooltip("update")
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

    jQuery(document).ready(function () {
        update_details(data.reference.name, data.reference.seq.length, data.meta.cli, data.meta.dir, data.meta.container)
        build_plots(data)
        build_dendrograms(data.gard.breakpointData)
        init_sequences()
        \$('[data-toggle="tooltip"]').tooltip()
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

    # points = list()
    with open(filepath) as fh:
        jdata = json.load(fh)
    #     for pos, support in jdata["breakpointData"].items():
    #         points.append([int(pos), support])
    # points.sort(key=lambda x: x[0])
        return jdata

alignments = "$msa"
json_input = "$json" if "$json" != "input.2" else False
window = $params.window
output = "idplot.html"
nextflow_command = "$workflow.commandLine"
launch_directory = "$workflow.launchDir"
workflow_container = "$workflow.container"

reference, queries = parse_alignments(alignments)
queries = process_queries(reference["seq"], queries)
gard_results = parse_gard(json_input)
data = {
        "reference": reference,
        "queries": queries,
        "gard": gard_results,
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
