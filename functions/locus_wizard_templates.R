# Helper templates for export_aligned_to_html_locus_wizard()
#
# Keeps the main exporter readable by isolating large CSS/JS blobs.

locus_wizard_build_js <- function(download_name_js,
                                 loci_js,
                                 seg_js,
                                 orig_js,
                                 colors_js,
                                 pvals_js) {
  paste0(
    "(function(){\n",
    "  // Data\n",
    "  var DOWNLOAD_NAME = ", download_name_js, ";\n",
    "  var LOCI = ", loci_js, ";\n",
    "  var SEGS = ", seg_js, ";\n",
    "  var ORIGS = ", orig_js, ";\n",
    "  var COLORS = ", colors_js, ";\n",
    "  var PVALS = ", pvals_js, ";\n",
    "  var idx = 0;\n",
    "  var dirty = false;\n",
                                "  function idEq(a, b) { return String(a) === String(b); }\n",
                                "  function getLocusRows(mergedId) {\n",
                                "    var out = [];\n",
                                "    for (var i = 0; i < SEGS.length; i++) {\n",
                                "      if (idEq(SEGS[i].MERGED_LOCUS_ID, mergedId)) out.push(SEGS[i]);\n",
                                "    }\n",
                                "    return out;\n",
                                "  }\n",
                                "  function getOrigRows(mergedId) {\n",
                                "    var out = [];\n",
                                "    for (var i = 0; i < ORIGS.length; i++) {\n",
                                "      if (idEq(ORIGS[i].MERGED_LOCUS_ID, mergedId)) out.push(ORIGS[i]);\n",
                                "    }\n",
                                "    return out;\n",
                                "  }\n",
                                "  function getPvalRows(mergedId) {\n",
                                "    var out = [];\n",
                                "    for (var i = 0; i < PVALS.length; i++) {\n",
                                "      if (idEq(PVALS[i].MERGED_LOCUS_ID, mergedId)) out.push(PVALS[i]);\n",
                                "    }\n",
                                "    return out;\n",
                                "  }\n",
                                "\n",
                                "  function updateDirty() {\n",
                                "    var el = byId('dirty');\n",
                                "    if (!el) return;\n",
                                "    if (dirty) {\n",
                                "      el.textContent = 'Unsaved';\n",
                                "      el.className = 'badge badge-warn';\n",
                                "    } else {\n",
                                "      el.textContent = 'Saved';\n",
                                "      el.className = 'badge badge-ok';\n",
                                "    }\n",
                                "  }\n",
                                "\n",
                                "  function renderLocusInfo(locus) {\n",
                                "    var elChr = byId('chr');\n",
                                "    if (elChr) elChr.textContent = locus.CHR || '';\n",
                                "    var elName = byId('locus_name_hdr');\n",
                                "    if (elName) {\n",
                                "      var nm = (locus.LOCUS_NAME || '').trim();\n",
                                "      if (!nm) nm = 'Locus ' + (locus.MERGED_LOCUS_ID || '');\n",
                                "      elName.textContent = nm;\n",
                                "    }\n",
                                "    var elCounter = byId('counter');\n",
                                "    if (elCounter) elCounter.textContent = (idx + 1) + ' / ' + LOCI.length;\n",
                                "    updateDirty();\n",
                                "  }\n",
                                "\n",
                                "  function renderSplits(arr) {\n",
                                "    var ul = byId('splits');\n",
                                "    if (!ul) return;\n",
                                "    while (ul.firstChild) ul.removeChild(ul.firstChild);\n",
                                "    if (!arr || arr.length === 0) {\n",
                                "      var li0 = document.createElement('li');\n",
                                "      li0.textContent = '(none)';\n",
                                "      li0.className = 'small';\n",
                                "      ul.appendChild(li0);\n",
                                "      return;\n",
                                "    }\n",
                                "    for (var i = 0; i < arr.length; i++) {\n",
                                "      (function(i) {\n",
                                "        var li = document.createElement('li');\n",
                                "        li.textContent = arr[i].start + '-' + arr[i].end + ' ';\n",
                                "        var btn = document.createElement('button');\n",
                                "        btn.textContent = 'x';\n",
                                "        btn.className = 'btn-small';\n",
                                "        btn.addEventListener('click', function() {\n",
                                "          var locus = LOCI[idx];\n",
                                "          var cur = parseSplits(locus.SPLITS);\n",
                                "          cur.splice(i, 1);\n",
                                "          locus.SPLITS = splitsToString(cur);\n",
                                "          dirty = true;\n",
                                "          updateDirty();\n",
                                "          render();\n",
                                "        });\n",
                                "        li.appendChild(btn);\n",
                                "        ul.appendChild(li);\n",
                                "      })(i);\n",
                                "    }\n",
                                "  }\n",
                                "\n",
                                "  function updateFromInputs() {\n",
                                "    var locus = LOCI[idx];\n",
                                "    if (!locus) return;\n",
                                "    var rs = byId('refined_start').value;\n",
                                "    var re = byId('refined_end').value;\n",
                                "    locus.REFINED_START = (rs === '' || rs === null) ? null : parseInt(rs, 10);\n",
                                "    locus.REFINED_END = (re === '' || re === null) ? null : parseInt(re, 10);\n",
                                "    locus.LOCUS_NAME = byId('locus_name').value || '';\n",
                                "    var elName = byId('locus_name_hdr');\n",
                                "    if (elName) {\n",
                                "      var nm = (locus.LOCUS_NAME || '').trim();\n",
                                "      if (!nm) nm = 'Locus ' + (locus.MERGED_LOCUS_ID || '');\n",
                                "      elName.textContent = nm;\n",
                                "    }\n",
                                "    locus.KEEP = !!byId('keep').checked;\n",
                                "    dirty = true;\n",
                                "    updateDirty();\n",
                                "  }\n",
                                "\n",
                                "  function validateInputs() {\n",
                                "    var locus = LOCI[idx];\n",
                                "    if (!locus) return;\n",
                                "    var a = locus.REFINED_START;\n",
                                "    var b = locus.REFINED_END;\n",
                                "    var ok = isFinite(a) && isFinite(b) && a <= b;\n",
                                "    var elA = byId('refined_start');\n",
                                "    var elB = byId('refined_end');\n",
                                "    if (elA) elA.classList.toggle('invalid', !ok);\n",
                                "    if (elB) elB.classList.toggle('invalid', !ok);\n",
                                "  }\n",
                                "\n",
    "  var selection = null;\n",
    "\n",
    "  function byId(id) { return document.getElementById(id); }\n",
    "\n",
    "  function showError(msg) {\n",
    "    var el = byId('js_error');\n",
    "    if (!el) return;\n",
    "    el.style.display = 'block';\n",
    "    el.textContent = String(msg);\n",
    "  }\n",
    "\n",
    "  function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }\n",
    "\n",
    "  function escTSV(x) {\n",
    "    if (x === null || x === undefined) return '';\n",
    "    var s = String(x);\n",
    "    return s.replace(/\\t/g, ' ').replace(/\\r?\\n/g, ' ');\n",
    "  }\n",
    "\n",
    "  function fmtInt(x) {\n",
    "    if (x === null || x === undefined || !isFinite(x)) return '';\n",
    "    var s = String(Math.round(x));\n",
    "    return s.replace(/\\B(?=(\\d{3})+(?!\\d))/g, ',');\n",
    "  }\n",
    "\n",
    "  function parseSplits(s) {\n",
    "    if (!s) return [];\n",
    "    if (Array.isArray(s)) return s;\n",
    "    var str = String(s).trim();\n",
    "    if (!str) return [];\n",
    "    var parts = str.split(/\\s*;\\s*/);\n",
    "    var out = [];\n",
    "    for (var i = 0; i < parts.length; i++) {\n",
    "      var p = parts[i];\n",
    "      if (!p) continue;\n",
    "      var m = p.split(/\\s*-\\s*/);\n",
    "      if (m.length < 2) continue;\n",
    "      var a = parseInt(m[0], 10);\n",
    "      var b = parseInt(m[1], 10);\n",
    "      if (!isFinite(a) || !isFinite(b)) continue;\n",
    "      out.push({ start: Math.min(a, b), end: Math.max(a, b) });\n",
    "    }\n",
    "    return out;\n",
    "  }\n",
    "\n",
    "  function splitsToString(arr) {\n",
    "    if (!arr || arr.length === 0) return '';\n",
    "    var parts = [];\n",
    "    for (var i = 0; i < arr.length; i++) {\n",
    "      if (!arr[i]) continue;\n",
    "      var a = arr[i].start;\n",
    "      var b = arr[i].end;\n",
    "      if (!isFinite(a) || !isFinite(b)) continue;\n",
    "      parts.push(String(Math.round(a)) + '-' + String(Math.round(b)));\n",
    "    }\n",
    "    return parts.join(';');\n",
    "  }\n",
    "\n",
    "  function renderForm(locus) {\n",
    "    byId('refined_start').value = isFinite(locus.REFINED_START) ? locus.REFINED_START : '';\n",
    "    byId('refined_end').value = isFinite(locus.REFINED_END) ? locus.REFINED_END : '';\n",
    "    byId('locus_name').value = locus.LOCUS_NAME || '';\n",
    "    byId('keep').checked = !!locus.KEEP;\n",
    "  }\n",
    "\n",
    "  function clearSvg(svg) { while (svg.firstChild) svg.removeChild(svg.firstChild); }\n",
    "\n",
    "  function svgEl(name) { return document.createElementNS('http://www.w3.org/2000/svg', name); }\n",
    "\n",
    "  function xScale(pos, start, end, width, padL, padR) {\n",
    "    var inner = width - padL - padR;\n",
    "    if (end <= start) return padL;\n",
    "    var frac = (pos - start) / (end - start);\n",
    "    frac = clamp(frac, 0, 1);\n",
    "    return padL + frac * inner;\n",
    "  }\n",
    "\n",
    "  function drawRuler(svg, start, end, width) {\n",
    "    var y = 35;\n",
    "    var line = svgEl('line');\n",
    "    line.setAttribute('x1', 60);\n",
    "    line.setAttribute('x2', width - 20);\n",
    "    line.setAttribute('y1', y);\n",
    "    line.setAttribute('y2', y);\n",
    "    line.setAttribute('class', 'ruler-line');\n",
    "    svg.appendChild(line);\n",
    "\n",
    "    var ticks = 5;\n",
    "    for (var i = 0; i <= ticks; i++) {\n",
    "      var p = start + (i / ticks) * (end - start);\n",
    "      var x = xScale(p, start, end, width, 60, 20);\n",
    "      var t = svgEl('line');\n",
    "      t.setAttribute('x1', x);\n",
    "      t.setAttribute('x2', x);\n",
    "      t.setAttribute('y1', y - 4);\n",
    "      t.setAttribute('y2', y + 4);\n",
    "      t.setAttribute('stroke', '#4a5568');\n",
    "      t.setAttribute('stroke-width', '1');\n",
    "      svg.appendChild(t);\n",
    "\n",
    "      var lbl = svgEl('text');\n",
    "      lbl.setAttribute('x', x);\n",
    "      lbl.setAttribute('y', y - 8);\n",
    "      lbl.setAttribute('text-anchor', 'middle');\n",
    "      lbl.setAttribute('class', 'ruler-label');\n",
    "      lbl.textContent = fmtInt(p);\n",
    "      svg.appendChild(lbl);\n",
    "    }\n",
    "  }\n",
    "\n",
    "  function drawInterval(svg, start, end, segStart, segEnd, y, height, color, label, width) {\n",
    "    var x1 = xScale(segStart, start, end, width, 60, 20);\n",
    "    var x2 = xScale(segEnd, start, end, width, 60, 20);\n",
    "    var rect = svgEl('rect');\n",
    "    rect.setAttribute('x', Math.min(x1, x2));\n",
    "    rect.setAttribute('y', y);\n",
    "    rect.setAttribute('width', Math.max(1, Math.abs(x2 - x1)));\n",
    "    rect.setAttribute('height', height);\n",
    "    rect.setAttribute('fill', color);\n",
    "    rect.setAttribute('opacity', '0.8');\n",
    "    svg.appendChild(rect);\n",
    "\n",
    "    if (label) {\n",
    "      var text = svgEl('text');\n",
    "      text.setAttribute('x', Math.min(x1, x2));\n",
    "      text.setAttribute('y', y - 2);\n",
    "      text.setAttribute('class', 'track-label');\n",
    "      text.textContent = label;\n",
    "      svg.appendChild(text);\n",
    "    }\n",
    "  }\n",
    "\n",
    "  function drawSelection(svg, start, end, width, height) {\n",
    "    if (!selection) return;\n",
    "    var x1 = xScale(selection.start, start, end, width, 60, 20);\n",
    "    var x2 = xScale(selection.end, start, end, width, 60, 20);\n",
    "    var rect = svgEl('rect');\n",
    "    rect.setAttribute('x', Math.min(x1, x2));\n",
    "    rect.setAttribute('y', 0);\n",
    "    rect.setAttribute('width', Math.abs(x2 - x1));\n",
    "    rect.setAttribute('height', isFinite(height) ? height : 220);\n",
    "    rect.setAttribute('fill', '#90cdf4');\n",
    "    rect.setAttribute('opacity', '0.35');\n",
    "    svg.appendChild(rect);\n",
    "  }\n",
    "\n",
    "  function drawPvalTrack(svg, locus, width, pRows, yBySource, laneH) {\n",
    "    if (!pRows || pRows.length === 0) return;\n",
    "    var THRESH = 7.3; // -log10(5e-8)\n",
    "\n",
    "    function fmtP(x) {\n",
    "      if (!isFinite(x)) return '';\n",
    "      return String(Math.round(x));\n",
    "    }\n",
    "    function getWs(r) { return (r.WIN_START !== undefined && r.WIN_START !== null) ? r.WIN_START : r.WINDOW_START; }\n",
    "    function getWe(r) { return (r.WIN_END !== undefined && r.WIN_END !== null) ? r.WIN_END : r.WINDOW_END; }\n",
    "\n",
    "    // Show only significant bins, but per source, as labels drawn on the\n",
    "    // corresponding dataset strip (no highlighting rectangles).\n",
    "    var bestBySourceWin = {};\n",
    "    for (var i = 0; i < pRows.length; i++) {\n",
    "      var r = pRows[i];\n",
    "      var source = r.SOURCE || '';\n",
    "      var ws = getWs(r);\n",
    "      var we = getWe(r);\n",
    "      if (!isFinite(ws) || !isFinite(we)) continue;\n",
    "      var v = r.MAX_NEG_LOG10_P;\n",
    "      if (!isFinite(v)) continue;\n",
    "      var key = String(ws) + '-' + String(we);\n",
    "      if (!bestBySourceWin[source]) bestBySourceWin[source] = {};\n",
    "      var cur = bestBySourceWin[source][key];\n",
    "      if (!cur || v > cur.v) bestBySourceWin[source][key] = { ws: ws, we: we, v: v };\n",
    "    }\n",
    "\n",
    "    Object.keys(bestBySourceWin).forEach(function(source) {\n",
    "      var y0 = yBySource && yBySource[source] !== undefined ? yBySource[source] : null;\n",
    "      if (y0 === null) return;\n",
    "      var winMap = bestBySourceWin[source];\n",
    "      Object.keys(winMap).forEach(function(key) {\n",
    "        var b = winMap[key];\n",
    "        if (!b || !isFinite(b.v) || b.v < THRESH) return;\n",
    "        var x1 = xScale(b.ws, locus.DISPLAY_START, locus.DISPLAY_END, width, 60, 20);\n",
    "        var x2 = xScale(b.we, locus.DISPLAY_START, locus.DISPLAY_END, width, 60, 20);\n",
    "        var xMid = (Math.min(x1, x2) + Math.max(x1, x2)) / 2;\n",
    "        var xMin = 60 + 8;\n",
    "        var xMax = (width - 20) - 8;\n",
    "\n",
    "        var txt = svgEl('text');\n",
    "        txt.setAttribute('y', y0 + (laneH || 12) - 2);\n",
    "        txt.setAttribute('text-anchor', 'middle');\n",
    "        txt.setAttribute('font-size', '6');\n",
    "        txt.setAttribute('fill', '#2d3748');\n",
    "        txt.setAttribute('stroke', '#ffffff');\n",
    "        txt.setAttribute('stroke-width', '3');\n",
    "        txt.setAttribute('paint-order', 'stroke');\n",
    "        txt.textContent = fmtP(b.v);\n",
    "        var tt = svgEl('title');\n",
    "        tt.textContent = (source ? (source + ' ') : '') + fmtInt(b.ws) + '-' + fmtInt(b.we) + '  max -log10(p)=' + fmtP(b.v);\n",
    "        txt.appendChild(tt);\n",
    "        svg.appendChild(txt);\n",
    "\n",
    "        // Clamp based on rendered text width so the label never spills\n",
    "        // outside the plotting area (stroke included).\n",
    "        try {\n",
    "          var len = (txt.getComputedTextLength ? txt.getComputedTextLength() : 0);\n",
    "          if (isFinite(len) && len > 0) {\n",
    "            var half = (len / 2) + 2;\n",
    "            var minX = xMin + half;\n",
    "            var maxX = xMax - half;\n",
    "            if (minX > maxX) xMid = (xMin + xMax) / 2;\n",
    "            else xMid = clamp(xMid, minX, maxX);\n",
    "          } else {\n",
    "            xMid = clamp(xMid, xMin, xMax);\n",
    "          }\n",
    "        } catch (e) {\n",
    "          xMid = clamp(xMid, xMin, xMax);\n",
    "        }\n",
    "        txt.setAttribute('x', xMid);\n",
    "      });\n",
    "    });\n",
    "  }\n",
    "\n",
    "  function renderViz(locus, rows, origRows, pRows) {\n",
    "    var svg = byId('viz');\n",
    "    if (!svg) return;\n",
    "    var width = 900;\n",
    "    clearSvg(svg);\n",
    "\n",
    "    // Title\n",
    "    var title = svgEl('text');\n",
    "    title.setAttribute('x', 10);\n",
    "    title.setAttribute('y', 18);\n",
    "    title.setAttribute('class', 'viz-title');\n",
    "    title.textContent = 'Display window: ' + locus.CHR + ':' + fmtInt(locus.DISPLAY_START) + '-' + fmtInt(locus.DISPLAY_END);\n",
    "    svg.appendChild(title);\n",
    "\n",
    "    drawRuler(svg, locus.DISPLAY_START, locus.DISPLAY_END, width);\n",
    "\n",
    "    // Draw original loci segments (colored by source)\n",
    "    var yOrig = 55;\n",
    "    var laneH = 12;\n",
    "    var laneGap = 10;\n",
    "    var laneStep = laneH + laneGap;\n",
    "    // group by source\n",
    "    var sources = {};\n",
    "    for (var i = 0; i < origRows.length; i++) {\n",
    "      var o = origRows[i];\n",
    "      if (!sources[o.SOURCE]) sources[o.SOURCE] = [];\n",
    "      sources[o.SOURCE].push(o);\n",
    "    }\n",
    "    var sourceKeys = Object.keys(sources).sort();\n",
    "    var yBySource = {};\n",
    "    sourceKeys.forEach(function(source, idxLane) {\n",
    "      yBySource[source] = yOrig + idxLane * laneStep;\n",
    "      sources[source].forEach(function(o) {\n",
    "        var col = COLORS && COLORS[source] ? COLORS[source] : '#2b6cb0';\n",
    "        drawInterval(svg, locus.DISPLAY_START, locus.DISPLAY_END, o.START, o.END, yOrig + idxLane * laneStep, laneH, col, null, width);\n",
    "      });\n",
    "      // label\n",
    "      var lbl = svgEl('text');\n",
    "      lbl.setAttribute('x', 10);\n",
    "      lbl.setAttribute('y', yOrig + idxLane * laneStep + 10);\n",
    "      lbl.setAttribute('class', 'track-label');\n",
    "      lbl.textContent = source;\n",
    "      svg.appendChild(lbl);\n",
    "    });\n",
    "\n",
    "    var yAfterLanes = yOrig + Math.max(1, sourceKeys.length) * laneStep;\n",
    "    var yMerged = yAfterLanes + 10;\n",
    "    var viewH = Math.max(260, yMerged + 16 + 60);\n",
    "    svg.setAttribute('viewBox', '0 0 ' + width + ' ' + viewH);\n",
    "    svg.style.height = Math.max(240, viewH) + 'px';\n",
    "\n",
    "    // Draw merged locus interval\n",
    "    drawInterval(svg, locus.DISPLAY_START, locus.DISPLAY_END, locus.MERGED_START, locus.MERGED_END, yMerged, 16, '#2b6cb0', 'Merged', width);\n",
    "\n",
    "    // Draw reported loci intervals\n",
    "    var yReported = yMerged + 28;\n",
    "    var nReported = 0;\n",
    "    for (var j = 0; j < rows.length; j++) {\n",
    "      var r = rows[j];\n",
    "      if (!isFinite(r.START) || !isFinite(r.END)) continue;\n",
    "      var st = String(r.ALIGNMENT_STATUS || '').toLowerCase();\n",
    "      var col = '#805ad5';\n",
    "      if (st === 'proximal') col = '#d69e2e';\n",
    "      drawInterval(svg, locus.DISPLAY_START, locus.DISPLAY_END, r.START, r.END, yReported, 10, col, null, width);\n",
    "      nReported++;\n",
    "    }\n",
    "    if (nReported > 0) {\n",
    "      var repLbl = svgEl('text');\n",
    "      repLbl.setAttribute('x', 10);\n",
    "      repLbl.setAttribute('y', yReported + 8);\n",
    "      repLbl.setAttribute('class', 'track-label');\n",
    "      repLbl.textContent = 'Reported';\n",
    "      svg.appendChild(repLbl);\n",
    "    }\n",
    "\n",
    "    // Draw selection overlay last so it's visible\n",
    "    drawSelection(svg, locus.DISPLAY_START, locus.DISPLAY_END, width, viewH);\n",
    "\n",
    "    // Draw per-source p-value labels on dataset strips\n",
    "    drawPvalTrack(svg, locus, width, pRows, yBySource, laneH);\n",
    "  }\n",
    "\n",
    "  function render() {\n",
    "    var locus = LOCI[idx];\n",
    "    if (!locus) return;\n",
    "\n",
    "    // normalize SPLITS to string for storage; render uses parsed\n",
    "    var rows = getLocusRows(locus.MERGED_LOCUS_ID);\n",
    "    var origRows = getOrigRows(locus.MERGED_LOCUS_ID);\n",
    "    var pRows = getPvalRows(locus.MERGED_LOCUS_ID);\n",
    "\n",
    "    renderLocusInfo(locus);\n",
    "    renderForm(locus);\n",
    "\n",
    "    var splitsArr = parseSplits(locus.SPLITS);\n",
    "    renderSplits(splitsArr);\n",
    "    renderViz(locus, rows, origRows, pRows);\n",
    "    validateInputs();\n",
    "  }\n",
    "\n",
    "  function downloadTSV() {\n",
    "    // Prepare a locus-level TSV\n",
    "    var headers = Object.keys(LOCI[0]);\n",
    "    // Ensure SPLITS is string\n",
    "    var lines = [];\n",
    "    lines.push(headers.join('\\t'));\n",
    "    for (var i = 0; i < LOCI.length; i++) {\n",
    "      var row = LOCI[i];\n",
    "      var vals = [];\n",
    "      for (var j = 0; j < headers.length; j++) {\n",
    "        vals.push(escTSV(row[headers[j]]));\n",
    "      }\n",
    "      lines.push(vals.join('\\t'));\n",
    "    }\n",
    "    var blob = new Blob([lines.join('\\n')], { type: 'text/tab-separated-values;charset=utf-8' });\n",
    "    var url = URL.createObjectURL(blob);\n",
    "    var a = document.createElement('a');\n",
    "    a.href = url;\n",
    "    a.download = DOWNLOAD_NAME || 'loci.tsv';\n",
    "    document.body.appendChild(a);\n",
    "    a.click();\n",
    "    document.body.removeChild(a);\n",
    "    URL.revokeObjectURL(url);\n",
    "    dirty = false;\n",
    "    updateDirty();\n",
    "  }\n",
    "\n",
    "  function maybeNext(delta) {\n",
    "    var next = idx + delta;\n",
    "    if (next < 0) next = 0;\n",
    "    if (next >= LOCI.length) next = LOCI.length - 1;\n",
    "    idx = next;\n",
    "    selection = null;\n",
    "    render();\n",
    "  }\n",
    "\n",
    "  function setupVizMouse() {\n",
    "    var svg = byId('viz');\n",
    "    if (!svg) return;\n",
    "\n",
    "    function evToPos(ev, locus) {\n",
    "      var rect = svg.getBoundingClientRect();\n",
    "      var x = ev.clientX - rect.left;\n",
    "      var width = rect.width;\n",
    "      var padL = 60;\n",
    "      var padR = 20;\n",
    "      var inner = width - padL - padR;\n",
    "      var frac = (x - padL) / inner;\n",
    "      frac = clamp(frac, 0, 1);\n",
    "      return locus.DISPLAY_START + frac * (locus.DISPLAY_END - locus.DISPLAY_START);\n",
    "    }\n",
    "\n",
    "    var dragging = false;\n",
    "    var startPos = null;\n",
    "\n",
    "    svg.addEventListener('mousedown', function(ev) {\n",
    "      var locus = LOCI[idx];\n",
    "      if (!locus) return;\n",
    "      dragging = true;\n",
    "      startPos = evToPos(ev, locus);\n",
    "      selection = { start: Math.round(startPos), end: Math.round(startPos) };\n",
    "      render();\n",
    "    });\n",
    "\n",
    "    svg.addEventListener('mousemove', function(ev) {\n",
    "      if (!dragging) return;\n",
    "      var locus = LOCI[idx];\n",
    "      if (!locus) return;\n",
    "      var p = evToPos(ev, locus);\n",
    "      selection.end = Math.round(p);\n",
    "      render();\n",
    "    });\n",
    "\n",
    "    window.addEventListener('mouseup', function() {\n",
    "      if (!dragging) return;\n",
    "      dragging = false;\n",
    "      if (!selection) return;\n",
    "      // Normalize ordering\n",
    "      var a = selection.start;\n",
    "      var b = selection.end;\n",
    "      selection.start = Math.min(a, b);\n",
    "      selection.end = Math.max(a, b);\n",
    "      render();\n",
    "    });\n",
    "  }\n",
    "\n",
    "  function addSplitFromSelection() {\n",
    "    var locus = LOCI[idx];\n",
    "    if (!locus) return;\n",
    "    if (!selection) { alert('Drag on the plot to select a region first.'); return; }\n",
    "    var a = Math.round(selection.start);\n",
    "    var b = Math.round(selection.end);\n",
    "    if (!isFinite(a) || !isFinite(b) || a >= b) { alert('Invalid selection.'); return; }\n",
    "    var arr = parseSplits(locus.SPLITS);\n",
    "    arr.push({ start: a, end: b });\n",
    "    locus.SPLITS = splitsToString(arr);\n",
    "    dirty = true;\n",
    "    updateDirty();\n",
    "    render();\n",
    "  }\n",
    "\n",
    "  function initEvents() {\n",
    "    byId('prev').addEventListener('click', function() { maybeNext(-1); });\n",
    "    byId('next').addEventListener('click', function() { maybeNext(1); });\n",
    "    byId('save').addEventListener('click', function() { downloadTSV(); });\n",
    "    byId('add_split').addEventListener('click', function() { addSplitFromSelection(); });\n",
    "\n",
    "    var inputs = ['refined_start', 'refined_end', 'locus_name', 'keep'];\n",
    "    inputs.forEach(function(id) {\n",
    "      var el = byId(id);\n",
    "      el.addEventListener('input', function() { \n",
    "        updateFromInputs(); \n",
    "        validateInputs(); \n",
    "      });\n",
    "      el.addEventListener('change', function() { \n",
    "        updateFromInputs(); \n",
    "        validateInputs(); \n",
    "      });\n",
    "    });\n",
    "    \n",
    "    setupVizMouse();\n",
    "    render();\n",
    "  }\n",
    "  \n",
    "  // Start everything\n",
    "  try {\n",
    "    initEvents();\n",
    "  } catch (e) {\n",
    "    showError(e && e.stack ? e.stack : e);\n",
    "  }\n",
    "})();\n"
  )
}

locus_wizard_build_css <- function() {
  "
    body { 
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif; 
      margin: 20px; 
      line-height: 1.5;
    }
    
    h1 { 
      color: #2d3748; 
      margin-bottom: 10px; 
      font-size: 24px;
    }
    
    .row { 
      display: flex; 
      gap: 20px; 
      align-items: flex-start; 
      margin-top: 20px;
    }
    
    .panel { 
      border: 1px solid #e2e8f0; 
      padding: 16px; 
      border-radius: 8px; 
      background: white;
      box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
    }
    
    .left { 
      flex: 1 1 auto; 
      min-width: 600px; 
    }
    
    .right { 
      width: 320px; 
      flex-shrink: 0;
    }
    
    .small { 
      color: #718096; 
      font-size: 12px; 
    }
    
    .controls { 
      display: flex; 
      gap: 10px; 
      align-items: center; 
      margin-bottom: 20px;
      padding: 12px;
      background: #f7fafc;
      border-radius: 6px;
    }
    
    .btn { 
      padding: 8px 16px; 
      cursor: pointer; 
      background: #4299e1;
      color: white;
      border: none;
      border-radius: 4px;
      font-size: 14px;
      font-weight: 500;
      transition: background-color 0.2s;
    }
    
    .btn:hover { 
      background: #3182ce; 
    }
    
    .btn-small { 
      padding: 2px 8px; 
      font-size: 12px; 
      cursor: pointer; 
      background: #e2e8f0;
      color: #4a5568;
      border: 1px solid #cbd5e0;
      border-radius: 3px;
    }
    
    .btn-small:hover {
      background: #cbd5e0;
    }
    
    .badge { 
      padding: 4px 10px; 
      border-radius: 12px; 
      font-size: 12px; 
      font-weight: 500;
      margin-left: auto;
    }
    
    .badge-ok { 
      background: #c6f6d5; 
      color: #22543d; 
    }
    
    .badge-warn { 
      background: #fed7d7; 
      color: #742a2a; 
    }
    
    label { 
      display: block; 
      margin-top: 12px; 
      font-weight: 600; 
      color: #4a5568;
      font-size: 14px;
    }
    
    input[type='text'], 
    input[type='number'] { 
      width: 100%; 
      box-sizing: border-box; 
      padding: 8px; 
      border: 1px solid #cbd5e0;
      border-radius: 4px;
      font-size: 14px;
      margin-top: 4px;
    }
    
    input.invalid { 
      outline: 2px solid #fc8181; 
      border-color: #fc8181;
    }
    
    input[type='checkbox'] {
      margin-top: 12px;
      width: 18px;
      height: 18px;
    }
    
    ul { 
      padding-left: 20px; 
      margin: 8px 0;
    }
    
    li {
      margin: 4px 0;
      font-size: 13px;
    }
    
    #viz { 
      width: 100%; 
      height: 240px; 
      border: 1px solid #e2e8f0; 
      border-radius: 6px; 
      background: white;
      margin: 10px 0;
    }
    
    .viz-title {
      font-weight: bold;
      fill: #2d3748;
    }
    
    .ruler-line {
      stroke: #4a5568;
      stroke-width: 1;
    }
    
    .ruler-label {
      font-size: 10px;
      fill: #718096;
    }
    
    .track-label {
      font-size: 11px;
      fill: #718096;
    }
    
    #js_error {
      display: none; 
      white-space: pre-wrap; 
      background: #fed7d7; 
      border: 1px solid #fc8181; 
      padding: 12px; 
      border-radius: 6px; 
      margin: 10px 0;
      color: #742a2a;
      font-size: 13px;
    }
    
    .locus-info-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 12px;
      margin: 15px 0;
      padding: 12px;
      background: #f7fafc;
      border-radius: 6px;
    }
    
    .info-item {
      font-size: 13px;
    }
    
    .info-label {
      font-weight: 600;
      color: #4a5568;
    }
    
    .info-value {
      color: #2d3748;
      font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
    }
  "
}
