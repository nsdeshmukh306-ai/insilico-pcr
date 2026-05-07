/**
 * plots.js — Plotly chart builders for the PCR dashboard.
 *
 * All functions accept processed data and a DOM element ID.
 * They return nothing; Plotly renders directly.
 */

const COLORS = {
  blue:   '#1e6ba8',
  teal:   '#0d7c6e',
  amber:  '#c4730a',
  red:    '#c0392b',
  green:  '#1e8449',
  gray:   '#6b7280',
  blue100:'#e8f4fd',
  seq: ['#1e6ba8','#0d7c6e','#c4730a','#8b5cf6','#e11d48'],
};

const LAYOUT_BASE = {
  paper_bgcolor: 'transparent',
  plot_bgcolor:  '#f9fafb',
  font: { family: 'Inter, sans-serif', size: 11, color: '#374151' },
  margin: { l: 48, r: 20, t: 30, b: 44 },
  hoverlabel: {
    bgcolor: '#111827',
    bordercolor: '#1e6ba8',
    font: { color: '#ffffff', size: 11 },
  },
  xaxis: { gridcolor: '#e5e7eb', zerolinecolor: '#d1d5db', linecolor: '#d1d5db' },
  yaxis: { gridcolor: '#e5e7eb', zerolinecolor: '#d1d5db', linecolor: '#d1d5db' },
};

const CONFIG = {
  displayModeBar: true,
  displaylogo: false,
  modeBarButtonsToRemove: ['select2d','lasso2d','autoScale2d'],
  responsive: true,
};

function _layout(overrides = {}) {
  return Plotly.react ? Object.assign({}, LAYOUT_BASE, overrides) : Object.assign({}, LAYOUT_BASE, overrides);
}

// ── Score bar chart ───────────────────────────────────────────────────────────
function plotScoreBreakdown(containerId, scoreComponents, pairName) {
  const labels = ['Binding', 'Tm compat.', 'GC content', 'Length'];
  const values = [
    (scoreComponents.s_bind || 0) * 100,
    (scoreComponents.s_tm   || 0) * 100,
    (scoreComponents.s_gc   || 0) * 100,
    (scoreComponents.s_len  || 0) * 100,
  ];
  const penLabels = ['Mismatch penalty', 'Off-target penalty'];
  const penValues = [
    (scoreComponents.p_mm   || 0) * 100,
    (scoreComponents.p_offt || 0) * 100,
  ];

  const traces = [
    {
      x: values,
      y: labels,
      type: 'bar',
      orientation: 'h',
      name: 'Positive',
      marker: { color: COLORS.blue, opacity: .85 },
      hovertemplate: '%{y}: %{x:.1f}%<extra></extra>',
    },
    {
      x: penValues,
      y: penLabels,
      type: 'bar',
      orientation: 'h',
      name: 'Penalty',
      marker: { color: COLORS.red, opacity: .75 },
      hovertemplate: '%{y}: %{x:.1f}%<extra></extra>',
    },
  ];

  const layout = _layout({
    title: { text: `Score breakdown — ${pairName}`, font: { size: 12 } },
    margin: { l: 120, r: 30, t: 36, b: 40 },
    xaxis: { title: 'Score (0–100)', range: [0, 105] },
    yaxis: { automargin: true },
    barmode: 'overlay',
    showlegend: true,
    legend: { orientation: 'h', y: -0.15, x: 0 },
    height: 220,
  });

  Plotly.newPlot(containerId, traces, layout, CONFIG);
}

// ── Tm comparison gauge-style bar ─────────────────────────────────────────────
function plotTmComparison(containerId, pairs) {
  const names = [], fwdTms = [], revTms = [];
  pairs.forEach(p => {
    const amps = p.amplicons || [];
    const best = amps.find(a => a.is_intended) || amps[0];
    if (!best) return;
    names.push(p.name);
    fwdTms.push(+(best.fwd_binding?.tm_celsius || 0).toFixed(1));
    revTms.push(+(best.rev_binding?.tm_celsius  || 0).toFixed(1));
  });

  // Filter out physically unrealistic Tm values from example data
  const validFwd = fwdTms.map(t => (t > 0 && t < 120) ? t : null);
  const validRev = revTms.map(t => (t > 0 && t < 120) ? t : null);

  const traces = [
    {
      x: names, y: validFwd, type: 'bar', name: 'Fwd Tm (°C)',
      marker: { color: COLORS.blue, opacity: .85 },
      hovertemplate: '%{x}<br>Fwd Tm: %{y:.1f}°C<extra></extra>',
    },
    {
      x: names, y: validRev, type: 'bar', name: 'Rev Tm (°C)',
      marker: { color: COLORS.teal, opacity: .85 },
      hovertemplate: '%{x}<br>Rev Tm: %{y:.1f}°C<extra></extra>',
    },
  ];

  const layout = _layout({
    title: { text: 'Primer Tm comparison', font: { size: 12 } },
    barmode: 'group',
    xaxis: { title: 'Primer pair' },
    yaxis: { title: 'Tm (°C)', range: [45, 80] },
    shapes: [
      { type: 'rect', xref: 'paper', x0: 0, x1: 1, y0: 55, y1: 68,
        fillcolor: 'rgba(30,132,73,.06)', line: { width: 0 }, layer: 'below' },
    ],
    annotations: [{
      xref: 'paper', yref: 'y', x: 1.01, y: 61.5,
      text: 'Optimal<br>55–68°C', showarrow: false,
      font: { size: 9, color: COLORS.green }, textangle: 0,
    }],
    showlegend: true,
    legend: { orientation: 'h', y: -0.2, x: 0 },
    height: 220,
  });

  Plotly.newPlot(containerId, traces, layout, CONFIG);
}

// ── Amplicon length distribution ──────────────────────────────────────────────
function plotAmpliconLengths(containerId, pairs) {
  const lengths = [];
  const labels  = [];
  pairs.forEach(p => {
    (p.amplicons || []).forEach(a => {
      lengths.push(a.length);
      labels.push(`${p.name} (${a.seq_id}:${a.start}–${a.end})`);
    });
  });

  if (!lengths.length) return;

  const trace = {
    x: lengths,
    type: 'histogram',
    nbinsx: Math.min(30, Math.max(5, Math.ceil(lengths.length / 3))),
    marker: { color: COLORS.blue, opacity: .8, line: { color: '#1a3a5c', width: 0.5 } },
    hovertemplate: 'Length: %{x} bp<br>Count: %{y}<extra></extra>',
  };

  const layout = _layout({
    title: { text: 'Amplicon length distribution', font: { size: 12 } },
    xaxis: { title: 'Amplicon length (bp)' },
    yaxis: { title: 'Count' },
    height: 200,
  });

  Plotly.newPlot(containerId, [trace], layout, CONFIG);
}

// ── GC content scatter ────────────────────────────────────────────────────────
function plotGCScatter(containerId, pairs) {
  const x = [], y = [], text = [], colors = [];
  pairs.forEach((p, pi) => {
    (p.amplicons || []).forEach(a => {
      x.push(a.length);
      y.push((a.gc_fraction * 100).toFixed(1));
      text.push(`${p.name}<br>${a.seq_id}:${a.start}–${a.end}<br>Score: ${a.final_score}`);
      colors.push(a.is_intended ? COLORS.blue : COLORS.gray);
    });
  });

  if (!x.length) return;

  const trace = {
    x, y: y.map(Number), text,
    mode: 'markers',
    type: 'scatter',
    marker: { color: colors, size: 8, opacity: .8,
              line: { color: 'rgba(255,255,255,.6)', width: 1 } },
    hovertemplate: '%{text}<br>Length: %{x} bp<br>GC: %{y:.1f}%<extra></extra>',
  };

  const layout = _layout({
    title: { text: 'Amplicon length vs GC content', font: { size: 12 } },
    xaxis: { title: 'Amplicon length (bp)' },
    yaxis: { title: 'GC content (%)', range: [20, 80] },
    shapes: [
      { type: 'rect', xref: 'paper', x0: 0, x1: 1, y0: 40, y1: 60,
        fillcolor: 'rgba(30,132,73,.05)', line: { width: 0 }, layer: 'below' },
    ],
    height: 220,
    showlegend: false,
  });

  Plotly.newPlot(containerId, [trace], layout, CONFIG);
}

// ── Genome position track ─────────────────────────────────────────────────────
function plotGenomeTrack(containerId, pairs) {
  const traces = [];

  pairs.forEach((p, pi) => {
    const color = COLORS.seq[pi % COLORS.seq.length];
    (p.amplicons || []).forEach((a, ai) => {
      traces.push({
        x: [a.start, a.end],
        y: [pi, pi],
        mode: 'lines+markers',
        type: 'scatter',
        name: a.is_intended ? p.name : `${p.name} (off-target)`,
        line: { color: a.is_intended ? color : COLORS.gray, width: a.is_intended ? 6 : 3 },
        marker: { size: 7, color: color, symbol: 'diamond' },
        hovertemplate:
          `${p.name}<br>${a.seq_id}: ${a.start}–${a.end}<br>` +
          `Length: ${a.length} bp | Score: ${a.final_score}<extra></extra>`,
        showlegend: ai === 0,
      });
    });
  });

  if (!traces.length) return;

  const pairNames = pairs.map(p => p.name);
  const layout = _layout({
    title: { text: 'Amplicon positions on genome', font: { size: 12 } },
    xaxis: { title: 'Genomic position (bp)', tickformat: ',d' },
    yaxis: {
      tickvals: pairs.map((_, i) => i),
      ticktext: pairNames,
      gridcolor: '#e5e7eb',
      range: [-0.7, pairs.length - 0.3],
    },
    height: Math.max(160, 80 + pairs.length * 50),
    showlegend: false,
  });

  Plotly.newPlot(containerId, traces, layout, CONFIG);
}

// ── Off-target score histogram ────────────────────────────────────────────────
function plotOfftargetScores(containerId, hits) {
  if (!hits || !hits.length) {
    document.getElementById(containerId).innerHTML =
      '<div class="empty-state"><div class="empty-state-title">No off-target hits</div></div>';
    return;
  }

  const scores = hits.map(h => h.offtarget_score);
  const trace = {
    x: scores,
    type: 'histogram',
    nbinsx: 15,
    marker: { color: COLORS.amber, opacity: .8 },
    hovertemplate: 'Score: %{x:.1f}<br>Count: %{y}<extra></extra>',
  };

  const layout = _layout({
    title: { text: 'Off-target score distribution', font: { size: 12 } },
    xaxis: { title: 'Off-target score' },
    yaxis: { title: 'Count' },
    height: 190,
  });

  Plotly.newPlot(containerId, [trace], layout, CONFIG);
}

// ── Thermodynamics radar (ΔH, ΔS, ΔG, Tm) ────────────────────────────────────
function plotThermoRadar(containerId, fwdPrimer, revPrimer, pairName) {
  const categories = ['Tm (norm.)', 'ΔG (norm.)', 'GC content', 'Length ok', 'Clamp ok'];

  function normalise(p) {
    const tm = p.tm_celsius || 0;
    const tmN = tm > 0 ? Math.min(1, Math.max(0, (tm - 45) / 35)) : 0;
    const dg  = p.delta_g_37 || 0;
    const dgN = Math.min(1, Math.max(0, (-dg) / 40));
    const gcN = p.gc_fraction || 0;
    const lenN = p.length_ok ? 1 : 0.3;
    const clN  = p.gc_clamp_ok ? 1 : 0.4;
    return [tmN, dgN, gcN, lenN, clN];
  }

  const fwdVals = normalise(fwdPrimer);
  const revVals = normalise(revPrimer);

  const traces = [
    {
      type: 'scatterpolar',
      r: [...fwdVals, fwdVals[0]],
      theta: [...categories, categories[0]],
      fill: 'toself',
      name: 'Forward',
      marker: { color: COLORS.blue },
      line: { color: COLORS.blue },
      fillcolor: 'rgba(30,107,168,.15)',
    },
    {
      type: 'scatterpolar',
      r: [...revVals, revVals[0]],
      theta: [...categories, categories[0]],
      fill: 'toself',
      name: 'Reverse',
      marker: { color: COLORS.teal },
      line: { color: COLORS.teal },
      fillcolor: 'rgba(13,124,110,.15)',
    },
  ];

  const layout = {
    polar: {
      radialaxis: { visible: true, range: [0, 1], tickvals: [0.25, 0.5, 0.75, 1] },
      bgcolor: '#f9fafb',
    },
    paper_bgcolor: 'transparent',
    font: { family: 'Inter, sans-serif', size: 11, color: '#374151' },
    showlegend: true,
    legend: { orientation: 'h', y: -0.15 },
    margin: { l: 40, r: 40, t: 20, b: 40 },
    height: 260,
    title: { text: `Primer quality radar — ${pairName}`, font: { size: 12 } },
  };

  Plotly.newPlot(containerId, traces, layout, CONFIG);
}

// ── Mismatch position heatmap ─────────────────────────────────────────────────
function plotMismatchPositions(containerId, bindingSite, primerLabel) {
  if (!bindingSite) return;
  const seq  = bindingSite.aligned_template || '';
  const mmps = bindingSite.mismatch_positions || [];
  const primer = bindingSite.primer_seq || '';
  const n = primer.length || seq.length;

  const x = Array.from({length: n}, (_, i) => i + 1);
  const y = [primerLabel];
  const z = [Array.from({length: n}, (_, i) => mmps.includes(i) ? 1 : 0)];

  const layout = _layout({
    title: { text: `Mismatch map — ${primerLabel}`, font: { size: 12 } },
    xaxis: { title: 'Position', dtick: 2 },
    yaxis: { visible: false },
    height: 100,
    margin: { l: 30, r: 20, t: 36, b: 40 },
    colorscale: [[0,'#e8f4fd'],[1,'#c0392b']],
  });

  Plotly.newPlot(containerId, [{
    z, x, y,
    type: 'heatmap',
    colorscale: [[0,'#e8f4fd'],[1,'#c0392b']],
    showscale: false,
    hovertemplate: 'Position %{x}: %{z}<extra></extra>',
  }], layout, CONFIG);
}

// ── Live param sweep result plot ──────────────────────────────────────────────
function plotLiveSweep(containerId, sweepResults) {
  // sweepResults = [{param_val, num_amplicons, top_score, top_length}]
  if (!sweepResults || !sweepResults.length) return;

  const x = sweepResults.map(r => r.param_val);
  const scores = sweepResults.map(r => r.top_score);
  const lengths = sweepResults.map(r => r.top_length);

  const traces = [
    {
      x, y: scores, type: 'scatter', mode: 'lines+markers',
      name: 'Top score', yaxis: 'y',
      line: { color: COLORS.blue, width: 2 },
      marker: { color: COLORS.blue, size: 7 },
      hovertemplate: 'Param: %{x}<br>Score: %{y:.1f}<extra></extra>',
    },
    {
      x, y: lengths, type: 'scatter', mode: 'lines+markers',
      name: 'Amplicon length', yaxis: 'y2',
      line: { color: COLORS.teal, width: 2, dash: 'dot' },
      marker: { color: COLORS.teal, size: 7 },
      hovertemplate: 'Param: %{x}<br>Length: %{y} bp<extra></extra>',
    },
  ];

  const layout = _layout({
    title: { text: 'Parameter sweep result', font: { size: 12 } },
    xaxis: { title: 'Parameter value' },
    yaxis: { title: 'Score (0–100)', side: 'left' },
    yaxis2: { title: 'Amplicon length (bp)', overlaying: 'y', side: 'right', showgrid: false },
    showlegend: true,
    legend: { orientation: 'h', y: -0.2 },
    height: 230,
  });

  Plotly.newPlot(containerId, traces, layout, CONFIG);
}
