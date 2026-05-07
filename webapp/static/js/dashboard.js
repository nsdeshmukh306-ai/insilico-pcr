/**
 * dashboard.js — Main controller for the in-silico PCR dashboard.
 *
 * Responsibilities:
 *   • Run form handling (demo load, custom run, genome input)
 *   • Data wiring to all 10 UI components
 *   • Tab and pair-selector navigation
 *   • Live parameter experiment panel
 *   • Export actions
 */

'use strict';

// ── State ─────────────────────────────────────────────────────────────────────
let _state = {
  result: null,          // full API JSON result
  activePairIdx: 0,      // which primer pair is selected
  activeAmpIdx: 0,       // which amplicon is selected
  genomeString: '',      // raw genome text for live param re-runs
};

// ── Boot ──────────────────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
  _bindEvents();
  _setupTabs();
});

// ── Event wiring ──────────────────────────────────────────────────────────────
function _bindEvents() {
  // Run button
  document.getElementById('btn-run')?.addEventListener('click', _handleRun);

  // Demo button
  document.getElementById('btn-demo')?.addEventListener('click', _handleDemo);

  // Genome source toggle
  document.querySelectorAll('input[name="genome-src"]').forEach(r => {
    r.addEventListener('change', _toggleGenomeInput);
  });

  // Live param sliders
  document.querySelectorAll('.live-slider').forEach(sl => {
    sl.addEventListener('input', _onSliderChange);
  });

  // Live run button
  document.getElementById('btn-live-run')?.addEventListener('click', _handleLiveRun);

  // Export buttons
  document.getElementById('btn-export-json')?.addEventListener('click', () => _export('json'));
  document.getElementById('btn-export-csv')?.addEventListener('click',  () => _export('csv'));
  document.getElementById('btn-export-html')?.addEventListener('click', _exportHtml);

  // Primer pair tabs (delegated)
  document.getElementById('pair-tabs')?.addEventListener('click', e => {
    const btn = e.target.closest('.pair-tab');
    if (!btn) return;
    _state.activePairIdx = +btn.dataset.idx;
    _state.activeAmpIdx  = 0;
    document.querySelectorAll('.pair-tab').forEach(b => b.classList.toggle('active', b === btn));
    _renderPairView();
  });

  // Amplicon table row click (delegated)
  document.getElementById('amplicon-table-body')?.addEventListener('click', e => {
    const row = e.target.closest('tr[data-amp-idx]');
    if (!row) return;
    _state.activeAmpIdx = +row.dataset.ampIdx;
    document.querySelectorAll('#amplicon-table-body tr').forEach(r =>
      r.classList.toggle('selected', r === row));
    _renderAmpliconDetail();
  });
}

function _setupTabs() {
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const group = btn.dataset.tabGroup;
      const target = btn.dataset.tab;
      document.querySelectorAll(`.tab-btn[data-tab-group="${group}"]`)
        .forEach(b => b.classList.toggle('active', b === btn));
      document.querySelectorAll(`.tab-pane[data-tab-group="${group}"]`)
        .forEach(p => p.classList.toggle('active', p.dataset.tab === target));
    });
  });
}

function _toggleGenomeInput() {
  const src = document.querySelector('input[name="genome-src"]:checked')?.value;
  document.getElementById('genome-string-wrap')?.classList.toggle('hidden', src !== 'string');
  document.getElementById('genome-fasta-wrap')?.classList.toggle('hidden', src !== 'fasta');
}

// ── Status helpers ────────────────────────────────────────────────────────────
function _status(msg, type = 'info', spinner = false) {
  const bar = document.getElementById('status-bar');
  if (!bar) return;
  bar.className = `show ${type}`;
  bar.innerHTML = (spinner ? '<div class="spinner"></div>' : '') + `<span>${msg}</span>`;
}

function _clearStatus() {
  const bar = document.getElementById('status-bar');
  if (bar) { bar.className = ''; bar.innerHTML = ''; }
}

let _dismissTimer = null;
function _autoDismissStatus(ms = 3000) {
  clearTimeout(_dismissTimer);
  _dismissTimer = setTimeout(_clearStatus, ms);
}

// ── Run handlers ──────────────────────────────────────────────────────────────
async function _handleDemo() {
  _status('Loading demo results…', 'info', true);
  document.getElementById('btn-demo').disabled = true;
  try {
    // Load results and genome string in parallel
    const [demoRes, genomeRes] = await Promise.all([
      fetch('/api/demo'),
      fetch('/api/demo/genome'),
    ]);
    if (!demoRes.ok) throw new Error(`Demo failed: ${demoRes.status}`);
    const data = await demoRes.json();

    // Store genome so the live panel can re-run without the user pasting it
    if (genomeRes.ok) {
      const { genome_string } = await genomeRes.json();
      _state.genomeString = genome_string || '';
    }

    _loadResult(data);
    _status(`Demo loaded — ${data.primer_pairs?.length || 0} primer pair(s)`, 'success');
    _autoDismissStatus(4000);
  } catch (err) {
    _status('Failed to load demo: ' + err.message, 'error');
  } finally {
    document.getElementById('btn-demo').disabled = false;
  }
}

async function _handleRun() {
  const fwd = document.getElementById('fwd-primer')?.value.trim();
  const rev = document.getElementById('rev-primer')?.value.trim();
  const pairName = document.getElementById('pair-name')?.value.trim() || 'pair_1';

  if (!fwd || !rev) {
    _status('Please enter both forward and reverse primer sequences.', 'error');
    return;
  }

  // Genome source
  const src = document.querySelector('input[name="genome-src"]:checked')?.value || 'string';
  let genomeString = '';

  if (src === 'string') {
    genomeString = document.getElementById('genome-string-input')?.value.trim() || '';
    // Strip FASTA headers if pasted
    genomeString = genomeString.split('\n').filter(l => !l.startsWith('>')).join('');
    _state.genomeString = genomeString;
  } else {
    const file = document.getElementById('genome-fasta-file')?.files[0];
    if (!file) { _status('Please select a FASTA file.', 'error'); return; }
    const text = await file.text();
    genomeString = text.split('\n').filter(l => !l.startsWith('>')).join('');
    _state.genomeString = genomeString;
  }

  if (!genomeString) {
    _status('Genome sequence is empty.', 'error');
    return;
  }

  const params = _collectRunParams();

  _status('Running pipeline…', 'info', true);
  document.getElementById('btn-run').disabled = true;

  try {
    const body = {
      primers: [{ name: pairName, forward: fwd, reverse: rev }],
      genome_string: genomeString,
      ...params,
    };
    const res = await fetch('/api/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      // FastAPI validation errors come as {detail: [{msg, loc, ...}]}
      const msg = Array.isArray(err.detail)
        ? err.detail.map(e => `${e.loc?.slice(-1)[0] ?? ''}: ${e.msg}`).join('; ')
        : (err.detail || err.error || res.statusText);
      throw new Error(msg);
    }
    const data = await res.json();
    _loadResult(data);
    _status(`Done in ${data.elapsed_seconds}s — ${data.primer_pairs?.length || 0} pair(s)`, 'success');
    _autoDismissStatus(5000);
  } catch (err) {
    _status('Run failed: ' + err.message, 'error');
  } finally {
    document.getElementById('btn-run').disabled = false;
  }
}

function _collectRunParams() {
  const val = id => parseFloat(document.getElementById(id)?.value) || undefined;
  const chk = id => document.getElementById(id)?.checked;
  return {
    max_mismatches:    parseInt(document.getElementById('param-max-mm')?.value) || 3,
    min_amplicon:      val('param-min-amp') || 50,
    max_amplicon:      val('param-max-amp') || 3000,
    na_conc_mm:        val('param-na')      || 50,
    mg_conc_mm:        val('param-mg')      || 0,
    dntp_conc_mm:      val('param-dntp')    || 0,
    primer_conc_nm:    val('param-primer-c')|| 250,
    three_prime_strict:chk('param-3p-strict') ?? true,
    run_hairpin:       chk('param-hairpin') ?? false,
    run_dimer:         chk('param-dimer')   ?? false,
  };
}

// ── Data loading ──────────────────────────────────────────────────────────────
function _loadResult(data) {
  _state.result       = data;
  _state.activePairIdx = 0;
  _state.activeAmpIdx  = 0;

  _renderSummaryStats(data);
  _renderPairTabs(data.primer_pairs || []);
  _renderPairView();
  _renderOverviewPlots(data.primer_pairs || []);
  _populateLivePanel();

  // JSON preview for the export tab
  const pre = document.getElementById('json-preview');
  if (pre) {
    const str = JSON.stringify(data, null, 2);
    pre.textContent = str.slice(0, 10000) + (str.length > 10000 ? '\n\n[truncated — download full JSON]' : '');
  }

  document.getElementById('results-section')?.classList.add('visible');
  window.scrollTo({ top: document.getElementById('results-section')?.offsetTop - 80, behavior: 'smooth' });
}

// ── Summary stats ─────────────────────────────────────────────────────────────
function _renderSummaryStats(data) {
  const pairs = data.primer_pairs || [];
  let totalAmp = 0, totalOT = 0, maxScore = 0;

  pairs.forEach(p => {
    (p.amplicons || []).forEach(a => {
      totalAmp++;
      if (a.final_score > maxScore) maxScore = a.final_score;
    });
    totalOT += (p.offtarget_summary?.total || 0);
  });

  _setEl('stat-pairs',   pairs.length);
  _setEl('stat-amps',    totalAmp);
  _setEl('stat-offtarget', totalOT);
  _setEl('stat-top-score', maxScore.toFixed(1));
  _setEl('stat-elapsed', (data.elapsed_seconds || 0).toFixed(2) + 's');
  _setEl('stat-version', data.run_info?.pipeline_version || '—');
}

// ── Primer pair tabs ──────────────────────────────────────────────────────────
function _renderPairTabs(pairs) {
  const container = document.getElementById('pair-tabs');
  if (!container) return;
  container.innerHTML = pairs.map((p, i) =>
    `<button class="pair-tab${i === 0 ? ' active' : ''}" data-idx="${i}">${p.name}</button>`
  ).join('');
}

// ── Pair-level view ───────────────────────────────────────────────────────────
function _renderPairView() {
  const pairs = _state.result?.primer_pairs || [];
  const pair  = pairs[_state.activePairIdx];
  if (!pair) return;

  _renderPrimerQuality(pair);
  _renderAmpliconTable(pair);
  _renderAmpliconDetail();
  _renderOfftargetPanel(pair);
  _renderThermodynamicsPanel(pair);
}

// ── Primer quality panel ──────────────────────────────────────────────────────
function _renderPrimerQuality(pair) {
  const fp = pair.forward_primer || {};
  const rp = pair.reverse_primer || {};

  function primerBlock(p, label) {
    const ok = cls => cls ? 'badge-green' : 'badge-red';
    const tm = (p.tm_celsius > 0 && p.tm_celsius < 120) ? p.tm_celsius.toFixed(1) : '—';
    const dg = p.delta_g_37?.toFixed(2) ?? '—';
    return `
      <div class="thermo-card">
        <div class="thermo-header">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
            <path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5"/>
          </svg>
          ${label} &nbsp;
          <span class="badge badge-gray text-mono">${p.sequence || '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Length</span>
          <span class="thermo-val">${p.length || '—'} bp
            <span class="badge ${ok(p.length_ok)}">${p.length_ok ? '✓' : '✗'}</span>
          </span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">GC content</span>
          <span class="thermo-val">${p.gc_fraction != null ? (p.gc_fraction*100).toFixed(1)+'%' : '—'}
            <span class="badge ${ok(p.gc_ok)}">${p.gc_ok ? '✓' : '✗'}</span>
          </span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Tm (SantaLucia)</span>
          <span class="thermo-val">${tm}°C
            <span class="badge ${ok(p.tm_ok)}">${p.tm_ok ? '✓' : '✗'}</span>
          </span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">ΔG₃₇ (kcal/mol)</span>
          <span class="thermo-val">${dg}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">GC clamp</span>
          <span class="thermo-val"><span class="badge ${ok(p.gc_clamp_ok)}">${p.gc_clamp_ok ? 'OK' : 'Missing'}</span></span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Max run</span>
          <span class="thermo-val">${p.max_run ?? '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Low complexity</span>
          <span class="thermo-val"><span class="badge ${p.low_complexity ? 'badge-red' : 'badge-green'}">${p.low_complexity ? 'Yes' : 'No'}</span></span>
        </div>
      </div>
    `;
  }

  _setElHtml('primer-quality-content',
    `<div class="thermo-grid">${primerBlock(fp, 'Forward')}&nbsp;${primerBlock(rp, 'Reverse')}</div>`
  );

  // Radar chart
  plotThermoRadar('primer-radar-plot', fp, rp, pair.name);
}

// ── Amplicon table ────────────────────────────────────────────────────────────
function _renderAmpliconTable(pair) {
  const amps  = pair.amplicons || [];
  const tbody = document.getElementById('amplicon-table-body');
  if (!tbody) return;

  if (!amps.length) {
    tbody.innerHTML = '<tr><td colspan="9" class="text-center text-muted" style="padding:20px">No amplicons found</td></tr>';
    return;
  }

  tbody.innerHTML = amps.map((a, i) => {
    const sc   = a.score_components || {};
    const pct  = (a.final_score || 0).toFixed(1);
    const fill = scoreColor(a.final_score || 0);
    const fb   = a.fwd_binding || {};
    const rb   = a.rev_binding || {};
    const fwdBadge = alignmentBadge(fb.mismatch_count || 0, fb.three_prime_mm);
    const revBadge = alignmentBadge(rb.mismatch_count || 0, rb.three_prime_mm);
    const intBadge = a.is_intended
      ? '<span class="badge badge-blue">Target</span>'
      : '<span class="badge badge-gray">Off-target</span>';

    return `
      <tr data-amp-idx="${i}" class="${a.is_intended ? 'intended' : ''}" style="cursor:pointer">
        <td>${a.rank ?? i+1}</td>
        <td>${a.seq_id || '—'}</td>
        <td class="mono">${a.start?.toLocaleString()}–${a.end?.toLocaleString()}</td>
        <td><strong>${a.length}</strong> bp</td>
        <td>${(a.gc_fraction != null ? (a.gc_fraction*100).toFixed(1) : '—')}%</td>
        <td>
          <div class="score-bar">
            <span>${pct}</span>
            <div class="score-track"><div class="score-fill" style="width:${pct}%;background:${fill}"></div></div>
          </div>
        </td>
        <td>${fwdBadge}</td>
        <td>${revBadge}</td>
        <td>${intBadge}</td>
      </tr>
    `;
  }).join('');

  // Score breakdown chart for this pair
  const bestAmp = amps.find(a => a.is_intended) || amps[0];
  if (bestAmp?.score_components) {
    plotScoreBreakdown('score-breakdown-plot', bestAmp.score_components, pair.name);
  }
}

// ── Amplicon detail (alignment + thermo) ──────────────────────────────────────
function _renderAmpliconDetail() {
  const pairs = _state.result?.primer_pairs || [];
  const pair  = pairs[_state.activePairIdx];
  if (!pair) return;
  const amps = pair.amplicons || [];
  const amp  = amps[_state.activeAmpIdx] || amps[0];
  if (!amp) return;

  // Alignment viewer
  renderAlignment('alignment-viewer-content', amp.fwd_binding, amp.rev_binding, amp.sequence);

  // Thermodynamics deep-dive
  _renderThermoDeepdive(amp, pair);
}

function _renderThermoDeepdive(amp, pair) {
  const fb = amp.fwd_binding || {};
  const rb = amp.rev_binding || {};
  const fp = pair.forward_primer || {};
  const rp = pair.reverse_primer || {};

  function tmRow(label, val) {
    const display = (val && val > 0 && val < 120) ? val.toFixed(2) + '°C' : '—';
    return `<div class="thermo-metric"><span class="thermo-key">${label}</span><span class="thermo-val">${display}</span></div>`;
  }

  const html = `
    <div class="thermo-grid">
      <div class="thermo-card">
        <div class="thermo-header">⇢ Forward binding site</div>
        ${tmRow('Tm at binding', fb.tm_celsius)}
        <div class="thermo-metric">
          <span class="thermo-key">ΔG₃₇</span>
          <span class="thermo-val">${fb.delta_g_kcal != null ? fb.delta_g_kcal.toFixed(2)+' kcal/mol' : '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">SW score</span>
          <span class="thermo-val">${fb.sw_score != null ? fb.sw_score.toFixed(1) : '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Binding score</span>
          <span class="thermo-val">${fb.binding_score != null ? fb.binding_score.toFixed(1)+'/100' : '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Mismatches</span>
          <span class="thermo-val">${alignmentBadge(fb.mismatch_count||0, fb.three_prime_mm)}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Strand</span>
          <span class="thermo-val">${fb.strand || '+'}</span>
        </div>
      </div>
      <div class="thermo-card">
        <div class="thermo-header">⇠ Reverse binding site</div>
        ${tmRow('Tm at binding', rb.tm_celsius)}
        <div class="thermo-metric">
          <span class="thermo-key">ΔG₃₇</span>
          <span class="thermo-val">${rb.delta_g_kcal != null ? rb.delta_g_kcal.toFixed(2)+' kcal/mol' : '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">SW score</span>
          <span class="thermo-val">${rb.sw_score != null ? rb.sw_score.toFixed(1) : '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Binding score</span>
          <span class="thermo-val">${rb.binding_score != null ? rb.binding_score.toFixed(1)+'/100' : '—'}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Mismatches</span>
          <span class="thermo-val">${alignmentBadge(rb.mismatch_count||0, rb.three_prime_mm)}</span>
        </div>
        <div class="thermo-metric">
          <span class="thermo-key">Strand</span>
          <span class="thermo-val">${rb.strand || '-'}</span>
        </div>
      </div>
    </div>
    <div class="mt-3 thermo-card">
      <div class="thermo-header">Primer thermodynamics (designed sequence)</div>
      <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px">
        <div>
          ${tmRow('Fwd Tm (perfect match)', fp.tm_celsius)}
          <div class="thermo-metric">
            <span class="thermo-key">Fwd ΔG₃₇</span>
            <span class="thermo-val">${fp.delta_g_37?.toFixed(2) ?? '—'} kcal/mol</span>
          </div>
        </div>
        <div>
          ${tmRow('Rev Tm (perfect match)', rp.tm_celsius)}
          <div class="thermo-metric">
            <span class="thermo-key">Rev ΔG₃₇</span>
            <span class="thermo-val">${rp.delta_g_37?.toFixed(2) ?? '—'} kcal/mol</span>
          </div>
        </div>
      </div>
    </div>
  `;
  _setElHtml('thermo-deepdive', html);
}

// ── Off-target panel ──────────────────────────────────────────────────────────
function _renderOfftargetPanel(pair) {
  const ot = pair.offtarget_summary || {};
  const hits = ot.hits || [];

  _setEl('ot-total',    ot.total ?? '—');
  _setEl('ot-high',     ot.high_risk ?? 0);
  _setEl('ot-spec',     ot.specificity_index != null ? ot.specificity_index.toFixed(1)+'%' : '—');

  // Hits list
  const container = document.getElementById('offtarget-list');
  if (container) {
    if (!hits.length) {
      container.innerHTML = '<div class="empty-state"><div class="empty-state-title">No off-target hits detected</div><div class="empty-state-sub">Excellent primer specificity</div></div>';
    } else {
      container.innerHTML = hits.map(h => `
        <div class="ot-item">
          <div class="ot-item-header">
            <span class="text-sm font-bold">${h.seq_id || '—'}: ${h.start?.toLocaleString()}–${h.end?.toLocaleString()}</span>
            <span class="badge ${h.offtarget_score > 60 ? 'badge-red' : h.offtarget_score > 30 ? 'badge-amber' : 'badge-gray'}">
              Score ${h.offtarget_score?.toFixed(1)}
            </span>
          </div>
          <div class="text-xs text-muted">Size: ${h.size} bp &nbsp;·&nbsp; GC: ${h.gc_fraction != null ? (h.gc_fraction*100).toFixed(1)+'%' : '—'}</div>
          <div class="ot-reasons">
            ${(h.reasons || []).map(r => `<span class="badge badge-amber">${r.replace(/_/g,' ')}</span>`).join('')}
          </div>
        </div>
      `).join('');
    }
  }

  plotOfftargetScores('offtarget-plot', hits);
}

// ── Thermodynamics panel (pair-level overview) ────────────────────────────────
function _renderThermodynamicsPanel(pair) {
  plotTmComparison('tm-comparison-plot', _state.result?.primer_pairs || []);
}

// ── Overview / global plots ───────────────────────────────────────────────────
function _renderOverviewPlots(pairs) {
  plotAmpliconLengths('amp-length-plot', pairs);
  plotGCScatter('gc-scatter-plot', pairs);
  plotGenomeTrack('genome-track-plot', pairs);
}

// ── Live parameters panel ─────────────────────────────────────────────────────
function _populateLivePanel() {
  const pairs = _state.result?.primer_pairs || [];
  const pair  = pairs[_state.activePairIdx];
  if (!pair) return;

  // Pre-fill primer sequences
  _setVal('live-fwd', pair.forward_primer?.sequence || '');
  _setVal('live-rev', pair.reverse_primer?.sequence || '');
  _setVal('live-pair-name', pair.name);
}

function _onSliderChange(e) {
  const sl  = e.target;
  const out = document.getElementById(sl.id + '-val');
  if (out) out.textContent = sl.value;
}

async function _handleLiveRun() {
  if (!_state.result) { _status('Run or load demo first.', 'error'); return; }

  // If genome came from demo we may not have stored the raw string; warn clearly.
  if (!_state.genomeString) {
    _status('Genome sequence unavailable for live re-run — paste it in the run panel first.', 'error');
    return;
  }

  const fwd = _val('live-fwd').trim().toUpperCase();
  const rev = _val('live-rev').trim().toUpperCase();

  if (!fwd || !rev) { _status('Enter both primer sequences in the Live panel.', 'error'); return; }

  const payload = {
    primer_name:   _val('live-pair-name') || 'live',
    fwd_primer:    fwd,
    rev_primer:    rev,
    genome_string: _state.genomeString,
    mismatches:    parseInt(_val('live-mm-slider'))  || 3,
    na_conc:       parseFloat(_val('live-na-slider')) || 50,
    mg_conc:       parseFloat(_val('live-mg-slider')) || 0,
    min_size:      50,
    max_size:      3000,
    run_hairpin:   false,
    run_dimer:     false,
  };

  console.log('[live-rerun] started, payload:', payload);

  const btn = document.getElementById('btn-live-run');
  btn.disabled = true;
  _status('Re-running with new parameters…', 'info', true);
  _setElHtml('live-result-summary', '');

  try {
    const res = await fetch('/api/params', {
      method:  'POST',
      headers: { 'Content-Type': 'application/json' },
      body:    JSON.stringify(payload),
    });

    // Parse JSON first — error body is also JSON from our envelope
    let data;
    try { data = await res.json(); }
    catch { throw new Error(`Server returned ${res.status} with non-JSON body`); }

    console.log('[live-rerun] response:', data);

    if (!res.ok || data.success === false) {
      // Clean error message from our envelope — never dump raw FastAPI validation JSON
      const msg = data?.error || data?.detail?.[0]?.msg || `Server error ${res.status}`;
      throw new Error(msg);
    }

    // ── Update live results summary panel ────────────────────────────────
    const { summary, elapsed } = data;
    const ampCount = summary?.amplicons_found ?? 0;
    const topScore = summary?.top_score != null ? summary.top_score.toFixed(1) : '—';
    const topLen   = summary?.top_length   ?? '—';
    const otCount  = summary?.offtargets   ?? 0;

    const summaryHtml = ampCount > 0
      ? `<div class="info-box" style="margin-top:12px;flex-wrap:wrap;gap:12px">
           <div><strong>${ampCount}</strong> amplicon(s) found</div>
           <div>Top score: <strong>${topScore}</strong></div>
           <div>Top length: <strong>${topLen} bp</strong></div>
           <div>Off-targets: <strong>${otCount}</strong></div>
           <div class="text-muted" style="font-size:.7rem">Elapsed: ${elapsed}s</div>
         </div>`
      : `<div class="info-box" style="margin-top:12px;border-color:var(--amber-600);color:var(--amber-600)">
           No amplicons found — try relaxing mismatches or amplicon size limits.
         </div>`;
    _setElHtml('live-result-summary', summaryHtml);

    // ── Update the main dashboard with the new result ─────────────────────
    // Merge the live pair into the result so all tabs refresh consistently
    if (data.results?.primer_pairs?.length) {
      const liveResult = {
        ...(_state.result),
        primer_pairs: data.results.primer_pairs,
        elapsed_seconds: elapsed,
      };
      // Use _loadResult so every component (table, plots, alignment) refreshes
      _loadResult(liveResult);
      // Re-populate live panel inputs (loadResult resets them to the pair values)
      _setVal('live-fwd', fwd);
      _setVal('live-rev', rev);
      _setVal('live-pair-name', payload.primer_name);
      // Restore summary (loadResult clears it)
      _setElHtml('live-result-summary', summaryHtml);
    }

    _status(`Live re-run complete — ${ampCount} amplicon(s) in ${elapsed}s`, 'success');
    _autoDismissStatus(4000);

  } catch (err) {
    console.error('[live-rerun] error:', err);
    _status('Live re-run failed: ' + err.message, 'error');
  } finally {
    btn.disabled = false;
  }
}

// ── Exports ───────────────────────────────────────────────────────────────────
async function _export(fmt) {
  if (!_state.result) { _status('No results to export.', 'error'); return; }
  try {
    const res = await fetch(`/api/export/${fmt}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ result: _state.result, format: fmt }),
    });
    if (!res.ok) throw new Error(await res.text());
    const blob = await res.blob();
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement('a');
    a.href     = url;
    a.download = `pcr_results.${fmt}`;
    a.click();
    URL.revokeObjectURL(url);
  } catch (err) {
    _status('Export failed: ' + err.message, 'error');
  }
}

function _exportHtml() {
  if (!_state.result) { _status('No results to export.', 'error'); return; }
  const data = JSON.stringify(_state.result, null, 2);
  const html = `<!DOCTYPE html><html><head><meta charset="utf-8">
<title>PCR Results</title>
<style>body{font-family:monospace;padding:20px;background:#f9fafb}pre{background:#fff;border:1px solid #e5e7eb;border-radius:8px;padding:16px;overflow-x:auto}</style>
</head><body><h2>In-Silico PCR Results</h2><pre>${data.replace(/</g,'&lt;')}</pre></body></html>`;
  const blob = new Blob([html], { type: 'text/html' });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href     = url; a.download = 'pcr_results.html'; a.click();
  URL.revokeObjectURL(url);
}

// ── Utility ───────────────────────────────────────────────────────────────────
function _setEl(id, val) {
  const el = document.getElementById(id);
  if (el) el.textContent = val ?? '—';
}

function _setElHtml(id, html) {
  const el = document.getElementById(id);
  if (el) el.innerHTML = html;
}

function _setVal(id, val) {
  const el = document.getElementById(id);
  if (el) el.value = val;
}

function _val(id) {
  return document.getElementById(id)?.value || '';
}

function scoreColor(score) {
  if (score >= 70) return '#1e8449';
  if (score >= 40) return '#c4730a';
  return '#c0392b';
}
