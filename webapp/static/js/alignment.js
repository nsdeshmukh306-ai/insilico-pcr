/**
 * alignment.js — Render nucleotide alignment viewer from PCR binding site data.
 *
 * Exposes:
 *   renderAlignment(containerId, fwdSite, revSite, ampliconSeq)
 */

function _complement(base) {
  return { A:'T', T:'A', C:'G', G:'C' }[base.toUpperCase()] || 'N';
}

function _span(text, cls) {
  return `<span class="${cls}">${text}</span>`;
}

/**
 * Build a coloured alignment row from primer and template strings.
 * mismatch_positions: array of 0-based primer positions that mismatch.
 */
function _buildAlignmentRows(label, primerSeq, templateSeq, mismatchPositions) {
  const mmSet = new Set(mismatchPositions || []);
  let primerHtml = '';
  let matchHtml  = '';
  let tmplHtml   = '';

  const n = Math.max(primerSeq.length, templateSeq.length);

  for (let i = 0; i < n; i++) {
    const p = primerSeq[i]   || '-';
    const t = templateSeq[i] || '-';
    const isMM = mmSet.has(i);

    primerHtml += isMM
      ? _span(p, 'aln-mm')
      : _span(p, 'aln-primer');

    matchHtml += isMM
      ? _span('✗', 'aln-mm')
      : _span('|', 'aln-match');

    tmplHtml += isMM
      ? _span(t, 'aln-mm')
      : _span(t, 'aln-match');
  }

  return { primerHtml, matchHtml, tmplHtml };
}

/**
 * Build position ruler string.
 */
function _ruler(start, length, step = 10) {
  let ruler = ' '.repeat(start % step === 0 ? 0 : (step - start % step));
  let pos   = Math.ceil(start / step) * step;
  while (pos - start < length) {
    const label = String(pos);
    ruler += label;
    pos += step;
    // pad to next tick
    const gap = step - label.length;
    if (gap > 0) ruler += ' '.repeat(gap);
  }
  return ruler.slice(0, length);
}

/**
 * Render the full alignment block into containerId.
 *
 * @param {string} containerId
 * @param {object} fwdSite  — fwd_binding from JSON
 * @param {object} revSite  — rev_binding from JSON
 * @param {string} ampliconSeq — full amplicon sequence (may be absent)
 */
function renderAlignment(containerId, fwdSite, revSite, ampliconSeq) {
  const el = document.getElementById(containerId);
  if (!el) return;

  if (!fwdSite && !revSite) {
    el.innerHTML = '<div class="empty-state"><div class="empty-state-title">No binding site data</div></div>';
    return;
  }

  let html = '<div class="alignment-viewer">';

  function renderSite(site, label, isRev) {
    if (!site) return '';
    const primer   = site.primer_seq || site.sequence || '';
    const template = site.aligned_template || primer.split('').map(_complement).join('');
    const mmPos    = site.mismatch_positions || [];
    const strand   = site.strand || (isRev ? '-' : '+');
    const start    = site.start || 0;
    const end      = site.end   || start + primer.length;

    const { primerHtml, matchHtml, tmplHtml } = _buildAlignmentRows(
      label, primer, template, mmPos
    );

    const mmCount  = site.mismatch_count || mmPos.length;
    const tm       = site.tm_celsius;
    const tmStr    = (tm && tm > 0 && tm < 120) ? `Tm=${tm.toFixed(1)}°C` : '';
    const swStr    = site.sw_score ? `SW=${site.sw_score.toFixed(0)}` : '';
    const infoStr  = [tmStr, swStr, `${mmCount}mm`, `strand ${strand}`].filter(Boolean).join(' · ');

    const ruler = _ruler(start, primer.length);

    let block = '';
    block += `<div class="aln-pos-row">${' '.repeat(10)}${ruler}</div>`;
    block += `<div class="aln-row"><span class="aln-label">${label} 5'</span><span class="aln-seq">${primerHtml}</span></div>`;
    block += `<div class="aln-row"><span class="aln-label"></span><span class="aln-seq">${matchHtml}</span></div>`;
    block += `<div class="aln-row"><span class="aln-label">Template</span><span class="aln-seq">${tmplHtml}</span></div>`;
    block += `<div class="aln-row mt-2"><span class="aln-label"></span><span class="text-xs text-muted">${infoStr}</span></div>`;
    block += `<div style="height:16px"></div>`;

    return block;
  }

  html += renderSite(fwdSite, 'Forward', false);
  html += renderSite(revSite, 'Reverse', true);

  // Amplicon sequence preview (first+last 40 bp)
  if (ampliconSeq && ampliconSeq.length > 0) {
    const preview = ampliconSeq.length <= 80
      ? ampliconSeq
      : ampliconSeq.slice(0, 40) + '···' + ampliconSeq.slice(-40);
    html += `<div class="aln-row"><span class="aln-label">Amplicon</span>`;
    html += `<span class="aln-seq" style="color:var(--gray-600);letter-spacing:.04em">${preview}</span></div>`;
    html += `<div class="aln-row"><span class="aln-label"></span><span class="text-xs text-muted">${ampliconSeq.length} bp</span></div>`;
  }

  html += '</div>';
  el.innerHTML = html;
}

/**
 * Render a compact inline alignment badge string (for table rows).
 */
function alignmentBadge(mismatchCount, threeEndMM) {
  if (mismatchCount === 0) return '<span class="badge badge-green">Perfect</span>';
  const cls = mismatchCount <= 1 ? 'badge-amber' : 'badge-red';
  let txt = `${mismatchCount}mm`;
  if (threeEndMM) txt += " 3'";
  return `<span class="badge ${cls}">${txt}</span>`;
}
