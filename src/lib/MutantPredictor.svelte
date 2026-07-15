<script>
  import { onMount } from 'svelte';
  // Accept the passed sequence from App.svelte
  let { initialPeptide = '' } = $props();

  const EXAMPLE_PEPTIDE = "LLFNKVTLADAGFMKQYGECLGDINA";

  const VARIABLE_MODELS = [
    { file: 'model-final-dpc-flex-svm.bin', label: 'LBtope_Variable', desc: 'original dataset' },
    { file: 'model-variable-non-redundant.bin', label: 'LBtope_Variable_non_redundant', desc: 'non-redundant dataset' },
    { file: 'model-final-dpc-flex-more2-svm.bin', label: 'LBtope_Confirm', desc: 'confirm dataset' },
  ];

  let selectedModelFile = $state('model-final-dpc-flex-svm.bin');
  let peptideInput = $state(EXAMPLE_PEPTIDE);
  let viewMode = $state('grid'); // 'grid' | 'ranked'

  // Tracker prevents Svelte 5 reactivity loops from calling predict() more than once
  let lastAutoPredictedPeptide = '';

  // Svelte 5 Effect: When initialPeptide arrives from Antigen Sequence, auto-populate & run!
  $effect(() => {
    const targetSeq = initialPeptide ? initialPeptide.trim() : '';
    // Only trigger if we have a sequence AND we haven't already auto-predicted for it!
    if (targetSeq && targetSeq !== lastAutoPredictedPeptide) {
      if (modelLoaded && !predicting) {
        lastAutoPredictedPeptide = targetSeq;
        peptideInput = targetSeq;
        predict();
      }
    }
  });

  const FEATURE_DIM = 400;
  const AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY';
  const EPITOPE_THRESHOLD = -0.3;

  let gamma = 0, bias = 0, numSV = 0;
  let svAlpha, svNormSq, svFeatOffset, svFeatIdx, svFeatVal;
  let modelLoaded = $state(false);
  let modelLoading = $state(false);
  let modelStatus = $state('');
  let modelError = $state('');

  let predicting = $state(false);
  let predictMs = $state(0);
  let wildTypeData = $state(null);
  let mutationMatrix = $state([]);
  let rankedMutations = $state([]);
  let tooltipData = $state(null);
  let tooltipX = $state(0);
  let tooltipY = $state(0);

  function modelUrl(file) {
    return `${import.meta.env.BASE_URL}models/${file}`;
  }

  function scoreToProbability(score) {
    return Math.max(0, Math.min(100, (score + 1) * 50));
  }

  function parseModelBuffer(buf) {
    const view = new DataView(buf);
    let offset = 0;
    const kernelType = view.getUint32(offset, true); offset += 4;
    gamma = view.getFloat64(offset, true); offset += 8;
    bias = view.getFloat64(offset, true); offset += 8;
    numSV = view.getUint32(offset, true); offset += 4;

    const alphaArr = new Float64Array(numSV);
    const normSqArr = new Float64Array(numSV);
    const offsetArr = new Uint32Array(numSV + 1);
    const idxParts = new Array(numSV);
    const valParts = new Array(numSV);
    let totalNNZ = 0;

    for (let i = 0; i < numSV; i++) {
      const alphaY = view.getFloat64(offset, true); offset += 8;
      const nf = view.getUint16(offset, true); offset += 2;
      const idxs = new Uint16Array(nf);
      const vals = new Float32Array(nf);
      let normSq = 0;
      for (let f = 0; f < nf; f++) {
        const idx = view.getUint16(offset, true); offset += 2;
        const val = view.getFloat64(offset, true); offset += 8;
        idxs[f] = idx; vals[f] = val;
        normSq += val * val;
      }
      alphaArr[i] = alphaY; normSqArr[i] = normSq;
      offsetArr[i] = totalNNZ;
      idxParts[i] = idxs; valParts[i] = vals;
      totalNNZ += nf;
    }
    offsetArr[numSV] = totalNNZ;
    const flatIdx = new Uint16Array(totalNNZ);
    const flatVal = new Float32Array(totalNNZ);
    let p = 0;
    for (let i = 0; i < numSV; i++) {
      flatIdx.set(idxParts[i], p); flatVal.set(valParts[i], p);
      p += idxParts[i].length;
    }
    svAlpha = alphaArr; svNormSq = normSqArr; svFeatOffset = offsetArr;
    svFeatIdx = flatIdx; svFeatVal = flatVal;
    return totalNNZ;
  }

  async function loadModel() {
    modelLoaded = false; modelError = '';
    wildTypeData = null; mutationMatrix = []; rankedMutations = [];
    lastAutoPredictedPeptide = ''; // Reset tracker so incoming peptides re-run on model change
    modelLoading = true;
    modelStatus = `Fetching ${selectedModelFile}...`;
    try {
      const res = await fetch(modelUrl(selectedModelFile));
      if (!res.ok) throw new Error(`${res.status} ${res.statusText} — is ${selectedModelFile} in public/models/?`);
      const buf = await res.arrayBuffer();
      modelStatus = `Parsing ${(buf.byteLength / 1e6).toFixed(1)} MB...`;
      await new Promise(r => setTimeout(r, 10));
      parseModelBuffer(buf);
      modelLoaded = true;
      modelStatus = `${numSV.toLocaleString()} support vectors loaded.`;
    } catch (err) {
      modelError = err.message || String(err);
      modelStatus = '';
    } finally { modelLoading = false; }
  }

  onMount(loadModel);
  function onModelChange() { loadModel(); }

  function cleanSequence(text) {
    const lines = text.split('\n').map(l => l.trim()).filter(Boolean);
    const seqLines = lines.filter(l => !l.startsWith('>'));
    return seqLines.join('').replace(/[^A-Z]/gi, '').toUpperCase();
  }

  function pasteExample() { peptideInput = EXAMPLE_PEPTIDE; }

  function extractDPCFeatures(sequence) {
    const seq = sequence.toUpperCase();
    const len = seq.length - 1;
    const features = new Array(400).fill(0);

    if (len <= 0) return features;

    const getIndex = (char) => {
      const idx = AA_ALPHABET.indexOf(char);
      return idx !== -1 ? idx : 19; // Map 'X'/unknown to Tyrosine ('Y') matching legacy Perl
    };

    for (let i = 0; i < len; i++) {
      const idx1 = getIndex(seq[i]);
      const idx2 = getIndex(seq[i + 1]);
      features[idx1 * 20 + idx2]++;
    }

    for (let i = 0; i < 400; i++) {
      if (features[i] > 0) {
        const perc = (features[i] * 100) / len;
        features[i] = Math.round(perc * 1000) / 1000;
      }
    }
    return features;
  }

  function svmScore(features) {
    let normX = 0;
    for (let i = 0; i < FEATURE_DIM; i++) normX += features[i] * features[i];
    let score = -bias;
    for (let s = 0; s < numSV; s++) {
      const start = svFeatOffset[s], end = svFeatOffset[s + 1];
      let dot = 0;
      for (let f = start; f < end; f++) dot += features[svFeatIdx[f] - 1] * svFeatVal[f];
      score += svAlpha[s] * Math.exp(-gamma * Math.max(0, normX - 2 * dot + svNormSq[s]));
    }
    return score;
  }

  async function predict() {
    const seq = cleanSequence(peptideInput);
    if (!modelLoaded || !seq) return;
    predicting = true; wildTypeData = null; mutationMatrix = []; rankedMutations = [];
    await new Promise(r => setTimeout(r, 10));

    const t0 = performance.now();
    
    // 1. Score Wild-Type
    const wtFeats = extractDPCFeatures(seq);
    const wtScore = svmScore(wtFeats);
    const wtProb = scoreToProbability(wtScore);
    wildTypeData = {
      sequence: seq,
      length: seq.length,
      score: wtScore,
      probability: wtProb,
      isEpitope: wtScore >= EPITOPE_THRESHOLD
    };

    // 2. Generate and score all L * 20 single-point mutants
    const matrix = [];
    const flatList = [];

    for (let i = 0; i < seq.length; i++) {
      const origAA = seq[i];
      const row = { position: i + 1, origAA, cells: {} };

      for (const aa of AA_ALPHABET) {
        if (aa === origAA) {
          row.cells[aa] = {
            aa, position: i + 1, origAA, sequence: seq,
            score: wtScore, probability: wtProb, deltaScore: 0, deltaProb: 0,
            isEpitope: wildTypeData.isEpitope, isWildType: true
          };
        } else {
          const mutSeq = seq.substring(0, i) + aa + seq.substring(i + 1);
          const feats = extractDPCFeatures(mutSeq);
          const score = svmScore(feats);
          const prob = scoreToProbability(score);
          const deltaScore = score - wtScore;
          const deltaProb = prob - wtProb;

          const mutData = {
            aa, position: i + 1, origAA, sequence: mutSeq,
            score, probability: prob, deltaScore, deltaProb,
            isEpitope: score >= EPITOPE_THRESHOLD, isWildType: false
          };
          row.cells[aa] = mutData;
          flatList.push(mutData);
        }
      }
      matrix.push(row);
    }

    // Sort ranked list by highest positive change in probability
    flatList.sort((a, b) => b.deltaProb - a.deltaProb);

    mutationMatrix = matrix;
    rankedMutations = flatList;
    predictMs = Math.round(performance.now() - t0);
    predicting = false;
  }

  function getCellColor(cell) {
    if (cell.isWildType) return 'var(--color-hairline)';
    if (cell.deltaProb > 15) return '#008800';      // strong green enhancement
    if (cell.deltaProb > 5) return '#4caf50';       // mild green enhancement
    if (cell.deltaProb < -15) return '#cc0000';     // strong red degradation
    if (cell.deltaProb < -5) return '#e65100';      // mild orange degradation
    return 'transparent';
  }

  function getCellTextColor(cell) {
    if (cell.isWildType) return 'var(--color-quiet)';
    if (Math.abs(cell.deltaProb) > 5) return '#ffffff';
    return 'var(--color-ink)';
  }

  function onCellEnter(e, cell) {
    const rect = e.target.getBoundingClientRect();
    tooltipX = rect.left + rect.width / 2;
    tooltipY = rect.top - 8;
    tooltipData = cell;
  }
  function onCellLeave() { tooltipData = null; }

  function downloadResults() {
    const lines = ['Position\tWT_AA\tMut_AA\tSequence\tScore\tDelta_Score\tProbability\tDelta_Prob\tPrediction\n'];
    lines.push(`0\t-\t-\t${wildTypeData.sequence}\t${wildTypeData.score.toFixed(4)}\t0.0000\t${wildTypeData.probability.toFixed(1)}%\t0.0%\t${wildTypeData.isEpitope ? 'Epitope' : 'Non-epitope'}\n`);
    
    for (const r of rankedMutations) {
      lines.push(`${r.position}\t${r.origAA}\t${r.aa}\t${r.sequence}\t${r.score.toFixed(4)}\t${r.deltaScore >= 0 ? '+' : ''}${r.deltaScore.toFixed(4)}\t${r.probability.toFixed(1)}%\t${r.deltaProb >= 0 ? '+' : ''}${r.deltaProb.toFixed(1)}%\t${r.isEpitope ? 'Epitope' : 'Non-epitope'}\n`);
    }
    const blob = new Blob(lines, { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'lbtope_peptide_mutants.txt';
    a.click(); URL.revokeObjectURL(url);
  }
</script>

<section class="predictor">
  <h2>Peptide Mutants Prediction</h2>
  <p class="intro">This module performs a systematic point-mutation scan across your peptide. It substitutes every position with all 20 standard amino acids to identify targeted mutations that enhance B-cell epitope affinity.</p>

  <div class="panel">
    <div class="field-row">
      <label class="field-label" for="peptideInput">Paste single peptide sequence</label>
      <button class="example-btn" onclick={pasteExample}>Example peptide sequence</button>
    </div>
    <input type="text" id="peptideInput" class="text-input" bind:value={peptideInput} placeholder="e.g. STAIHADQLTPAWRIYSTGNN" />
  </div>

  <div class="panel model-panel">
    <p class="field-label">Select variable-length model for prediction</p>
    <fieldset>
      <legend>Models developed on variable length epitopes</legend>
      {#each VARIABLE_MODELS as m}
        <label class="radio-label">
          <input type="radio" bind:group={selectedModelFile} value={m.file} onchange={onModelChange} />
          {m.label} <span class="model-desc">({m.desc})</span>
        </label>
      {/each}
    </fieldset>
    {#if modelError}
      <p class="error">Error: {modelError}</p>
    {:else if modelStatus}
      <p class="status" class:ready={modelLoaded}>{modelStatus}</p>
    {/if}
  </div>

  <div class="actions">
    <button class="predict-btn" onclick={predict} disabled={!modelLoaded || predicting || !cleanSequence(peptideInput)}>
      {predicting ? 'Scanning mutations...' : 'Generate mutants & predict epitopes'}
    </button>
  </div>

  {#if wildTypeData}
    <div class="results-section">
      <div class="results-header">
        <div>
          <h2>Mutation Scan Results ({mutationMatrix.length * 20} evaluated, {predictMs}ms)</h2>
          <p class="wt-summary">
            Wild-Type: <span class="mono font-bold">{wildTypeData.sequence}</span> | 
            Score: <span class="mono">{wildTypeData.score.toFixed(4)}</span> | 
            Probability: <span class="mono font-bold">{wildTypeData.probability.toFixed(1)}%</span> 
            <span class="badge" class:epitope={wildTypeData.isEpitope}>
              {wildTypeData.isEpitope ? 'Epitope' : 'Non-epitope'}
            </span>
          </p>
        </div>
        <div class="header-controls">
          <div class="toggle-group">
            <button class:active={viewMode === 'grid'} onclick={() => viewMode = 'grid'}>Heatmap Grid</button>
            <button class:active={viewMode === 'ranked'} onclick={() => viewMode = 'ranked'}>Top Enhancements</button>
          </div>
          <button class="download-btn" onclick={downloadResults}>Download TSV</button>
        </div>
      </div>

      {#if viewMode === 'grid'}
        <div class="legend-bar">
          <span class="legend-item"><span class="box" style="background:#cc0000"></span> &lt; -15% Loss</span>
          <span class="legend-item"><span class="box" style="background:#e65100"></span> -5% to -15%</span>
          <span class="legend-item"><span class="box" style="background:var(--color-hairline)"></span> Neutral / Wild-Type</span>
          <span class="legend-item"><span class="box" style="background:#4caf50"></span> +5% to +15%</span>
          <span class="legend-item"><span class="box" style="background:#008800"></span> &gt; +15% Gain</span>
        </div>

        <div class="table-wrap">
          <table class="grid-table">
            <thead>
              <tr>
                <th class="sticky-col">Pos : WT</th>
                {#each AA_ALPHABET as aa}
                  <th class="text-center">{aa}</th>
                {/each}
              </tr>
            </thead>
            <tbody>
              {#each mutationMatrix as row}
                <tr>
                  <td class="sticky-col font-bold">{row.position} : {row.origAA}</td>
                  {#each AA_ALPHABET as aa}
                    {@const cell = row.cells[aa]}
                    <td 
                      class="cell text-center mono"
                      style="background: {getCellColor(cell)}; color: {getCellTextColor(cell)}"
                      onmouseenter={(e) => onCellEnter(e, cell)}
                      onmouseleave={onCellLeave}
                    >
                      {cell.isWildType ? '•' : `${cell.deltaProb >= 0 ? '+' : ''}${Math.round(cell.deltaProb)}`}
                    </td>
                  {/each}
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      {:else}
        <div class="table-wrap">
          <table>
            <thead>
              <tr>
                <th>Rank</th>
                <th>Mutation</th>
                <th>Position</th>
                <th>Mutant Sequence</th>
                <th>Score</th>
                <th>Probability</th>
                <th>Δ Probability</th>
                <th>Prediction</th>
              </tr>
            </thead>
            <tbody>
              {#each rankedMutations.slice(0, 50) as r, index}
                <tr class:positive={r.deltaProb > 0}>
                  <td>#{index + 1}</td>
                  <td class="font-bold mono">{r.origAA} → {r.aa}</td>
                  <td>{r.position}</td>
                  <td class="mono">{r.sequence}</td>
                  <td class="mono">{r.score.toFixed(4)}</td>
                  <td class="mono font-bold">{r.probability.toFixed(1)}%</td>
                  <td class="mono font-bold" style="color: {r.deltaProb > 0 ? '#008800' : r.deltaProb < 0 ? '#cc0000' : 'inherit'}">
                    {r.deltaProb >= 0 ? '+' : ''}{r.deltaProb.toFixed(1)}%
                  </td>
                  <td>
                    <span class="badge" class:epitope={r.isEpitope}>
                      {r.isEpitope ? 'Epitope' : 'Non-epitope'}
                    </span>
                  </td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      {/if}
    </div>
  {/if}
</section>

{#if tooltipData}
  <div class="tooltip" style="left:{tooltipX}px; top:{tooltipY}px">
    {#if tooltipData.isWildType}
      <span>Pos {tooltipData.position} ({tooltipData.origAA}) | Wild-Type | Prob: {tooltipData.probability.toFixed(1)}%</span>
    {:else}
      <span>Pos {tooltipData.position}: {tooltipData.origAA} → {tooltipData.aa} | Prob: {tooltipData.probability.toFixed(1)}% ({tooltipData.deltaProb >= 0 ? '+' : ''}{tooltipData.deltaProb.toFixed(1)}%)</span>
    {/if}
  </div>
{/if}

<style>
  .predictor { max-width: 880px; margin: 0 auto; padding: 3rem 2rem 5rem; width: 100%; }
  h2 { font-family: var(--font-display); font-weight: 600; font-size: 1.75rem; margin: 0 0 0.5rem; }
  .intro { color: var(--color-quiet); font-size: 0.95rem; line-height: 1.5; margin: 0 0 1.5rem; }
  .panel { background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 10px; padding: 1.25rem 1.5rem; margin-bottom: 1.25rem; }
  .field-row { display: flex; align-items: center; justify-content: space-between; margin-bottom: 0.6rem; flex-wrap: wrap; gap: 0.5rem; }
  .field-label { font-weight: 500; font-size: 0.9rem; }
  .example-btn { font-family: var(--font-body); font-size: 0.8rem; background: var(--color-paper); border: 1px solid var(--color-hairline); border-radius: 6px; padding: 0.35rem 0.7rem; cursor: pointer; }
  .example-btn:hover { border-color: var(--color-structural); }
  .text-input { width: 100%; border: 1px solid var(--color-hairline); border-radius: 8px; padding: 0.75rem; font-family: var(--font-mono); font-size: 0.95rem; letter-spacing: 0.05em; text-transform: uppercase; }
  .model-panel fieldset { border: 1px solid var(--color-hairline); border-radius: 8px; padding: 1rem; margin: 0.75rem 0 0; }
  .model-panel legend { font-weight: 500; font-size: 0.85rem; color: var(--color-structural-dim); padding: 0 0.4rem; }
  .radio-label { display: block; font-size: 0.9rem; padding: 0.3rem 0; cursor: pointer; }
  .model-desc { color: var(--color-quiet); font-size: 0.8rem; }
  .status { font-family: var(--font-mono); font-size: 0.8rem; color: var(--color-quiet); margin: 0.75rem 0 0; }
  .status.ready { color: var(--color-structural); }
  .error { font-family: var(--font-mono); font-size: 0.8rem; color: #b3261e; margin: 0.75rem 0 0; }
  .actions { margin-bottom: 2rem; }
  .predict-btn { font-family: var(--font-body); font-weight: 500; font-size: 0.95rem; background: var(--color-structural); color: var(--color-paper-raised); border: none; border-radius: 8px; padding: 0.8rem 1.4rem; cursor: pointer; }
  .predict-btn:disabled { background: var(--color-hairline); color: var(--color-quiet); cursor: not-allowed; }
  .results-section { margin-top: 1.5rem; }
  .results-header { display: flex; align-items: flex-start; justify-content: space-between; flex-wrap: wrap; gap: 1rem; margin-bottom: 1rem; }
  .results-header h2 { font-size: 1.4rem; color: var(--color-structural); margin: 0 0 0.3rem; }
  .wt-summary { font-size: 0.9rem; color: var(--color-ink); margin: 0; }
  .header-controls { display: flex; align-items: center; gap: 0.75rem; flex-wrap: wrap; }
  .toggle-group { display: inline-flex; border: 1px solid var(--color-hairline); border-radius: 6px; overflow: hidden; background: var(--color-paper-raised); }
  .toggle-group button { border: none; background: none; font-family: var(--font-body); font-size: 0.8rem; padding: 0.45rem 0.8rem; cursor: pointer; color: var(--color-quiet); transition: all 0.15s; }
  .toggle-group button.active { background: var(--color-structural); color: #fff; font-weight: 500; }
  .download-btn { font-family: var(--font-body); font-size: 0.85rem; padding: 0.45rem 0.9rem; background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 6px; cursor: pointer; }
  .download-btn:hover { border-color: var(--color-structural); }
  .legend-bar { display: flex; align-items: center; gap: 1.25rem; font-size: 0.8rem; color: var(--color-quiet); margin-bottom: 0.75rem; flex-wrap: wrap; }
  .legend-item { display: inline-flex; align-items: center; gap: 0.35rem; }
  .legend-item .box { width: 12px; height: 12px; border-radius: 2px; display: inline-block; }
  .table-wrap { overflow-x: auto; border: 1px solid var(--color-hairline); border-radius: 8px; background: var(--color-paper-raised); max-height: 600px; }
  table { width: 100%; border-collapse: collapse; font-size: 0.85rem; text-align: left; }
  th { font-weight: 600; color: var(--color-quiet); border-bottom: 2px solid var(--color-hairline); padding: 0.6rem 0.5rem; background: #f8f9fa; position: sticky; top: 0; z-index: 2; }
  td { padding: 0.5rem; border-bottom: 1px solid var(--color-hairline); }
  tr:last-child td { border-bottom: none; }
  tr.positive { background: rgba(0, 136, 0, 0.04); }
  .grid-table th, .grid-table td { padding: 0.4rem 0.25rem; }
  .sticky-col { position: sticky; left: 0; background: var(--color-paper-raised); z-index: 1; border-right: 1px solid var(--color-hairline); padding-left: 0.75rem !important; padding-right: 0.75rem !important; white-space: nowrap; }
  th.sticky-col { background: #f8f9fa; z-index: 3; }
  .cell { cursor: default; transition: opacity 0.1s; font-size: 0.78rem; border-radius: 2px; }
  .cell:hover { opacity: 0.8; outline: 1px solid var(--color-ink); }
  .text-center { text-align: center; }
  .mono { font-family: var(--font-mono); letter-spacing: 0.02em; }
  .font-bold { font-weight: 600; }
  .badge { display: inline-block; font-size: 0.72rem; font-weight: 500; padding: 0.2rem 0.5rem; border-radius: 999px; background: #eef1f0; color: var(--color-quiet); }
  .badge.epitope { background: var(--color-epitope-bg, #e8f5e9); color: var(--color-epitope, #2e7d32); }
  .tooltip { position: fixed; transform: translate(-50%, -100%); background: #222; color: #fff; font-family: var(--font-mono); font-size: 0.78rem; font-weight: 500; padding: 0.4rem 0.75rem; border-radius: 4px; white-space: nowrap; pointer-events: none; z-index: 1000; box-shadow: 0 2px 8px rgba(0,0,0,0.25); }
  @media (max-width: 640px) { .predictor { padding: 2rem 1.25rem 3rem; } .header-controls { width: 100%; justify-content: space-between; } }
</style>