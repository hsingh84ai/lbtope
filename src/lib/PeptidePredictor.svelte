<script>
  import { onMount } from 'svelte';

  const EXAMPLE_FASTA = `>pep1
STAIHADQLTPAWRIYSTGNN
>pep2
FQQFGRDVSDFTDSVRDPKT
>pep3
IFGFNALVDR
>pep4
IFGFNALVDR
>pep5
LLFNKVTLADAGFMKQYGECLGDINA`;

  const VARIABLE_MODELS = [
    { file: 'model-final-dpc-flex-svm.bin', label: 'LBtope_Variable', desc: 'original dataset' },
    { file: 'model-variable-non-redundant.bin', label: 'LBtope_Variable_non_redundant', desc: 'non-redundant dataset' },
    { file: 'model-final-dpc-flex-more2-svm.bin', label: 'LBtope_Confirm', desc: 'confirm dataset' },
  ];

  let selectedModelFile = $state('model-final-dpc-flex-svm.bin');
  let probabilityThreshold = $state(60);
  let fastaInput = $state('');

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
  let results = $state([]);

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
    results = [];
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

  function parseFasta(text) {
    const sequences = [];
    let header = 'Unnamed';
    let seqLines = [];
    for (const raw of text.split('\n')) {
      const line = raw.trim();
      if (line.startsWith('>')) {
        if (seqLines.length) { sequences.push({ header, sequence: seqLines.join('') }); seqLines = []; }
        header = line.substring(1).trim() || 'Unnamed';
      } else if (line) {
        seqLines.push(line.replace(/\s/g, '').toUpperCase());
      }
    }
    if (seqLines.length) sequences.push({ header, sequence: seqLines.join('') });
    return sequences;
  }

  function pasteExample() { fastaInput = EXAMPLE_FASTA; }

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
    if (!modelLoaded || !fastaInput.trim()) return;
    predicting = true; results = [];
    await new Promise(r => setTimeout(r, 10));

    const t0 = performance.now();
    const peptides = parseFasta(fastaInput);
    const batchOutput = [];

    for (const pep of peptides) {
      const feats = extractDPCFeatures(pep.sequence);
      const score = svmScore(feats);
      const prob = scoreToProbability(score);

      batchOutput.push({
        id: pep.header,
        sequence: pep.sequence,
        length: pep.sequence.length,
        score,
        probability: prob,
        isEpitope: score >= EPITOPE_THRESHOLD
      });
    }

    results = batchOutput;
    predictMs = Math.round(performance.now() - t0);
    predicting = false;
  }

  function downloadResults() {
    const lines = ['ID\tSequence\tLength\tScore\tProbability\tPrediction\n'];
    for (const r of results) {
      lines.push(`${r.id}\t${r.sequence}\t${r.length}\t${r.score.toFixed(4)}\t${r.probability.toFixed(1)}%\t${r.isEpitope ? 'Epitope' : 'Non-epitope'}\n`);
    }
    const blob = new Blob(lines, { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'lbtope_multiple_peptides.txt';
    a.click(); URL.revokeObjectURL(url);
  }
</script>

<section class="predictor">
  <h2>Multiple Peptides Submission</h2>
  <p class="intro">This module allows users to predict B-cell epitopes in peptides of various lengths simultaneously. Submit multiple peptides in FASTA format.</p>

  <div class="panel">
    <div class="field-row">
      <label class="field-label" for="fastaInput">Paste peptide sequences in FASTA format</label>
      <button class="example-btn" onclick={pasteExample}>Example peptide sequences</button>
    </div>
    <textarea id="fastaInput" rows="6" bind:value={fastaInput} placeholder=">pep1&#10;STAIHADQLTPAWRIYSTGNN"></textarea>
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

  <div class="panel inline-panel">
    <label class="field-label" for="probThreshold">Highlight predictions with probability more than:</label>
    <select id="probThreshold" bind:value={probabilityThreshold}>
      <option value={20}>20%</option>
      <option value={40}>40%</option>
      <option value={60}>60%</option>
      <option value={80}>80%</option>
    </select>
  </div>

  <div class="actions">
    <button class="predict-btn" onclick={predict} disabled={!modelLoaded || predicting || !fastaInput.trim()}>
      {predicting ? 'Scoring peptides...' : 'Submit peptides for prediction'}
    </button>
  </div>

  {#if results.length > 0}
    <div class="results-section">
      <div class="results-header">
        <h2>Batch Prediction Results ({results.length} peptides, {predictMs}ms)</h2>
        <button class="download-btn" onclick={downloadResults}>Download as TSV</button>
      </div>

      <div class="table-wrap">
        <table>
          <thead>
            <tr>
              <th>Peptide ID</th>
              <th>Sequence</th>
              <th>Length</th>
              <th>Score</th>
              <th>Probability</th>
              <th>Prediction</th>
            </tr>
          </thead>
          <tbody>
            {#each results as r}
              <tr class:highlighted={r.probability > probabilityThreshold}>
                <td class="font-bold">{r.id}</td>
                <td class="mono">{r.sequence}</td>
                <td>{r.length}</td>
                <td class="mono">{r.score.toFixed(4)}</td>
                <td class="mono font-bold">{r.probability.toFixed(1)}%</td>
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
    </div>
  {/if}
</section>

<style>
  .predictor { max-width: 840px; margin: 0 auto; padding: 3rem 2rem 5rem; width: 100%; }
  h2 { font-family: var(--font-display); font-weight: 600; font-size: 1.75rem; margin: 0 0 0.5rem; }
  .intro { color: var(--color-quiet); font-size: 0.95rem; line-height: 1.5; margin: 0 0 1.5rem; }
  .panel { background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 10px; padding: 1.25rem 1.5rem; margin-bottom: 1.25rem; }
  .field-row { display: flex; align-items: center; justify-content: space-between; margin-bottom: 0.6rem; flex-wrap: wrap; gap: 0.5rem; }
  .field-label { font-weight: 500; font-size: 0.9rem; }
  .example-btn { font-family: var(--font-body); font-size: 0.8rem; background: var(--color-paper); border: 1px solid var(--color-hairline); border-radius: 6px; padding: 0.35rem 0.7rem; cursor: pointer; }
  .example-btn:hover { border-color: var(--color-structural); }
  textarea { width: 100%; border: 1px solid var(--color-hairline); border-radius: 8px; padding: 0.75rem; font-family: var(--font-mono); font-size: 0.85rem; resize: vertical; }
  select { border: 1px solid var(--color-hairline); border-radius: 6px; padding: 0.4rem 0.6rem; font-family: var(--font-body); font-size: 0.9rem; background: var(--color-paper-raised); }
  .model-panel fieldset { border: 1px solid var(--color-hairline); border-radius: 8px; padding: 1rem; margin: 0.75rem 0 0; }
  .model-panel legend { font-weight: 500; font-size: 0.85rem; color: var(--color-structural-dim); padding: 0 0.4rem; }
  .radio-label { display: block; font-size: 0.9rem; padding: 0.3rem 0; cursor: pointer; }
  .model-desc { color: var(--color-quiet); font-size: 0.8rem; }
  .inline-panel { display: flex; align-items: center; gap: 0.75rem; flex-wrap: wrap; }
  .status { font-family: var(--font-mono); font-size: 0.8rem; color: var(--color-quiet); margin: 0.75rem 0 0; }
  .status.ready { color: var(--color-structural); }
  .error { font-family: var(--font-mono); font-size: 0.8rem; color: #b3261e; margin: 0.75rem 0 0; }
  .actions { margin-bottom: 2rem; }
  .predict-btn { font-family: var(--font-body); font-weight: 500; font-size: 0.95rem; background: var(--color-structural); color: var(--color-paper-raised); border: none; border-radius: 8px; padding: 0.8rem 1.4rem; cursor: pointer; }
  .predict-btn:disabled { background: var(--color-hairline); color: var(--color-quiet); cursor: not-allowed; }
  .results-section { margin-top: 1rem; }
  .results-header { display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 1rem; margin-bottom: 1rem; }
  .results-header h2 { font-size: 1.4rem; color: var(--color-structural); margin: 0; }
  .download-btn { font-family: var(--font-body); font-size: 0.85rem; padding: 0.5rem 1rem; background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 6px; cursor: pointer; }
  .download-btn:hover { border-color: var(--color-structural); }
  .table-wrap { overflow-x: auto; border: 1px solid var(--color-hairline); border-radius: 8px; background: var(--color-paper-raised); }
  table { width: 100%; border-collapse: collapse; font-size: 0.88rem; text-align: left; }
  th { font-weight: 600; color: var(--color-quiet); border-bottom: 2px solid var(--color-hairline); padding: 0.75rem 1rem; background: rgba(0,0,0,0.02); }
  td { padding: 0.75rem 1rem; border-bottom: 1px solid var(--color-hairline); }
  tr:last-child td { border-bottom: none; }
  tr.highlighted { background: rgba(0, 136, 0, 0.06); }
  .mono { font-family: var(--font-mono); letter-spacing: 0.02em; }
  .font-bold { font-weight: 600; }
  .badge { display: inline-block; font-size: 0.75rem; font-weight: 500; padding: 0.25rem 0.6rem; border-radius: 999px; background: #eef1f0; color: var(--color-quiet); }
  .badge.epitope { background: var(--color-epitope-bg, #e8f5e9); color: var(--color-epitope, #2e7d32); }
  @media (max-width: 640px) { .predictor { padding: 2rem 1.25rem 3rem; } }
</style>