<script>
  import { onMount } from 'svelte';
  // Accept the callback prop from App.svelte
  let { onAnalyzeMutant = () => {} } = $props();

  const EXAMPLE_FASTA = `>1525
MARFGDEMPARYGGGGGAGAAAGVVVGAAGGRGAGGSRQGGQPGAQRMYKQSMAQRARTMALYNPIPVRQNCLTVNRSLFLFSE
DNVVRKYAKKITEWPPFEYMILATIIANCIVIALEQHLPDDDKTPMSERLDDDTEPYFIGIFCFEAGIKIALGFAFHKGSYLRNG
WNVMDFVVVLTGILATVGTEFDLRTLRAVRVLRPLKLVSGIPSLQVVLKSIMKAMIPLLQIGLLLFFAILIFAIIGLEMFYMG
KFHTTCFEEGTDDIQGESPAPCGTEEARTCPNGTRCQPYWEGPNNGITQFDNILFAVLTVFQCITMEGWTDLLYNSNDASGNT
WNWLYFIPLIIIGSFFMLNLVLGVLSGEFAKERERVENRRAFLKLRRQQQIERELNGYMEWISKAEEVIILAEDEETDVEQRHP
FDGALRRATIKKSSKTDLLHPEEAEDQLADIASVGSPFARASIKSKLENSSFHKKERRRMRFYIRRMVKTQAFYWTVLSLVALN
TLCVAIVHYNQPEWLSDFLYYAEFIIFLGLFMSEMFIKMYGLGTRPYFHSSFNCFDCGVIIGSIFEVIWAVIKPGTSFGISVLR
ALRLLRIFKVTKYWASLRNLVVSLLNSMKSISLILFLLFLFIVVFALLGMQLFGGQFNFDEGTPPPTNFDTFPAAIMTVFQILT
GEDWNEVMYDGIKSSQGGVQGGMVFSIYFIIVLTLFGNYTLLNVFLAIAVDNLANAQELTKDEQEEEEAVNQKLALQKAKEVAE
VSPLSAANMSIAMKEQQKNQKPAKSVWEQRTSSEMRKQNLLASREALYSEMDPEERWKASYARHLRPDMKTHLDRRPLVVDPQE
NRNNNNTNKSRVAEPTVDQRLGQQRAEEDFLRKQARRHHDRARDDPSAHAAAGLDARRPWAGSQEAELSREGPYGRESDDHQAR
EGGLEPPGFWEEGEAERGKAGDPHRRHAHRRQGVGGSGGSSRSGSPPRTGTADGEPRRHRVHRRPGEDGPDDKAERRGRHREGS
RPARSGGEGEAEGPDGGGGGGGGERRRRRHRHGPPPPAYDPDARRDD
>129246
MKKYLLGIGLILALIACKQNVSSLDEKNSVSVDLPGEMKVLVSKEKNKDGKYDLIATVDKLELKGTSDKNNGSGVLEGVKADKS
KVKLTISDDLGQTTLEVFKEDGKTLVSKKVTSKDKSSTEEKFNEKGEVSEKIITRADGTRLEYTGIKSDGSGKAKEVLKGYVL
EGTLTAEKTTLVVKEGTVTLSKNISKSGEVSVELNDTDSSAATKKTTAAWNSGTSTLTITIVNSSKKTKDLVFTKENTTITIQQ
YDSNGTKLEGSAVEITKLDEIKNALK`;

  const FIXED_MODELS = [
    { file: 'model-20mer-complete.bin', label: 'LBtope_Fixed', desc: 'original dataset', window: 20 },
    { file: 'model-20mer-non-variant.bin', label: 'LBtope_Fixed_non_redundant', desc: 'non-redundant dataset', window: 20 },
  ];
  const VARIABLE_MODELS = [
    { file: 'model-final-dpc-flex-svm.bin', label: 'LBtope_Variable', desc: 'original dataset' },
    { file: 'model-variable-non-redundant.bin', label: 'LBtope_Variable_non_redundant', desc: 'non-redundant dataset' },
    { file: 'model-final-dpc-flex-more2-svm.bin', label: 'LBtope_Confirm', desc: 'confirm dataset' },
  ];

  let selectedModelFile = $state('model-final-dpc-flex-svm.bin');
  let epitopeLength = $state(15);
  let probabilityThreshold = $state(60);
  let fastaInput = $state('');

  function getModelConfig(file) {
    const fixed = FIXED_MODELS.find(m => m.file === file);
    if (fixed) return { ...fixed, type: 'fixed' };
    const variable = VARIABLE_MODELS.find(m => m.file === file);
    if (variable) return { ...variable, type: 'variable' };
    return null;
  }

  let currentWindow = $derived.by(() => {
    const cfg = getModelConfig(selectedModelFile);
    if (!cfg) return 15;
    return cfg.type === 'fixed' ? cfg.window : epitopeLength;
  });

  const FEATURE_DIM = 400;
  const ALPHABET = 'ACDEFGHIKLMNPQRSTVWY';
  const EPITOPE_THRESHOLD = -0.3;

  let gamma = 0, bias = 0, numSV = 0;
  let svAlpha, svNormSq, svFeatOffset, svFeatIdx, svFeatVal;
  let modelLoaded = $state(false);
  let modelLoading = $state(false);
  let modelStatus = $state('');
  let modelError = $state('');

  let predicting = $state(false);
  let predictMs = $state(0);
  let allResults = $state([]);
  let scoreStats = $state(null);

  let tooltipData = $state(null);
  let pinnedData = $state(null);
  let tooltipX = $state(0);
  let tooltipY = $state(0);

  function modelUrl(file) {
    return `${import.meta.env.BASE_URL}models/${file}`;
  }

  // ---------------------------------------------------------------------------
  // Score-to-probability: PLACEHOLDER — replace with the real mapping from
  // your original Perl/PHP once you find it.
  //
  // This version uses the score's rank within the observed distribution for
  // THIS sequence, so it's self-calibrating rather than using a made-up
  // sigmoid.  It won't match the original tool's exact numbers, but it will
  // at least produce a reasonable color distribution.
  // ---------------------------------------------------------------------------
  function scoreToProbability(score) {
    // Linear: -1 = 0%, +1 = 100%
    return Math.max(0, Math.min(100, (score + 1) * 50));
  }

  function probabilityToColor(prob) {
      // If the probability is below or equal to the selected threshold, do not color code it
      if (prob <= probabilityThreshold) return '#222222'; // black: non-epitope / below threshold
      
      if (prob > 80) return '#cc0000';      // red: 81-100%
      if (prob > 60) return '#008800';      // green: 61-80%
      if (prob > 40) return '#0066cc';      // blue: 41-60%
      if (prob > 20) return '#cc8800';      // orange: 21-40%
      return '#222222';                     // black: 0-20%
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
        idxs[f] = idx; vals[f] = val; normSq += val * val;
      }
      alphaArr[i] = alphaY; normSqArr[i] = normSq;
      offsetArr[i] = totalNNZ; idxParts[i] = idxs; valParts[i] = vals;
      totalNNZ += nf;
    }
    offsetArr[numSV] = totalNNZ;
    const flatIdx = new Uint16Array(totalNNZ);
    const flatVal = new Float32Array(totalNNZ);
    let p = 0;
    for (let i = 0; i < numSV; i++) {
      flatIdx.set(idxParts[i], p); flatVal.set(valParts[i], p); p += idxParts[i].length;
    }
    svAlpha = alphaArr; svNormSq = normSqArr; svFeatOffset = offsetArr;
    svFeatIdx = flatIdx; svFeatVal = flatVal;
    return totalNNZ;
  }

  async function loadModel() {
    modelLoaded = false; modelError = ''; allResults = []; scoreStats = null;
    modelLoading = true;
    modelStatus = `Fetching ${selectedModelFile}...`;
    try {
      const res = await fetch(modelUrl(selectedModelFile));
      if (!res.ok) throw new Error(`${res.status} ${res.statusText} — is ${selectedModelFile} in public/models/?`);
      const buf = await res.arrayBuffer();
      modelStatus = `Parsing ${(buf.byteLength / 1e6).toFixed(1)} MB...`;
      await new Promise(r => setTimeout(r, 10));
      const totalNNZ = parseModelBuffer(buf);
      modelLoaded = true;
      modelStatus = `${numSV.toLocaleString()} support vectors loaded.`;
    } catch (err) {
      modelError = err.message || String(err);
      modelStatus = '';
    } finally { modelLoading = false; }
  }

  onMount(() => {
    loadModel();
    // Close pinned popup when clicking outside
    window.addEventListener('click', closePinned);
    return () => window.removeEventListener('click', closePinned);
  });
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

  const AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY";

  /**
   * Computes 400 Dipeptide Composition (DPC) features matching legacy pro2dpc Perl behavior.
   * @param {string} sequence - The amino acid window (padded with 'X' if applicable).
   * @returns {number[]} Array of 400 float64 features.
   */
  export function extractDPCFeatures(sequence) {
    const seq = sequence.toUpperCase();
    const len = seq.length - 1;
    const features = new Array(400).fill(0);

    if (len <= 0) return features;

    // Helper: Map residue to 0-19 index. Legacy Perl mapped 'X'/unknown to index -1, 
    // which resolved to index 20 ('Y') in 1-based indexing (index 19 in 0-based JS).
    const getIndex = (char) => {
      const idx = AA_ALPHABET.indexOf(char);
      return idx !== -1 ? idx : 19; // Map 'X' and any invalid residue to 'Y' (19)
    };

    // Count all adjacent pairs
    for (let i = 0; i < len; i++) {
      const idx1 = getIndex(seq[i]);
      const idx2 = getIndex(seq[i + 1]);
      const featureIndex = idx1 * 20 + idx2;
      features[featureIndex]++;
    }

    // Convert to percentages rounded to 3 decimal places to match legacy %5.3f formatting
    for (let i = 0; i < 400; i++) {
      if (features[i] > 0) {
        const perc = (features[i] * 100) / len;
        // Match C/Perl %5.3f float rounding precision exactly
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
    predicting = true; allResults = [];
    await new Promise(r => setTimeout(r, 10));

    const t0 = performance.now();
    const sequences = parseFasta(fastaInput);
    const win = currentWindow;
    
    // Use floor for left and ceil for right so even window lengths (like 20-mers) 
    // always generate exactly 1 window per residue (windowResults.length === seqLen)
    const padLeft = Math.floor((win - 1) / 2);
    const padRight = Math.ceil((win - 1) / 2);
    const output = [];
    let globalMin = Infinity, globalMax = -Infinity;

    for (const seq of sequences) {
      const padded = 'X'.repeat(padLeft) + seq.sequence + 'X'.repeat(padRight);
      const windowResults = [];
      
      for (let i = 0; i <= padded.length - win; i++) {
        const frag = padded.substring(i, i + win);
        const feats = extractDPCFeatures(frag);
        const score = svmScore(feats);
        windowResults.push({
          position: i + 1,
          fragment: frag,
          score,
          probability: 0,
          isEpitope: score >= EPITOPE_THRESHOLD
        });
        if (score < globalMin) globalMin = score;
        if (score > globalMax) globalMax = score;
      }

      // Assign probability: linear map from SVM score [-1, +1] → [0%, 100%]
      for (const r of windowResults) {
        r.probability = scoreToProbability(r.score);
      }

      const seqLen = seq.sequence.length;
      const residueData = [];
      
      // 1-to-1 Centered Mapping: Each residue r gets the score of Window r (centered directly on it)
      for (let r = 0; r < seqLen; r++) {
        const centerWindow = windowResults[r];
        const prob = centerWindow ? centerWindow.probability : 0;

        residueData.push({
          aa: seq.sequence[r],
          position: r + 1,
          probability: prob,
          color: probabilityToColor(prob),
          bestFragment: centerWindow ? centerWindow.fragment : '',
          bestScore: centerWindow ? centerWindow.score : 0,
          bestProbability: prob,
          bestPosition: centerWindow ? centerWindow.position : 0,
        });
      }
      output.push({ header: seq.header, seqLen, residueData, windowResults });
    }

    allResults = output;
    scoreStats = { min: globalMin.toFixed(4), max: globalMax.toFixed(4) };
    predictMs = Math.round(performance.now() - t0);
    predicting = false;
  }

  let leaveTimeout;
  let enterTimeout;

  function onResidueEnter(e, rd) {
    // If a popup is already pinned by a click, don't distract the user with hover tooltips
    if (pinnedData) return; 
    
    const rect = e.target.getBoundingClientRect();
    tooltipX = rect.left + rect.width / 2;
    tooltipY = rect.top - 8;
    tooltipData = rd;
  }

  function onResidueLeave() {
    tooltipData = null;
  }

  // Clicking a residue pins the popup so it stays open permanently until acted upon!
  function onResidueClick(e, rd) {
    e.stopPropagation(); // Prevent immediate closing from window click listeners
    const rect = e.target.getBoundingClientRect();
    tooltipX = rect.left + rect.width / 2;
    tooltipY = rect.top - 8;
    
    // Toggle: if clicking the same residue again, unpin it. Otherwise, pin the new one.
    if (pinnedData && pinnedData.position === rd.position) {
      pinnedData = null;
    } else {
      tooltipData = null; // Clear hover state
      pinnedData = rd;
    }
  }

  // Allow users to close the pinned popup by clicking anywhere outside of it
  function closePinned() {
    pinnedData = null;
  }

  function onTooltipEnter() {
    clearTimeout(leaveTimeout);
    clearTimeout(enterTimeout);
  }

  function onTooltipLeave() {
    leaveTimeout = setTimeout(() => { tooltipData = null; }, 200);
  }

  function downloadResults(seqResult) {
    const lines = [`>${seqResult.header}\tSequence length=${seqResult.seqLen}\n`];
    lines.push('Position\tFragment\tScore\tProbability\tPrediction\n');
    for (const r of seqResult.windowResults) {
      lines.push(`${r.position}\t${r.fragment}\t${r.score.toFixed(4)}\t${r.probability.toFixed(1)}\t${r.isEpitope ? 'Epitope' : 'Non-epitope'}\n`);
    }
    const blob = new Blob(lines, { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = `lbtope_${seqResult.header.replace(/[^a-zA-Z0-9]/g, '_')}.txt`;
    a.click(); URL.revokeObjectURL(url);
  }
</script>

<section class="predictor">
  <h2>Antigen or Protein Submission</h2>
  <p class="intro">This module allows you to predict B-cell epitopes in antigen sequences. Submit one or more protein/antigen sequences in FASTA format.</p>

  <div class="panel">
    <div class="field-row">
      <label class="field-label" for="fastaInput">Paste antigen sequences in FASTA format</label>
      <button class="example-btn" onclick={pasteExample}>Example antigen sequences</button>
    </div>
    <textarea id="fastaInput" rows="6" bind:value={fastaInput} placeholder=">sequence_name&#10;MFTKQVGKRLDGRPFV..."></textarea>
  </div>

  <div class="panel model-panel">
    <p class="field-label">Select model for epitope prediction</p>
    <fieldset>
      <legend>Fixed length epitopes (20-mer)</legend>
      {#each FIXED_MODELS as m}
        <label class="radio-label">
          <input type="radio" bind:group={selectedModelFile} value={m.file} onchange={onModelChange} />
          {m.label} <span class="model-desc">({m.desc})</span>
        </label>
      {/each}
    </fieldset>
    <fieldset>
      <legend>Variable length epitopes</legend>
      {#each VARIABLE_MODELS as m}
        <label class="radio-label">
          <input type="radio" bind:group={selectedModelFile} value={m.file} onchange={onModelChange} />
          {m.label} <span class="model-desc">({m.desc})</span>
        </label>
      {/each}
      <div class="length-picker">
        <label for="epitopeLength">Select length of epitopes:</label>
        <select id="epitopeLength" bind:value={epitopeLength} disabled={getModelConfig(selectedModelFile)?.type === 'fixed'}>
          {#each Array.from({length: 26}, (_, i) => i + 5) as len}
            <option value={len}>{len}</option>
          {/each}
        </select>
      </div>
    </fieldset>
    {#if modelError}
      <p class="error">Error: {modelError}</p>
    {:else if modelStatus}
      <p class="status" class:ready={modelLoaded}>{modelStatus}</p>
    {/if}
  </div>

  <div class="panel inline-panel">
    <label class="field-label" for="probThreshold">Probability threshold for epitope display</label>
    <select id="probThreshold" bind:value={probabilityThreshold}>
      <option value={20}>20%</option>
      <option value={40}>40%</option>
      <option value={60}>60%</option>
      <option value={80}>80%</option>
    </select>
  </div>

  <div class="actions">
    <button class="predict-btn" onclick={predict} disabled={!modelLoaded || predicting || !fastaInput.trim()}>
      {predicting ? 'Scoring windows...' : 'Submit antigen for prediction'}
    </button>
  </div>

  {#if allResults.length > 0}
    <div class="results-section">
      <h2>Your sequences have been processed</h2>



      <p class="legend-label">Probability scale in percentage (0-100%) for each epitope [Above {probabilityThreshold}% motif can be considered as epitope]</p>
      <div class="legend-bar">
        <span class="legend-swatch" style="background:#222222"><small>0-20%</small></span>
        <span class="legend-swatch" style="background:#cc8800"><small>21-40%</small></span>
        <span class="legend-swatch" style="background:#0066cc"><small>41-60%</small></span>
        <span class="legend-swatch" style="background:#008800"><small>61-80%</small></span>
        <span class="legend-swatch" style="background:#cc0000"><small>81-100%</small></span>
      </div>
      <p class="legend-hint">Mouse over the amino acids of a sequence to view the detail result in a popup: first value shows the epitope, second shows the probability score of that epitope.</p>

      {#each allResults as seqResult}
        <div class="seq-result-block">
          <div class="seq-header-row">
            <span class="seq-id">>{seqResult.header}</span>
            <span class="seq-len">Sequence length={seqResult.seqLen}</span>
          </div>

          <button class="download-btn" onclick={() => downloadResults(seqResult)}>Download results as a text file</button>
          {#if pinnedData}
            <div class="selected-action-bar">
              <div class="sel-info">
                <span class="sel-label">Selected Window (Pos {pinnedData.bestPosition}):</span>
                <code class="sel-frag">{pinnedData.bestFragment}</code>
                <span class="sel-prob">({pinnedData.bestProbability.toFixed(1)}% prob)</span>
              </div>
              <button 
                type="button" 
                class="action-bar-btn" 
                onclick={() => onAnalyzeMutant(pinnedData.bestFragment)}
              >
                ↗ Analyze in Peptide Mutants
              </button>
            </div>
          {/if}
          <div class="seq-colorview">
            {#each seqResult.residueData as rd}
              <span
                class="aa"
                class:pinned={pinnedData && pinnedData.position === rd.position}
                style="color: {probabilityToColor(rd.probability)}"
                onmouseenter={(e) => onResidueEnter(e, rd)}
                onmouseleave={onResidueLeave}
                onclick={(e) => onResidueClick(e, rd)}
                title="Click to lock popup for {rd.bestFragment}"
              >{rd.aa}</span>
            {/each}
          </div>

          <details class="detail-table">
            <summary>Show detailed window-by-window results ({seqResult.windowResults.length} windows, {predictMs}ms)</summary>
            <div class="table-wrap">
              <table>
                <thead>
                  <tr><th>Pos</th><th>{currentWindow}-mer</th><th>Score</th><th>Prob %</th><th>Prediction</th></tr>
                </thead>
                <tbody>
                  {#each seqResult.windowResults as r}
                    <tr>
                      <td>{r.position}</td>
                      <td class="mono">{r.fragment}</td>
                      <td class="mono">{r.score.toFixed(4)}</td>
                      <td class="mono">{r.probability.toFixed(1)}</td>
                      <td><span class="badge" class:epitope={r.isEpitope}>{r.isEpitope ? 'Epitope' : 'Non-epitope'}</span></td>
                    </tr>
                  {/each}
                </tbody>
              </table>
            </div>
          </details>
        </div>
      {/each}
    </div>
  {/if}
</section>

{#if tooltipData && !pinnedData}
  <div class="tooltip hover-tooltip" style="left:{tooltipX}px; top:{tooltipY}px">
    <span class="tt-frag">{tooltipData.bestFragment}</span>
    <span class="tt-sep">|</span>
    <span class="tt-prob">{tooltipData.bestProbability.toFixed(1)}%</span>
    <span class="tt-sep">|</span>
    <span class="tt-pos">Pos {tooltipData.bestPosition}</span>
    <div class="tt-hint">Click residue to lock &amp; analyze</div>
  </div>
{/if}

{#if pinnedData}
  <div 
    class="tooltip pinned-popover" 
    style="left:{tooltipX}px; top:{tooltipY}px"
    onclick={(e) => e.stopPropagation()}
  >
    <div class="popover-header">
      <span>Locked Window (Pos {pinnedData.bestPosition})</span>
      <button type="button" class="close-pop-btn" onclick={closePinned}>&times;</button>
    </div>
    <div class="popover-body">
      <button 
        type="button" 
        class="tt-frag-btn" 
        onclick={(e) => {
          e.stopPropagation();
          const frag = pinnedData.bestFragment;
          closePinned();
          onAnalyzeMutant(frag);
        }}
      >
        <span class="frag-text">{pinnedData.bestFragment}</span>
        <span class="tt-action">↗ Analyze Mutants</span>
      </button>
    </div>
    <div class="popover-footer">
      <span>Probability: <strong>{pinnedData.bestProbability.toFixed(1)}%</strong></span>
    </div>
  </div>
{/if}

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
  .length-picker { margin-top: 0.75rem; display: flex; align-items: center; gap: 0.5rem; font-size: 0.9rem; }
  .inline-panel { display: flex; align-items: center; gap: 0.75rem; flex-wrap: wrap; }
  .status { font-family: var(--font-mono); font-size: 0.8rem; color: var(--color-quiet); margin: 0.75rem 0 0; }
  .status.ready { color: var(--color-structural); }
  .error { font-family: var(--font-mono); font-size: 0.8rem; color: #b3261e; margin: 0.75rem 0 0; }
  .actions { margin-bottom: 2rem; }
  .predict-btn { font-family: var(--font-body); font-weight: 500; font-size: 0.95rem; background: var(--color-structural); color: var(--color-paper-raised); border: none; border-radius: 8px; padding: 0.8rem 1.4rem; cursor: pointer; }
  .predict-btn:disabled { background: var(--color-hairline); color: var(--color-quiet); cursor: not-allowed; }
  .results-section { margin-top: 1rem; }
  .results-section h2 { color: var(--color-structural); font-size: 1.4rem; margin-bottom: 1rem; }
  .score-debug { font-family: var(--font-mono); font-size: 0.78rem; color: #92660f; background: #fef9ec; border: 1px solid #e8d5a0; border-radius: 6px; padding: 0.6rem 0.85rem; margin: 0 0 1rem; line-height: 1.5; }
  .legend-label { font-size: 0.85rem; color: var(--color-ink); margin: 0 0 0.5rem; }
  .legend-bar { display: flex; gap: 2px; margin-bottom: 0.5rem; }
  .legend-swatch { display: flex; align-items: center; justify-content: center; padding: 0.35rem 0.8rem; border-radius: 3px; }
  .legend-swatch small { color: #fff; font-family: var(--font-mono); font-size: 0.7rem; font-weight: 500; }
  .legend-hint { font-size: 0.8rem; color: var(--color-quiet); margin: 0 0 1.5rem; line-height: 1.4; }
  .seq-result-block { margin-bottom: 2.5rem; }
  .seq-header-row { font-family: var(--font-mono); font-size: 0.95rem; padding: 0.5rem 0.75rem; background: #e8f0ec; border: 1px solid var(--color-hairline); border-bottom: none; border-radius: 6px 6px 0 0; display: flex; justify-content: space-between; flex-wrap: wrap; gap: 0.5rem; }
  .seq-id { color: #cc0000; font-weight: 600; }
  .seq-len { color: var(--color-structural); }
  .download-btn { font-family: var(--font-body); font-size: 0.8rem; margin: 0.75rem 0; padding: 0.4rem 0.8rem; background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 6px; cursor: pointer; }
  .download-btn:hover { border-color: var(--color-structural); }
  .seq-colorview { font-family: var(--font-mono); font-size: 0.9rem; line-height: 1.65; letter-spacing: 0.01em; word-break: break-all; background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 0 0 6px 6px; padding: 0.75rem; user-select: text; }
  .aa { cursor: default; font-weight: 600; transition: background 0.1s; border-radius: 1px; padding: 0 0.5px; }
  .aa:hover { background: rgba(0, 0, 0, 0.08); }
  .tooltip { position: fixed; transform: translate(-50%, -100%); background: #cc0000; color: #fff; font-family: var(--font-mono); font-size: 0.78rem; font-weight: 500; padding: 0.35rem 0.65rem; border-radius: 4px; white-space: nowrap; pointer-events: none; z-index: 1000; box-shadow: 0 2px 8px rgba(0,0,0,0.25); }
  .tooltip::after { content: ''; position: absolute; top: 100%; left: 50%; transform: translateX(-50%); border: 5px solid transparent; border-top-color: #cc0000; }
  .tt-sep { opacity: 0.5; margin: 0 0.3rem; }
  .detail-table { margin-top: 1rem; }
  .detail-table summary { font-family: var(--font-body); font-size: 0.85rem; font-weight: 500; color: var(--color-structural); cursor: pointer; padding: 0.5rem 0; }
  .table-wrap { overflow-x: auto; margin-top: 0.5rem; }
  table { width: 100%; border-collapse: collapse; font-size: 0.82rem; }
  th { text-align: left; font-weight: 500; color: var(--color-quiet); border-bottom: 1px solid var(--color-hairline); padding: 0.45rem 0.5rem; }
  td { padding: 0.45rem 0.5rem; border-bottom: 1px solid var(--color-hairline); }
  .mono { font-family: var(--font-mono); letter-spacing: 0.02em; }
  .badge { display: inline-block; font-size: 0.72rem; font-weight: 500; padding: 0.2rem 0.5rem; border-radius: 999px; background: #eef1f0; color: var(--color-quiet); }
  .badge.epitope { background: var(--color-epitope-bg); color: var(--color-epitope); }
  @media (max-width: 640px) { .predictor { padding: 2rem 1.25rem 3rem; } .legend-bar { flex-wrap: wrap; } }
.tooltip { 
    position: fixed; 
    transform: translate(-50%, -100%); 
    background: #222222; 
    color: #fff; 
    font-family: var(--font-mono); 
    font-size: 0.8rem; 
    font-weight: 500; 
    padding: 0.4rem 0.7rem; 
    border-radius: 6px; 
    white-space: nowrap; 
    /* Removed pointer-events: none; so it can be clicked! */
    z-index: 1000; 
    box-shadow: 0 4px 12px rgba(0,0,0,0.3); 
    display: flex;
    align-items: center;
    border: 1px solid #444;
  }
  .tooltip::after { 
    content: ''; 
    position: absolute; 
    top: 100%; 
    left: 50%; 
    transform: translateX(-50%); 
    border: 5px solid transparent; 
    border-top-color: #222222; 
  }
  .tt-frag-btn {
    background: rgba(255, 255, 255, 0.1);
    border: 1px solid rgba(255, 255, 255, 0.2);
    color: #ffffff;
    font-family: var(--font-mono);
    font-size: 0.82rem;
    font-weight: 600;
    padding: 0.2rem 0.5rem;
    border-radius: 4px;
    cursor: pointer;
    transition: all 0.15s ease;
    display: inline-flex;
    align-items: center;
    gap: 0.35rem;
  }
  .tt-frag-btn:hover {
    background: var(--color-structural, #006600);
    border-color: #008800;
    transform: scale(1.03);
  }
  .tt-action {
    font-family: var(--font-body);
    font-size: 0.68rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    opacity: 0.85;
    background: rgba(0, 0, 0, 0.3);
    padding: 0.1rem 0.35rem;
    border-radius: 3px;
  }
  .tt-sep { opacity: 0.4; margin: 0 0.4rem; }
.aa { 
    cursor: pointer; 
    font-weight: 600; 
    transition: all 0.1s; 
    border-radius: 2px; 
    padding: 0 1px;
    display: inline-block;
  }
  .aa:hover { 
    background: rgba(0, 0, 0, 0.12); 
    transform: scale(1.15); 
  }
  .aa.pinned {
    background: var(--color-structural, #006600);
    color: #ffffff !important;
    border-radius: 3px;
    box-shadow: 0 0 0 2px rgba(0, 102, 0, 0.3);
  }

  /* Inline Action Bar Styling */
  .selected-action-bar {
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 0.75rem;
    background: #f0f7f2;
    border: 1px solid var(--color-structural, #006600);
    border-radius: 6px 6px 0 0;
    padding: 0.6rem 1rem;
    margin-bottom: -1px; /* Connects flush with sequence box below */
  }
  .sel-info { font-size: 0.85rem; color: #222; }
  .sel-label { font-weight: 600; color: var(--color-structural, #006600); }
  .sel-frag { font-family: var(--font-mono); font-weight: 700; font-size: 0.95rem; background: #fff; padding: 0.15rem 0.4rem; border-radius: 4px; border: 1px solid #ccc; margin: 0 0.3rem; }
  .sel-prob { color: #555; font-weight: 500; }
  .action-bar-btn {
    background: var(--color-structural, #006600);
    color: #fff;
    border: none;
    font-family: var(--font-body);
    font-size: 0.85rem;
    font-weight: 600;
    padding: 0.45rem 0.9rem;
    border-radius: 6px;
    cursor: pointer;
    transition: all 0.15s ease;
  }
  .action-bar-btn:hover { background: #008800; transform: translateY(-1px); box-shadow: 0 2px 6px rgba(0,0,0,0.15); }

  /* Base Tooltip Styling */
  .tooltip { 
    position: fixed; 
    transform: translate(-50%, -100%); 
    background: #222222; 
    color: #fff; 
    font-family: var(--font-mono); 
    font-size: 0.8rem; 
    border-radius: 6px; 
    white-space: nowrap; 
    z-index: 1000; 
    box-shadow: 0 4px 14px rgba(0,0,0,0.35); 
    border: 1px solid #444;
  }
  .tooltip::after { 
    content: ''; 
    position: absolute; 
    top: 100%; 
    left: 50%; 
    transform: translateX(-50%); 
    border: 5px solid transparent; 
    border-top-color: #222222; 
  }

  /* Hover Preview Mode (Never catches the mouse) */
  .hover-tooltip {
    padding: 0.4rem 0.7rem;
    pointer-events: none !important;
    text-align: center;
  }
  .tt-hint { font-family: var(--font-body); font-size: 0.65rem; color: #aaa; margin-top: 0.2rem; border-top: 1px solid #333; padding-top: 0.15rem; }
  .tt-sep { opacity: 0.4; margin: 0 0.4rem; }

  /* Pinned Popover Mode (Solid surface, safe to click) */
  .pinned-popover {
    background: #1a1a1a;
    border: 1px solid var(--color-structural, #006600);
    box-shadow: 0 8px 24px rgba(0,0,0,0.6);
    padding: 0;
    min-width: 220px;
    cursor: default;
    pointer-events: auto !important;
  }
  .pinned-popover::after { border-top-color: var(--color-structural, #006600); }
  .popover-header { display: flex; justify-content: space-between; align-items: center; background: #2a2a2a; padding: 0.4rem 0.7rem; font-size: 0.72rem; color: #ccc; border-top-left-radius: 5px; border-top-right-radius: 5px; border-bottom: 1px solid #333; }
  .close-pop-btn { background: none; border: none; color: #888; font-size: 1.1rem; line-height: 1; cursor: pointer; padding: 0; }
  .close-pop-btn:hover { color: #fff; }
  .popover-body { padding: 0.7rem; text-align: center; }
  .popover-footer { background: #222; padding: 0.35rem 0.7rem; font-size: 0.72rem; color: #aaa; border-bottom-left-radius: 5px; border-bottom-right-radius: 5px; border-top: 1px solid #333; text-align: center; }
  .popover-footer strong { color: #fff; }

  /* Clickable Action Button inside popover */
  .tt-frag-btn {
    background: var(--color-structural, #006600);
    border: 1px solid #008800;
    color: #ffffff;
    font-family: var(--font-mono);
    font-size: 0.9rem;
    font-weight: 600;
    padding: 0.5rem 0.8rem;
    border-radius: 6px;
    cursor: pointer;
    transition: all 0.15s ease;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.3rem;
    width: 100%;
    box-sizing: border-box;
  }
  .tt-frag-btn:hover { background: #008800; transform: translateY(-1px); box-shadow: 0 4px 10px rgba(0, 136, 0, 0.3); }
  .frag-text { letter-spacing: 0.05em; }
  .tt-action { font-family: var(--font-body); font-size: 0.7rem; text-transform: uppercase; letter-spacing: 0.04em; background: rgba(0, 0, 0, 0.25); padding: 0.15rem 0.5rem; border-radius: 999px; width: 100%; }
</style>