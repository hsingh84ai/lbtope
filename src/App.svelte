<script>
  import SequenceStrip from './lib/SequenceStrip.svelte';
  import AntigenPredictor from './lib/AntigenPredictor.svelte';
  import PeptidePredictor from './lib/PeptidePredictor.svelte';
  import MutantPredictor from './lib/MutantPredictor.svelte';  
  import AlgorithmView from './lib/AlgorithmView.svelte';
  import HelpView from './lib/HelpView.svelte';
  let view = $state('home');
  let mutantTargetPeptide = $state('');

  function imgUrl(file) {
    return `${import.meta.env.BASE_URL}images/${file}`;
  }

  // Triggered when a user clicks "Analyze Mutants" in Antigen Sequence
  function handleAnalyzeMutant(peptideSeq) {
    mutantTargetPeptide = peptideSeq;
    view = 'mutant'; // Automatically switch to the Peptide Mutants tab!
  }
</script>

<div class="page">
  <header class="site-header">
    <button class="wordmark" onclick={() => (view = 'home')} aria-label="LBTope home">
      <img src={imgUrl('lbtope1.gif')} alt="LBtope Logo" class="legacy-logo" onerror={(e) => e.target.style.display='none'} />
      <span>LB<span class="accent">Tope</span></span>
    </button>
    <nav>
      <button class="nav-link" class:active={view === 'antigen'} onclick={() => (view = 'antigen')}>Antigen Sequence</button>
      <button class="nav-link" class:active={view === 'peptide'} onclick={() => (view = 'peptide')}>Multiple Peptides</button>
      <button class="nav-link" class:active={view === 'mutant'} onclick={() => (view = 'mutant')}>Peptide Mutants</button>
      <button class="nav-link" class:active={view === 'algorithm'} onclick={() => (view = 'algorithm')}>Algorithm</button>
      <button class="nav-link" class:active={view === 'help'} onclick={() => (view = 'help')}>Help</button>
    </nav>
  </header>
  <main>
    {#if view === 'home'}
      <section class="hero">
        <h1>LBtope: Prediction of Linear B-cell Epitopes</h1>
        <p class="hero-subtitle">
          An improved method for predicting and designing B-cell epitopes using antigen primary sequence data, running entirely client-side in your browser.
        </p>

        <div class="bio-context-card">
          <p>
            Predicting B-cell epitopes (antigenic regions) with high accuracy is one of the major challenges in designing subunit/peptide vaccines and immunotherapies. Earlier computational methods (such as <a href="https://webs.iiitd.edu.in/raghava/abcpred/" target="_blank" rel="noopener">ABCpred</a>) suffered from two major limitations: they were developed on small datasets, and they relied on random peptides as assumed non-epitopes.
          </p>
          <p>
            <strong>LBtope</strong> overcomes these limitations by leveraging machine learning models (SVM, IBk) trained on a massive dataset from the <a href="https://www.iedb.org/" target="_blank" rel="noopener">IEDB database</a>—featuring <strong>12,063 experimentally validated epitopes</strong> and <strong>20,589 experimentally validated non-epitopes</strong>. For the first time, confirmed non-B-cell epitopes were utilized as negative controls to build highly selective prediction models.
          </p>
        </div>
      </section>

      <section class="modules">
        <h2 class="section-label">Prediction &amp; Design Modules</h2>
        <div class="card-grid">
          <button class="module-card" onclick={() => (view = 'antigen')}>
            <h3>Antigen Sequence</h3>
            <p>Submit one or more protein/antigen sequences in FASTA format. Predict B-cell epitopes using fixed-length (20-mer) or variable-length models.</p>
          </button>
          <button class="module-card" onclick={() => (view = 'peptide')}>
            <h3>Multiple Peptides</h3>
            <p>Submit a list of discrete peptides and classify each as an epitope or non-epitope in a single high-speed batch.</p>
          </button>
          <button class="module-card" onclick={() => (view = 'mutant')}>
            <h3>Peptide Mutants</h3>
            <p>Perform a systematic point-mutation scan across a peptide to identify single-residue substitutions that optimize epitope affinity.</p>
          </button>
        </div>
      </section>

      <section class="signature"><SequenceStrip /></section>

      <footer class="cite-footer">
        <div class="reference-box">
          <strong>Reference:</strong> Singh H, Ansari HR, Raghava GPS (2013). Improved Method for Linear B-Cell Epitope Prediction Using Antigen's Primary Sequence. 
          <a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0062216" target="_blank" rel="noopener">
            <em>PLoS ONE</em> 8(5): e62216
          </a>
        </div>
      </footer>
    {:else if view === 'antigen'}
      <AntigenPredictor onAnalyzeMutant={handleAnalyzeMutant} />
    {:else if view === 'peptide'}
      <PeptidePredictor />
    {:else if view === 'mutant'}
      <MutantPredictor initialPeptide={mutantTargetPeptide} />
    {:else if view === 'algorithm'}
      <AlgorithmView />
    {:else if view === 'help'}
      <HelpView />
    {/if}
  </main>
</div>

<style>
  .page { min-height: 100vh; display: flex; flex-direction: column; }
  .site-header { display: flex; align-items: center; justify-content: space-between; padding: 1.25rem 2rem; border-bottom: 1px solid var(--color-hairline); flex-wrap: wrap; gap: 0.75rem; }
  .wordmark { background: none; border: none; cursor: pointer; padding: 0; font-family: var(--font-display); font-weight: 700; font-size: 1.15rem; letter-spacing: 0.02em; color: var(--color-ink); }
  .wordmark .accent { color: var(--color-structural); }
  nav { display: flex; gap: 0.25rem; flex-wrap: wrap; }
  .nav-link { font-family: var(--font-body); font-weight: 500; font-size: 0.85rem; background: none; color: var(--color-quiet); border: none; border-radius: 999px; padding: 0.45rem 0.9rem; cursor: pointer; transition: color 0.15s ease, background 0.15s ease; }
  .nav-link:hover { color: var(--color-ink); background: var(--color-hairline); }
  .nav-link.active { color: var(--color-paper-raised); background: var(--color-structural); }
  main { flex: 1; display: flex; flex-direction: column; }
  .hero { padding: 5rem 2rem 2rem; max-width: 780px; margin: 0 auto; width: 100%; }
  .hero h1 { font-family: var(--font-display); font-weight: 600; font-size: clamp(2rem, 4.5vw, 3rem); line-height: 1.15; margin: 0 0 1.25rem; color: var(--color-ink); }
  .hero-desc { font-size: 1.05rem; line-height: 1.6; color: var(--color-quiet); max-width: 60ch; margin: 0; }
  .modules { max-width: 780px; margin: 0 auto; padding: 2.5rem 2rem 1rem; width: 100%; }
  .section-label { font-family: var(--font-display); font-weight: 600; font-size: 1.1rem; color: var(--color-structural); margin: 0 0 1rem; }
  .card-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; }
  .module-card { text-align: left; background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-radius: 12px; padding: 1.5rem; cursor: pointer; transition: border-color 0.15s ease, box-shadow 0.15s ease; }
  .module-card:hover { border-color: var(--color-structural); box-shadow: 0 2px 8px rgba(0,0,0,0.06); }
  .module-card h3 { font-family: var(--font-display); font-weight: 600; font-size: 1rem; margin: 0 0 0.5rem; }
  .module-card p { font-size: 0.85rem; line-height: 1.5; color: var(--color-quiet); margin: 0; }
  .signature { max-width: 780px; margin: 0 auto; padding: 2rem 2rem 1rem; width: 100%; }
  .cite-footer { max-width: 780px; margin: 0 auto; padding: 1rem 2rem 4rem; width: 100%; }
  .cite-footer a { font-family: var(--font-mono); font-size: 0.8rem; color: var(--color-quiet); text-decoration: none; border-bottom: 1px solid var(--color-hairline); }
  .cite-footer a:hover { color: var(--color-structural); border-color: var(--color-structural); }
  .placeholder { max-width: 780px; margin: 0 auto; padding: 4rem 2rem; width: 100%; }
  .placeholder h2 { font-family: var(--font-display); font-weight: 600; font-size: 1.75rem; margin: 0 0 0.75rem; }
  .placeholder p { color: var(--color-quiet); margin: 0 0 1.5rem; }
  .back-link { font-family: var(--font-body); font-size: 0.9rem; background: none; border: none; color: var(--color-structural); cursor: pointer; padding: 0; }
  .back-link:hover { text-decoration: underline; }
  @media (max-width: 640px) { .site-header { padding: 1rem 1.25rem; } .hero { padding: 3rem 1.25rem 1.5rem; } .modules, .signature, .cite-footer, .placeholder { padding-left: 1.25rem; padding-right: 1.25rem; } }
.hero-subtitle {
    font-size: 1.1rem;
    line-height: 1.5;
    color: var(--color-structural);
    font-weight: 500;
    margin: 0 0 1.5rem;
  }
  .bio-context-card {
    background: var(--color-paper-raised);
    border: 1px solid var(--color-hairline);
    border-left: 4px solid var(--color-structural);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    font-size: 0.92rem;
    line-height: 1.65;
    color: var(--color-quiet);
    text-align: justify;
    box-shadow: 0 1px 4px rgba(0, 0, 0, 0.03);
  }
  .bio-context-card p {
    margin: 0 0 0.85rem;
  }
  .bio-context-card p:last-child {
    margin: 0;
  }
  .bio-context-card strong {
    color: var(--color-ink);
  }
  .bio-context-card a {
    color: var(--color-structural);
    text-decoration: underline;
  }
  .reference-box {
    background: rgba(0, 0, 0, 0.02);
    border: 1px solid var(--color-hairline);
    border-radius: 6px;
    padding: 0.85rem 1.25rem;
    font-size: 0.82rem;
    line-height: 1.5;
    color: var(--color-quiet);
    text-align: center;
  }
  .reference-box strong {
    color: var(--color-ink);
  }
:global(:root) {
    /* Match legacy LBtope PHP theme colors */
    --color-structural: #006600;
    --color-structural-dim: #004400;
    --color-epitope-bg: #e8f5e9;
    --color-epitope: #2e7d32;
  }
.wordmark {
    display: flex;
    align-items: center;
    gap: 0.75rem;
  }
  .legacy-logo {
    height: 42px;
    width: auto;
    border-radius: 6px;
    /* Solid green badge background to make transparent GIF visible */
    background: var(--color-structural);
    padding: 4px 8px;
    border: 1px solid var(--color-structural-dim);
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.15);
    object-fit: contain;
  }
</style>
