<script>
  function imgUrl(file) {
    return `${import.meta.env.BASE_URL}images/${file}`;
  }

  let activeSection = $state('antigen'); // 'antigen' | 'peptides' | 'mutants'
  let selectedImage = $state(null);

  function openModal(imgSrc, title) {
    selectedImage = { src: imgSrc, title };
  }

  function closeModal() {
    selectedImage = null;
  }
</script>

<section class="help-view">
  <h2>LBtope Help &amp; Documentation</h2>
  <p class="intro">
    LBtope discriminates linear B-cell epitopes from non-epitopes using Support Vector Machine (SVM) models trained on dipeptide composition (DPC) features. Because this modernized version runs 100% client-side in your browser, predictions are instantaneous without server queues or file uploads.
  </p>

  <div class="nav-pills">
    <button type="button" class:active={activeSection === 'antigen'} onclick={() => activeSection = 'antigen'}>
      1. Antigen Sequence
    </button>
    <button type="button" class:active={activeSection === 'peptides'} onclick={() => activeSection = 'peptides'}>
      2. Multiple Peptides
    </button>
    <button type="button" class:active={activeSection === 'mutants'} onclick={() => activeSection = 'mutants'}>
      3. Peptide Mutants
    </button>
  </div>

  <div class="content-card">
    {#if activeSection === 'antigen'}
      <h3>Antigen Sequence Prediction</h3>
      <p>
        The Antigen Sequence module scans protein or antigen sequences using a sliding window approach to identify continuous B-cell epitope motifs. The overall accuracy of the underlying SVM models is approximately 81%.
      </p>

      <h4>Sequence Input &amp; Format</h4>
      <p>
        Amino acid sequences must be entered in standard <strong>FASTA format</strong> (starting with a <code>&gt;header</code> line followed by the sequence). Non-standard characters (such as numbers or punctuation like <code>*&amp;^%$@#!()_+~=</code>) are automatically ignored during processing.
      </p>

      <h4>Available Prediction Models</h4>
      <ul class="model-list">
        <li><strong>LBtope_Fixed:</strong> A 20-mer model developed on the LBtope_Fixed dataset. A fixed window length of 20 residues is slid across the query sequence.</li>
        <li><strong>LBtope_Fixed_non_redundant:</strong> Trained on an 80% non-redundant dataset generated via CD-HIT. Uses a fixed 20-residue scanning window.</li>
        <li><strong>LBtope_Variable:</strong> Developed on variable-length epitope data. It allows users to select a custom sliding window size between 5 and 30 residues (defaulting to 15 residues as established by Kringelum et al., 2013).</li>
        <li><strong>LBtope_Variable_non_redundant:</strong> Trained on the 80% non-redundant variable-length dataset. Supports custom window sizing.</li>
        <li><strong>LBtope_Confirm:</strong> Developed exclusively on gold-standard epitopes verified by at least two independent experimental studies. Supports custom window sizing.</li>
      </ul>

      <h4>Probability Threshold &amp; Visual Output</h4>
      <p>
        During prediction, each overlapping window receives a raw SVM decision score, which is linearly mapped to a 0–100% probability scale. Users can adjust the <strong>Probability Threshold</strong> dropdown (defaulting to 60%) to control visual sensitivity. Setting a higher threshold yields higher specificity and confidence.
      </p>
      <p>
        The results view displays a color-coded sequence block where each amino acid is colored according to its centered window probability. Hovering over any residue reveals an interactive tooltip displaying the exact epitope fragment, raw probability percentage, and position.
      </p>

      <div class="screenshot-grid">
        <button type="button" class="img-preview-btn" onclick={() => openModal(imgUrl('result_protein.jpg'), 'Legacy Output: Color-Coded Protein Sequence')}>
          <img src={imgUrl('result_protein.jpg')} alt="Protein Sequence Output" />
          <span>View Legacy Sequence Output</span>
        </button>
        <button type="button" class="img-preview-btn" onclick={() => openModal(imgUrl('result_protein_popup.jpg'), 'Legacy Output: Interactive Hover Tooltip')}>
          <img src={imgUrl('result_protein_popup.jpg')} alt="Hover Tooltip Output" />
          <span>View Legacy Tooltip Popup</span>
        </button>
      </div>

    {:else if activeSection === 'peptides'}
      <h3>Multiple Peptides Batch Classification</h3>
      <p>
        The Multiple Peptides module evaluates discrete peptides of varying lengths simultaneously, classifying each sequence as an epitope or non-epitope without sliding windows.
      </p>

      <h4>Sequence Input &amp; Format</h4>
      <p>
        Submit multiple peptides in standard <strong>FASTA format</strong>. Each sequence should represent a single discrete peptide (e.g., ranging from 8 to 30 amino acids). Non-standard symbols and spacing are cleaned automatically.
      </p>

      <h4>Available Prediction Models</h4>
      <p>
        Because submitted peptides vary in native length, this module exclusively utilizes variable-length models:
      </p>
      <ul class="model-list">
        <li><strong>LBtope_Variable:</strong> Trained on the complete variable-length dataset.</li>
        <li><strong>LBtope_Variable_non_redundant:</strong> Trained on the CD-HIT 80% non-redundant dataset.</li>
        <li><strong>LBtope_Confirm:</strong> Trained on epitopes confirmed by two or more independent studies.</li>
      </ul>

      <h4>Tabular Batch Output</h4>
      <p>
        Results are generated instantly and displayed in a structured data table showing the Peptide ID, Sequence, Length, Raw SVM Score, Probability Percentage, and Classification Badge. Rows meeting or exceeding your selected probability threshold are highlighted automatically. The entire result set can be exported to a TSV file for downstream analysis.
      </p>

      <div class="screenshot-grid">
        <button type="button" class="img-preview-btn" onclick={() => openModal(imgUrl('result_multiple.JPG'), 'Legacy Output: Multiple Peptides Table')}>
          <img src={imgUrl('result_multiple.JPG')} alt="Multiple Peptides Output" onerror={(e) => { e.target.src = imgUrl('result_multiple.jpg'); }} />
          <span>View Legacy Tabular Output</span>
        </button>
      </div>

    {:else if activeSection === 'mutants'}
      <h3>Peptide Mutants Optimization Scan</h3>
      <p>
        The Peptide Mutants module performs systematic single-point amino acid substitutions across a query peptide. It is designed to help researchers identify the minimal mutations required to engineer a weak peptide into a high-affinity B-cell epitope.
      </p>

      <h4>Sequence Input &amp; Execution</h4>
      <p>
        Enter a single peptide sequence in plain text format (e.g., <code>STAIHADQLTPAWRIYSTGNN</code>). The application generates all L × 20 possible single-residue variations and computes dipeptide composition features for each mutant simultaneously.
      </p>

      <h4>Available Prediction Models</h4>
      <p>
        Like the batch peptide tool, the mutant scanner utilizes variable-length models (<strong>LBtope_Variable</strong>, <strong>LBtope_Variable_non_redundant</strong>, and <strong>LBtope_Confirm</strong>) to accurately score mutated peptides of any length.
      </p>

      <h4>Heatmap &amp; Ranked Outputs</h4>
      <p>
        To make mutation data actionable, results are presented in two complementary views:
      </p>
      <ul>
        <li><strong>Heatmap Grid:</strong> A 2D matrix comparing sequence positions (rows) against all 20 amino acid substitutions (columns). Cells are color-coded by the change in probability (ΔProb), where green indicates improved binding affinity and red indicates degradation.</li>
        <li><strong>Top Enhancements:</strong> A ranked table prioritizing substitutions that produce the greatest positive gain in epitope probability over the wild-type sequence.</li>
      </ul>

      <div class="screenshot-grid">
        <button type="button" class="img-preview-btn" onclick={() => openModal(imgUrl('result_mutants.jpg'), 'Legacy Output: Mutant Prediction Table')}>
          <img src={imgUrl('result_mutants.jpg')} alt="Mutants Output" />
          <span>View Legacy Mutants Table</span>
        </button>
      </div>
    {/if}
  </div>
</section>

{#if selectedImage}
  <div class="modal-backdrop" onclick={closeModal}>
    <div class="modal-content" onclick={(e) => e.stopPropagation()}>
      <div class="modal-header">
        <h4>{selectedImage.title}</h4>
        <button type="button" class="close-btn" onclick={closeModal}>&times;</button>
      </div>
      <div class="modal-body">
        <img src={selectedImage.src} alt={selectedImage.title} />
      </div>
    </div>
  </div>
{/if}

<style>
  .help-view { max-width: 880px; margin: 0 auto; padding: 3rem 2rem 5rem; width: 100%; }
  h2 { font-family: var(--font-display); font-weight: 600; font-size: 1.75rem; color: var(--color-structural); margin: 0 0 0.75rem; }
  .intro { color: var(--color-quiet); font-size: 1rem; line-height: 1.6; margin: 0 0 2rem; }
  
  .nav-pills { display: flex; gap: 0.5rem; margin-bottom: 1.5rem; flex-wrap: wrap; }
  .nav-pills button { font-family: var(--font-body); font-size: 0.9rem; font-weight: 500; background: var(--color-paper-raised); color: var(--color-quiet); border: 1px solid var(--color-hairline); border-radius: 999px; padding: 0.55rem 1.25rem; cursor: pointer; transition: all 0.15s ease; }
  .nav-pills button:hover { border-color: var(--color-structural); color: var(--color-ink); }
  .nav-pills button.active { background: var(--color-structural); color: #fff; border-color: var(--color-structural); }

  .content-card { background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-top: 4px solid var(--color-structural); border-radius: 8px; padding: 2rem; box-shadow: 0 1px 4px rgba(0, 0, 0, 0.03); }
  .content-card h3 { font-family: var(--font-display); font-size: 1.35rem; color: var(--color-ink); margin: 0 0 1rem; border-bottom: 1px solid var(--color-hairline); padding-bottom: 0.5rem; }
  .content-card h4 { font-family: var(--font-display); font-size: 1.05rem; color: var(--color-structural); margin: 1.5rem 0 0.5rem; }
  .content-card p { font-size: 0.95rem; line-height: 1.65; color: var(--color-quiet); margin: 0 0 1rem; text-align: justify; }
  .content-card ul { margin: 0 0 1rem; padding-left: 1.5rem; color: var(--color-quiet); font-size: 0.95rem; line-height: 1.65; }
  .content-card li { margin-bottom: 0.5rem; }
  .content-card strong { color: var(--color-ink); }
  .content-card code { font-family: var(--font-mono); font-size: 0.85rem; background: rgba(0,0,0,0.05); padding: 0.15rem 0.35rem; border-radius: 4px; }

  .screenshot-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 1rem; margin-top: 1.5rem; }
  .img-preview-btn { background: #fff; border: 1px solid var(--color-hairline); border-radius: 6px; padding: 0.5rem; cursor: pointer; display: flex; flex-direction: column; align-items: center; gap: 0.5rem; transition: all 0.15s ease; }
  .img-preview-btn:hover { border-color: var(--color-structural); transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.08); }
  .img-preview-btn img { width: 100%; height: 140px; object-fit: cover; object-position: top; border-radius: 4px; border: 1px solid rgba(0,0,0,0.05); }
  .img-preview-btn span { font-family: var(--font-body); font-size: 0.8rem; font-weight: 500; color: var(--color-structural); }

  /* Modal Zoom Styles */
  .modal-backdrop { position: fixed; top: 0; left: 0; width: 100vw; height: 100vh; background: rgba(0, 0, 0, 0.75); display: flex; align-items: center; justify-content: center; z-index: 2000; padding: 1.5rem; }
  .modal-content { background: var(--color-paper-raised); border-radius: 10px; max-width: 900px; width: 100%; max-height: 90vh; display: flex; flex-direction: column; overflow: hidden; box-shadow: 0 10px 30px rgba(0,0,0,0.5); }
  .modal-header { display: flex; align-items: center; justify-content: space-between; padding: 1rem 1.5rem; border-bottom: 1px solid var(--color-hairline); background: #f8f9fa; }
  .modal-header h4 { margin: 0; font-family: var(--font-display); font-size: 1.1rem; color: var(--color-structural); }
  .close-btn { background: none; border: none; font-size: 1.75rem; line-height: 1; color: var(--color-quiet); cursor: pointer; }
  .close-btn:hover { color: #cc0000; }
  .modal-body { padding: 1rem; overflow-y: auto; text-align: center; background: #fff; }
  .modal-body img { max-width: 100%; height: auto; border-radius: 4px; }

  @media (max-width: 640px) { .help-view { padding: 2rem 1.25rem 3rem; } .content-card { padding: 1.25rem; } }
</style>