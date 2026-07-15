<script>
  function imgUrl(file) {
    return `${import.meta.env.BASE_URL}images/${file}`;
  }

  let selectedImage = $state(null);

  function openModal(imgSrc, title) {
    selectedImage = { src: imgSrc, title };
  }

  function closeModal() {
    selectedImage = null;
  }
</script>

<section class="algorithm-view">
  <h2>Algorithms &amp; Datasets in LBtope</h2>
  <p class="intro">
    To predict linear B-cell epitopes with high accuracy, LBtope implements models trained on five distinct datasets using Support Vector Machines (SVM) and Dipeptide Composition (DPC) features.
  </p>

  <div class="dataset-grid">
    <div class="dataset-card" id="fix">
      <h3>1. LBtope_Fixed</h3>
      <p>
        Developed on the <strong>LBtope_Fixed</strong> dataset, this model is optimized for predicting fixed-length linear B-cell epitopes of exactly 20 residues. It contains <strong>12,063 unique positive patterns</strong> and <strong>20,589 unique negative patterns</strong>. Extension and termination techniques were applied to standardize experimental epitopes to 20-residue lengths.
      </p>
      <button class="view-diagram-btn" onclick={() => openModal(imgUrl('algotithm1.jpg'), 'Algorithm: LBtope_Fixed')}>
        View Pipeline Schematic
      </button>
    </div>

    <div class="dataset-card" id="var">
      <h3>2. LBtope_Variable</h3>
      <p>
        Designed for predicting variable-length B-cell epitopes, this model was trained on the <strong>LBtope_Variable</strong> dataset containing <strong>14,876 unique B-cell epitopes</strong> and <strong>23,321 unique non-B-cell epitopes</strong> of varying native lengths.
      </p>
      <button class="view-diagram-btn" onclick={() => openModal(imgUrl('algotithm2.jpg'), 'Algorithm: LBtope_Variable')}>
        View Pipeline Schematic
      </button>
    </div>

    <div class="dataset-card" id="red">
      <h3>3. LBtope_Confirm</h3>
      <p>
        Built on the highly rigorous <strong>LBtope_Confirm</strong> dataset, which includes only those epitopes and non-epitopes that have been experimentally validated by <strong>two or more independent studies</strong>. This gold-standard dataset contains <strong>1,042 unique B-cell epitopes</strong> and <strong>1,795 unique non-B-cell epitopes</strong>.
      </p>
      <button class="view-diagram-btn" onclick={() => openModal(imgUrl('algotithm3.jpg'), 'Algorithm: LBtope_Confirm')}>
        View Pipeline Schematic
      </button>
    </div>

    <div class="dataset-card" id="fixnr">
      <h3>4. LBtope_Fixed_non_redundant</h3>
      <p>
        Developed on the <strong>LBtope_Fixed_non_redundant</strong> dataset, which was reduced to an 80% non-redundant threshold using CD-HIT to remove sequence bias. It contains <strong>7,824 unique positive patterns</strong> and <strong>7,853 unique negative patterns</strong> (each 20 residues long).
      </p>
    </div>

    <div class="dataset-card" id="fixnrv">
      <h3>5. LBtope_Variable_non_redundant</h3>
      <p>
        Trained on the <strong>LBtope_Variable_non_redundant</strong> dataset, also refined to an 80% non-redundant threshold via CD-HIT. This variable-length dataset contains <strong>8,011 unique positive patterns</strong> and <strong>10,868 unique negative patterns</strong>.
      </p>
    </div>
  </div>
</section>

{#if selectedImage}
  <div class="modal-backdrop" onclick={closeModal}>
    <div class="modal-content" onclick={(e) => e.stopPropagation()}>
      <div class="modal-header">
        <h4>{selectedImage.title}</h4>
        <button class="close-btn" onclick={closeModal}>&times;</button>
      </div>
      <div class="modal-body">
        <img src={selectedImage.src} alt={selectedImage.title} />
      </div>
    </div>
  </div>
{/if}

<style>
  .algorithm-view { max-width: 880px; margin: 0 auto; padding: 3rem 2rem 5rem; width: 100%; }
  h2 { font-family: var(--font-display); font-weight: 600; font-size: 1.75rem; color: var(--color-structural); margin: 0 0 0.75rem; }
  .intro { color: var(--color-quiet); font-size: 1rem; line-height: 1.6; margin: 0 0 2rem; }
  .dataset-grid { display: flex; flex-direction: column; gap: 1.5rem; }
  .dataset-card { background: var(--color-paper-raised); border: 1px solid var(--color-hairline); border-left: 4px solid var(--color-structural); border-radius: 8px; padding: 1.5rem; box-shadow: 0 1px 4px rgba(0, 0, 0, 0.03); }
  .dataset-card h3 { font-family: var(--font-display); font-size: 1.15rem; color: var(--color-ink); margin: 0 0 0.75rem; }
  .dataset-card p { font-size: 0.92rem; line-height: 1.65; color: var(--color-quiet); margin: 0 0 1rem; }
  .dataset-card p:last-child { margin-bottom: 0; }
  .dataset-card strong { color: var(--color-ink); }
  
  .view-diagram-btn { font-family: var(--font-body); font-size: 0.85rem; font-weight: 500; color: var(--color-structural); background: rgba(0, 102, 0, 0.06); border: 1px solid var(--color-structural); border-radius: 6px; padding: 0.45rem 0.9rem; cursor: pointer; transition: all 0.15s ease; }
  .view-diagram-btn:hover { background: var(--color-structural); color: #fff; }

  /* Modal Zoom Styles */
  .modal-backdrop { position: fixed; top: 0; left: 0; width: 100vw; height: 100vh; background: rgba(0, 0, 0, 0.75); display: flex; align-items: center; justify-content: center; z-index: 2000; padding: 1.5rem; }
  .modal-content { background: var(--color-paper-raised); border-radius: 10px; max-width: 900px; width: 100%; max-height: 90vh; display: flex; flex-direction: column; overflow: hidden; box-shadow: 0 10px 30px rgba(0,0,0,0.5); }
  .modal-header { display: flex; align-items: center; justify-content: space-between; padding: 1rem 1.5rem; border-bottom: 1px solid var(--color-hairline); background: #f8f9fa; }
  .modal-header h4 { margin: 0; font-family: var(--font-display); font-size: 1.1rem; color: var(--color-structural); }
  .close-btn { background: none; border: none; font-size: 1.75rem; line-height: 1; color: var(--color-quiet); cursor: pointer; }
  .close-btn:hover { color: #cc0000; }
  .modal-body { padding: 1rem; overflow-y: auto; text-align: center; background: #fff; }
  .modal-body img { max-width: 100%; height: auto; border-radius: 4px; }
  
  @media (max-width: 640px) { .algorithm-view { padding: 2rem 1.25rem 3rem; } }
</style>