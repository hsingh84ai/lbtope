# LBtope: Client-Side Linear B-Cell Epitope Prediction & Design

[![Built with Svelte 5](https://img.shields.io/badge/Svelte_5-%23f1413d.svg?style=flat-square&logo=svelte&logoColor=white)](https://svelte.dev/)
[![Built with Vite](https://img.shields.io/badge/Vite-%23646CFF.svg?style=flat-square&logo=vite&logoColor=white)](https://vitejs.dev/)
[![Client-Side ML](https://img.shields.io/badge/100%25_Client--Side-Zero_Server_Uploads-006600.svg?style=flat-square)](#)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg?style=flat-square)](http://www.gnu.org/licenses/gpl-3.0)

A modernized, high-performance, 100% client-side web application for predicting and designing linear B-cell epitopes. 

This repository is a browser-native port of the legacy **LBtope** server originally developed by Dr. G.P.S. Raghava's group at IMTECH. By converting the legacy SVM-light models into a custom optimized binary format and implementing feature extraction directly in modern JavaScript/Svelte 5, **all sequence scanning, dipeptide feature extraction, and RBF kernel classification run locally in your browser**—ensuring instantaneous results with zero server queues or data privacy concerns.

---

## 🧬 Biological Context & Why LBtope?

Predicting B-cell epitopes (antigenic regions) with high accuracy is one of the primary challenges in designing subunit/peptide vaccines and immunotherapies. Earlier computational tools (like ABCpred) suffered from two major limitations:
1. They were trained on small datasets.
2. They relied on random peptides as assumed non-epitopes.

**LBtope** overcame these limitations by leveraging Support Vector Machine (SVM) models trained on a massive dataset from the [IEDB database](https://www.iedb.org/)—featuring **12,063 experimentally validated epitopes** and **20,589 experimentally validated non-epitopes**. For the first time, confirmed non-B-cell epitopes were utilized as negative controls to build highly selective prediction models with an overall accuracy of **~81%**.

---

## ✨ Key Features & Modules

### 1. Antigen Sequence Prediction
* **Sliding Window Scanning:** Evaluates full-length protein/antigen sequences submitted in FASTA format.
* **1-to-1 Centered Mapping:** Each amino acid is assigned the prediction score of the sliding window centered directly on it (padded with `'X'` residues at sequence boundaries).
* **Interactive Visualization:** Sequence residues are color-coded by probability tier (0–20% Black, 21–40% Orange, 41–60% Blue, 61–80% Green, 81–100% Red).
* **Pinned Popovers & Action Bar:** Click any amino acid to lock its prediction data in place. Seamlessly transfer any selected window directly into the **Peptide Mutants** tool for rational design.
* **Data Export:** Download detailed window-by-window scores as TSV/text files.

### 2. Multiple Peptides Batch Classification
* **High-Speed Batch Scoring:** Submit lists of discrete peptides of varying lengths simultaneously in FASTA format.
* **Instant Classification:** Evaluates dipeptide composition features without sliding windows, categorizing each sequence as an **Epitope** or **Non-epitope**.
* **Dynamic Thresholding:** Highlight results exceeding customizable probability thresholds (20% to 80%).

### 3. Peptide Mutants (Rational Vaccine Design)
* **Systematic Point-Mutation Scan:** Performs a complete $L \times 20$ substitution matrix scan across a submitted peptide, evaluating every possible single-amino-acid replacement.
* **Heatmap Grid View:** A 2D matrix comparing sequence positions against all 20 amino acids, color-coded by $\Delta\text{Probability}$ (green for binding enhancements, red for degradations).
* **Ranked Enhancements Table:** Sorts substitutions by the highest gain in epitope affinity over the wild-type sequence.

### 4. Algorithm & Help Documentation
* **Embedded Methodology:** Built-in interactive documentation detailing the 5 core datasets, legacy pipeline schematics, and interpretation guides.

---

## 🧠 The Client-Side ML Pipeline

Instead of relying on a Perl/PHP backend calling a C-compiled `svm_classify` binary, this Svelte app implements the entire classification pipeline in JavaScript:
[FASTA Input] ➔ [extractDPCFeatures() (400 DPC %)] ➔ [classifySVM() (RBF Kernel)] ➔ [scoreToProbability() (%)]
### 1. Dipeptide Composition (DPC) Extraction (`extractDPCFeatures`)
For any peptide sequence of length $L$, the algorithm computes the percentage composition of all 400 ordered amino acid pairs ($\text{AA}, \text{AC}, \dots, \text{YY}$):
$$\text{Feature}_{i,j} = \frac{\text{Count}(AA_i, AA_j) \times 100}{L - 1}$$
* **Legacy Perl Fidelity:** To match the exact behavior of the original `pro2dpc` Perl script, padding `'X'` characters and unrecognized symbols map to index `-1`, which resolves to **Tyrosine (`'Y'`)** in feature vector weighting.
* Precision is rounded to 3 decimal places (`%5.3f`) before classification.

### 2. RBF Support Vector Machine Classifier (`classifySVM`)
The raw SVM decision score is calculated over the stored support vectors using a Radial Basis Function (RBF) kernel:
$$\text{Score} = \left[ \sum_{k=1}^{N_{sv}} \alpha_k \cdot \exp\left(-\gamma \Vert{}x - sv_k\Vert{}^2\right) \right] - \text{bias}$$
To maintain precision across sparse feature representations without dropping query magnitude terms, Euclidean distances are computed via the quadratic decomposition:
$$\Vert{}x - sv\Vert{}^2 = \Vert{}x\Vert{}^2 - 2(x \cdot sv) + \Vert{}sv\Vert{}^2$$

### 3. Binary Model Architecture (`.bin` format)
Legacy text models (`model.dat`) are pre-converted via a custom Python script (`convertsvmtobin.py`) into a compact, little-endian binary format served statically from `public/models/`:
* **Header (24 Bytes):** `uint32` kernelType, `float64` gamma, `float64` bias, `uint32` numSV.
* **Per Support Vector:** `float64` $\alpha_y$ weight, `uint16` active feature count ($N_f$).
* **Per Active Feature:** `uint16` feature index (1–400), `float64` value (near-zero floating-point noise $< 10^{-9}$ is stripped during conversion).

---

## 📦 Available Models & Datasets

The application serves 5 distinct models optimized for different predictive needs:

| Model File | Website Label | Window Size | Training Dataset Description |
| :--- | :--- | :--- | :--- |
| `model-20mer-complete.bin` | **LBtope_Fixed** | 20 (Fixed) | 12,063 positive & 20,589 negative 20-mer patterns. |
| `model-20mer-non-variant.bin` | **LBtope_Fixed_non_redundant** | 20 (Fixed) | CD-HIT 80% non-redundant dataset (7,824 pos / 7,853 neg). |
| `model-final-dpc-flex-svm.bin` | **LBtope_Variable** | 5–30 (User selected, def: 15) | 14,876 positive & 23,321 negative variable-length peptides. |
| `model-variable-non-redundant.bin` | **LBtope_Variable_non_redundant** | 5–30 (User selected, def: 15) | CD-HIT 80% non-redundant variable-length dataset. |
| `model-final-dpc-flex-more2-svm.bin` | **LBtope_Confirm** | 5–30 (User selected, def: 15) | Gold-standard set confirmed by $\ge 2$ independent studies. |

---

## 🚀 Local Development & Setup

### Prerequisites
* [Node.js](https://nodejs.org/) (v18 or higher recommended)
* [npm](https://www.npmjs.com/) or [pnpm](https://pnpm.io/)

### Installation
1. **Clone the repository:**
   ```bash
   git clone [https://github.com/your-username/lbtope.git](https://github.com/your-username/lbtope.git)
   cd lbtope
Install dependencies:
  npm install
Start the local development server:
  npm run dev
Open your browser to http://localhost:5173.

Build for Production
To generate an optimized static production build:
  npm run build
The compiled assets will be output to the dist/ folder.

🌐 GitHub Pages Deployment
This repository is configured to deploy seamlessly as a static site to GitHub Pages. In vite.config.js, the base public path is set dynamically:

JavaScript
export default defineConfig(({ command }) => ({
  base: command === 'build' ? '/lbtope/' : '/',
  // ...
}));
Manual Deployment
If you use the gh-pages branch deployment method:
  npm run build
  npx gh-pages -d dist
Automated GitHub Actions Deployment
If deploying via GitHub Actions, ensure your repository Settings ➔ Pages is set to deploy from GitHub Actions, and include the standard Vite static workflow pointing to your dist artifact.

📁 Project Structure
Plaintext
lbtope/
├── public/
│   ├── images/                # Legacy server logos, flowcharts, & help screenshots
│   │   ├── lbtope1.gif
│   │   ├── algotithm1.jpg ...
│   │   └── result_protein.jpg ...
│   └── models/                # Pre-converted binary SVM models (.bin)
│       ├── model-final-dpc-flex-svm.bin
│       └── ...
├── src/
│   ├── lib/
│   │   ├── AntigenPredictor.svelte    # Sliding-window protein scanner & popovers
│   │   ├── PeptidePredictor.svelte    # Multiple peptide batch classification table
│   │   ├── MutantPredictor.svelte     # L×20 substitution heatmap & ranker
│   │   ├── AlgorithmView.svelte       # Dataset methodology & flowcharts
│   │   ├── HelpView.svelte            # User manual & legacy references
│   │   └── SequenceStrip.svelte       # Header visual ornament
│   ├── App.svelte                     # Root application, navigation & routing state
│   ├── main.js                        # Svelte app mounting entry point
│   └── app.css                        # Global design system & legacy green theme
├── convertsvmtobin.py                 # Python utility to convert SVM-light text -> .bin
├── index.html                         # HTML root
├── package.json
└── vite.config.js                     # Vite bundling & GitHub Pages base config
🔬 Citation & References
If you use LBtope or this modernized software in your research, please cite the original publication:

Singh H, Ansari HR, Raghava GPS (2013). Improved Method for Linear B-Cell Epitope Prediction Using Antigen's Primary Sequence. PLoS ONE 8(5): e62216. https://doi.org/10.1371/journal.pone.0062216

Original Server & Datasets: IMTECH / OSDD Raghava Lab

IEDB Database: Immune Epitope Database and Analysis Resource

📄 License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL v3) as published by the Free Software Foundation.