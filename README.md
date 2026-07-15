# LBTope Predictor (Svelte, static, serverless)

## What this is

A minimal Vite + Svelte app (no SvelteKit, no router) with two views toggled
by a single `view` variable in `App.svelte`:

- **home** — header nav + hero/center text (your request: menu with just
  "Predict Epitopes", and center text introducing the tool)
- **predict** — the epitope predictor, corrected from the prototype you
  uploaded (see "What changed" below)

## Local setup

```bash
npm install
npm run dev       # local dev server
npm run build     # outputs static site to dist/
npm run preview   # preview the production build locally
```

## Deploying to GitHub Pages

1. In `vite.config.js`, set `base: '/your-repo-name/'` to match your repo.
2. Build: `npm run build` — this produces `dist/`.
3. Either:
   - Push `dist/` to a `gh-pages` branch (e.g. via the `gh-pages` npm
     package), and enable Pages on that branch in repo settings, **or**
   - Add a GitHub Actions workflow that runs `npm run build` and deploys
     `dist/` on every push to `main` (recommended — no local build step to
     remember).

## The 76MB model file

This app keeps the prototype's "load model" file picker rather than
auto-fetching the `.bin` on page load, so you're not committing a 76MB
binary into git history (GitHub's hard limit is 100MB/file, but Pages will
still serve it slowly on first load for every visitor). If you want
auto-loading later:
- Re-check the sparsity of your support vectors — from the sample line
  you shared, most features per SV were zero, so a properly sparse binary
  encoding may shrink this well below 76MB.
- If you keep the file large, use Git LFS and `fetch()` it once on the
  predict view's mount, with a visible progress indicator.

## What changed from the prototype you uploaded

**Bug fix (correctness):** the original `predict15mer` only summed the
squared distance over each support vector's *stored* (nonzero) feature
indices. Any dimension where the query 15-mer was nonzero but a given SV
was zero (and therefore not stored in the sparse format) was silently
skipped — undercounting `||x - sv||²` for essentially every support
vector. Fixed via the standard decomposition:

```
||x - sv||² = ||x||² - 2(x·sv) + ||sv||²
```

`||x||²` is computed once per query window (dense, 400 dims). `||sv||²` is
precomputed once per support vector at model-load time. `x·sv` still only
touches each SV's nonzero entries, so this is both correct *and* no slower
than the buggy version.

**Performance/memory (not a behavior change):** support vectors are now
stored as typed arrays (structure-of-arrays: `Float64Array` for alphas,
`Uint16Array`/`Float32Array` for concatenated feature indices/values, with
an offset array marking each SV's slice) instead of ~31,000 individual JS
objects each holding an array of small feature objects. This avoids a
large number of small allocations when loading a model this size.

**Everything else** (DPC feature extraction, the 7×'X' padding, the
`-0.3` epitope threshold) was carried over unchanged from your prototype —
those encode your actual published method and I have no basis to second-guess
them. Worth double-checking the threshold against your real `lbtope.pl`
if predictions look off, since I only had the prototype's comment
claiming it matches line 50.

## Validating before you retire the live PHP tool

Run a batch of real sequences through both the old PHP/SVM-light pipeline
and this app, and diff the scores. The fix above should change some scores
(any window where the old bug undercounted distance) — that's expected,
not a regression.
