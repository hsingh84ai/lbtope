<script>
  import { onMount, onDestroy } from 'svelte';

  // Hen egg-white lysozyme (a standard reference sequence) — real biological
  // data rather than placeholder text, fitting for the tool's own domain.
  const sequence =
    'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL';
  const windowSize = 15;

  let start = 0;
  let timer;

  const prefersReducedMotion =
    typeof window !== 'undefined' &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  onMount(() => {
    if (prefersReducedMotion) return;
    timer = setInterval(() => {
      start = (start + 1) % (sequence.length - windowSize);
    }, 550);
  });

  onDestroy(() => timer && clearInterval(timer));
</script>

<div class="strip" aria-hidden="true">
  <div class="seq">
    {#each sequence.split('') as residue, i}
      <span class="res" class:hit={i >= start && i < start + windowSize}>{residue}</span>
    {/each}
  </div>
  <div class="caption">15-residue scoring window, sliding across the sequence</div>
</div>

<style>
  .strip {
    width: 100%;
    overflow: hidden;
  }

  .seq {
    font-family: var(--font-mono);
    font-size: 0.95rem;
    line-height: 1.9;
    letter-spacing: 0.02em;
    word-break: break-all;
    color: var(--color-quiet);
  }

  .res {
    display: inline-block;
    transition: color 0.15s ease, background 0.15s ease;
    border-radius: 2px;
  }

  .res.hit {
    color: var(--color-ink);
    background: var(--color-epitope-bg);
    box-shadow: 0 0 0 1px var(--color-epitope);
  }

  .caption {
    margin-top: 0.6rem;
    font-family: var(--font-mono);
    font-size: 0.75rem;
    color: var(--color-quiet);
  }

  @media (max-width: 640px) {
    .seq {
      font-size: 0.8rem;
      line-height: 1.7;
    }
  }
</style>
