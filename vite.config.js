import { defineConfig } from 'vite';
import { svelte } from '@sveltejs/vite-plugin-svelte';

export default defineConfig(({ command }) => ({
  plugins: [svelte()],
  // Only applied for `npm run build` (GitHub Pages) — dev server stays at '/'
  base: command === 'build' ? '/lbtope/' : '/',
  build: {
    target: 'es2020'
  }
}));