import { syncEmbedded } from './assets.mjs';

const args = process.argv.slice(2);
if (args.some(arg => arg !== '--check')) throw new Error('Usage: sync-assets.mjs [--check]');
await syncEmbedded({ check: args.includes('--check') });
console.log(args.includes('--check') ? 'Embedded assets match pinned dependencies.' : 'Embedded assets synchronized.');
