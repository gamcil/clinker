// The packaged editor uses D3's ESM entry for colour parsing. The application
// already loads the complete offline UMD bundle for clustermap, so expose the
// same function through a local module rather than shipping D3 twice.
export const color = globalThis.d3.color;
