# clinker web application

Extract the archive and serve this directory through an HTTP server or static
host. For a local preview, with Python installed:

```sh
python3 -m http.server 8000
```

Open `http://localhost:8000`. Opening `index.html` directly does not work:
JavaScript modules, Web Workers and WebAssembly require HTTP serving.

For public hosting, use HTTPS and serve `.js` and `.mjs` as JavaScript and
`.wasm` as `application/wasm`. Keep the directory structure intact, including
`pkg/` and `examples/`. No backend service, API key or database is required.
The application loads its assets from this directory and processes input files
locally in the browser.

Use a current browser with WebAssembly SIMD and module-worker support. SVG
rendering remains available when GPU rendering is unavailable. The automated
browser checks cover Chromium, Firefox and WebKit. A final check in actual
Safari remains useful because Playwright WebKit is not the Safari application.

Save work with **Export → Project** and reopen it with **Open project**.
SVG and Similarity CSV are output formats, not restorable projects. Python
session files and standalone CLI analysis sessions cannot be opened here.

The archive includes the clinker, clustermap.js and D3 licence notices.
