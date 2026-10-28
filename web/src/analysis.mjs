const WORKER_URL = new URL("./worker.js", import.meta.url);

export function parseInputFiles(
  files,
  settings,
  { onProgress, onPrefilterStart, onPrefilterProgress },
  workflow,
) {
  const worker = workflow.worker(WORKER_URL);
  return new Promise((resolve, reject) => {
    let nextFileIndex = 0;
    let settled = false;
    const startupTimeout = setTimeout(() => {
      finish(reject, new Error("The analysis worker did not start. Reload the page and try again."));
    }, 30000);
    const finish = (callback, value) => {
      if (settled) return;
      settled = true;
      clearTimeout(startupTimeout);
      workflow.signal.removeEventListener("abort", cancel);
      workflow.releaseWorker(worker);
      callback(value);
    };
    const cancel = () => finish(reject, workflow.signal.reason);
    workflow.signal.addEventListener("abort", cancel, { once: true });
    const sendNextFile = () => {
      if (settled) return;
      try {
        if (nextFileIndex === files.length) {
          if (settings.enabled) onPrefilterStart();
          worker.postMessage({ type: "parse-finish" });
          return;
        }
        const fileIndex = nextFileIndex;
        const file = files[nextFileIndex++];
        // Transfer each input buffer once; it belongs to this run.
        worker.postMessage(
          { type: "parse-file", file, fileIndex },
          [file.bytes.buffer],
        );
      } catch (error) { finish(reject, error); }
    };
    worker.onmessage = ({ data }) => {
      if (settled) return;
      try {
        if (data.type === "error") finish(reject, new Error(data.message));
        else if (data.type === "parse-ready") {
          clearTimeout(startupTimeout);
          sendNextFile();
        }
        else if (data.type === "file-parsed") {
          onProgress(data.fileIndex + 1, files.length);
          setTimeout(sendNextFile, 0);
        } else if (data.type === "filter-progress") {
          onPrefilterProgress(data.completed, data.total);
        } else if (data.type === "parsed") {
          finish(resolve, { ...data.parsed, wasmModule: data.wasmModule });
        }
      } catch (error) { finish(reject, error); }
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.onmessageerror = () => finish(reject, new Error("Could not read the analysis worker response."));
    try { worker.postMessage({ type: "parse-start", prefilter: settings }); }
    catch (error) { finish(reject, error); }
  });
}

export function openAnalysisSession(session, workflow) {
  return workflow.request(WORKER_URL, { type: "open-session", session });
}

export function postProcess({ links, ...options }, onProgress, workflow) {
  return workflow.request(WORKER_URL, {
    type: "post-process",
    ...options,
    reportProgress: true,
    links: links.map(({ query, target, identity, similarity }) => ({ query, target, identity, similarity })),
  }, onProgress);
}
