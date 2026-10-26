// One user operation owns its reads, requests and workers. Cancelling rejects
// awaited work immediately and terminates synchronous WASM in its workers.
export class Workflow {
  #controller = new AbortController();
  #workers = new Set();

  get signal() { return this.#controller.signal; }

  throwIfCancelled() { this.signal.throwIfAborted(); }

  cancel() {
    this.dispose();
  }

  dispose() {
    this.#controller.abort(new DOMException("Cancelled", "AbortError"));
    for (const worker of this.#workers) this.releaseWorker(worker);
  }

  worker(url) {
    this.throwIfCancelled();
    const worker = new Worker(url, { type: "module" });
    this.#workers.add(worker);
    return worker;
  }

  releaseWorker(worker) {
    worker.onmessage = worker.onerror = worker.onmessageerror = null;
    worker.terminate();
    this.#workers.delete(worker);
  }

  request(url, message, onProgress = () => {}) {
    const worker = this.worker(url);
    return new Promise((resolve, reject) => {
      let settled = false;
      const finish = (callback, value) => {
        if (settled) return;
        settled = true;
        this.signal.removeEventListener("abort", cancel);
        this.releaseWorker(worker);
        callback(value);
      };
      const cancel = () => finish(reject, this.signal.reason);
      this.signal.addEventListener("abort", cancel, { once: true });
      worker.onmessage = ({ data }) => {
        if (settled) return;
        try {
          if (data.type === "error") finish(reject, new Error(data.message));
          else if (data.type === "post-process") finish(resolve, data.result);
          else if (data.type === "post-process-progress") onProgress(data);
        } catch (error) { finish(reject, error); }
      };
      worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
      worker.onmessageerror = () => finish(reject, new Error("Could not read the analysis worker response."));
      try { worker.postMessage(message); }
      catch (error) { finish(reject, error); }
    });
  }

  wait(promise) {
    return new Promise((resolve, reject) => {
      const abort = () => reject(this.signal.reason);
      const finish = (callback, value) => {
        this.signal.removeEventListener("abort", abort);
        callback(value);
      };
      // Observe the underlying promise even after cancellation so a later
      // failure cannot become an unhandled rejection.
      Promise.resolve(promise).then(
        value => finish(resolve, value),
        error => finish(reject, error),
      );
      if (this.signal.aborted) abort();
      else this.signal.addEventListener("abort", abort, { once: true });
    });
  }

  readFile(file, asText = false) {
    this.throwIfCancelled();
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      const abort = () => reader.abort();
      const finish = (callback, value) => {
        this.signal.removeEventListener("abort", abort);
        callback(value);
      };
      reader.onload = () => finish(resolve, reader.result);
      reader.onerror = () => finish(reject, reader.error);
      reader.onabort = () => finish(reject, this.signal.reason);
      this.signal.addEventListener("abort", abort, { once: true });
      if (asText) reader.readAsText(file);
      else reader.readAsArrayBuffer(file);
    });
  }
}
