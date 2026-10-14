// One user operation owns its reads, requests and workers. Cancelling rejects
// awaited work immediately and terminates synchronous WASM in its workers.
export class Workflow {
  #controller = new AbortController();
  #workers = new Set();

  get signal() { return this.#controller.signal; }

  throwIfCancelled() { this.signal.throwIfAborted(); }

  cancel() {
    this.#controller.abort(new DOMException("Cancelled", "AbortError"));
    this.dispose();
  }

  dispose() {
    for (const worker of this.#workers) worker.terminate();
    this.#workers.clear();
  }

  worker(url) {
    this.throwIfCancelled();
    const worker = new Worker(url, { type: "module" });
    this.#workers.add(worker);
    return worker;
  }

  releaseWorker(worker) {
    worker.terminate();
    this.#workers.delete(worker);
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
