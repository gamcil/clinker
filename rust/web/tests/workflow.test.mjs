import test from "node:test";
import assert from "node:assert/strict";
import { Workflow } from "../workflow.mjs";

test("cancellation interrupts preparation even if the underlying promise never settles", async () => {
  const workflow = new Workflow();
  const waiting = workflow.wait(new Promise(() => {}));
  workflow.cancel();
  await assert.rejects(waiting, { name: "AbortError" });
  assert.throws(() => workflow.throwIfCancelled(), { name: "AbortError" });
});

test("a cancelled wait still observes failures from abandoned work", async () => {
  const workflow = new Workflow();
  let fail;
  const waiting = workflow.wait(new Promise((_, reject) => { fail = reject; }));
  workflow.cancel();
  await assert.rejects(waiting, { name: "AbortError" });
  fail(new Error("late failure"));
  await new Promise(resolve => setImmediate(resolve));
});

test("cancelling disposes every owned worker, and a new workflow starts independently", t => {
  class FakeWorker {
    terminations = 0;
    terminate() { this.terminations += 1; }
  }
  const originalWorker = globalThis.Worker;
  globalThis.Worker = FakeWorker;
  t.after(() => {
    if (originalWorker) globalThis.Worker = originalWorker;
    else delete globalThis.Worker;
  });
  const workflow = new Workflow();
  const finished = workflow.worker("worker.js");
  const active = workflow.worker("worker.js");
  workflow.releaseWorker(finished);
  workflow.cancel();
  assert.equal(finished.terminations, 1);
  assert.equal(active.terminations, 1);
  assert.throws(() => workflow.worker("worker.js"), { name: "AbortError" });
  const retry = new Workflow();
  assert.equal(retry.signal.aborted, false);
  retry.dispose();
});

test("cancelling aborts an in-progress file read", async t => {
  let reader;
  class FakeFileReader {
    constructor() { reader = this; }
    readAsArrayBuffer() { this.reading = true; }
    abort() { this.reading = false; this.onabort(); }
  }
  const originalReader = globalThis.FileReader;
  globalThis.FileReader = FakeFileReader;
  t.after(() => {
    if (originalReader) globalThis.FileReader = originalReader;
    else delete globalThis.FileReader;
  });
  const workflow = new Workflow();
  const read = workflow.readFile({});
  assert.equal(reader.reading, true);
  workflow.cancel();
  await assert.rejects(read, { name: "AbortError" });
  assert.equal(reader.reading, false);
});
