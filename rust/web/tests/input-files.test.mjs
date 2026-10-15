import test from "node:test";
import assert from "node:assert/strict";
import { inspectInputFiles } from "../input-files.mjs";

const file = (name, text = "") => ({ name, text: async () => text });

test("GFF and companion FASTA share one removable input-file row, regardless of selection order", async () => {
  const validation = await inspectInputFiles([file("cluster.FNA"), file("cluster.gff3"), file("other.gbk")]);
  assert.equal(validation.ready, true);
  assert.equal(validation.rows.length, 2);
  assert.deepEqual(validation.rows[0].files, ["cluster.gff3", "cluster.FNA"]);
  assert.equal(validation.rows[0].description, "+ cluster.FNA");
});

test("missing and ambiguous companions remain attached to their input-file error", async () => {
  const missing = await inspectInputFiles([file("cluster.gff3")]);
  assert.equal(missing.ready, false);
  assert.equal(missing.rows[0].detail, "Missing FASTA");
  const ambiguous = await inspectInputFiles([file("cluster.fa"), file("cluster.gff"), file("cluster.fna")]);
  assert.equal(ambiguous.ready, false);
  assert.equal(ambiguous.rows.length, 1);
  assert.deepEqual(ambiguous.rows[0].files, ["cluster.gff", "cluster.fa", "cluster.fna"]);
});

test("embedded FASTA needs no companion, but orphan FASTA blocks analysis", async () => {
  const embedded = await inspectInputFiles([file("cluster.gff3", "##gff-version 3\n##FASTA\n>chr\nATG")]);
  assert.equal(embedded.ready, true);
  assert.equal(embedded.rows[0].description, "GFF3 with embedded FASTA");
  const orphan = await inspectInputFiles([file("cluster.fa")]);
  assert.equal(orphan.ready, false);
  assert.equal(orphan.rows[0].detail, "Missing GFF3");
});

test("empty selection and unsupported files cannot start a run", async () => {
  assert.equal((await inspectInputFiles([])).ready, false);
  assert.equal((await inspectInputFiles([file("notes.txt")])).ready, false);
});

test("inspection joins companion metadata, preserves coordinates and reuses unchanged files", async () => {
  const files = [
    file("cluster.gff3", "##gff-version 3"),
    file("cluster.fna", ">scaffold_1 description\nATG\nAAA\n>scaffold_2\nCC"),
    file("other.gbk", "LOCUS       chromosome 42 bp DNA"),
    file("offset.gff3", "##sequence-region part 101 200"),
    file("offset.fa", ">part\nATG"),
  ];
  const cache = new WeakMap();
  const reads = [];
  let reading = false;
  const readText = async file => {
    assert.equal(reading, false, "file reads must not overlap");
    reading = true;
    reads.push(file.name);
    await new Promise(resolve => setImmediate(resolve));
    reading = false;
    return file.text();
  };
  const result = await inspectInputFiles(files, { readText, cache });
  assert.equal(result.ready, true);
  assert.deepEqual(result.filesWithRecords.map(({ records }) => records), [
    [{ name: "scaffold_1", start: 1, end: 6, length: 6 }, { name: "scaffold_2", start: 1, end: 2, length: 2 }],
    [{ name: "chromosome", start: 1, end: 42, length: 42 }],
    [{ name: "part", start: 101, end: 200, length: 100 }],
  ]);
  const added = file("added.gbk", "LOCUS       new 10 bp DNA");
  await inspectInputFiles([...files, added], { readText, cache });
  assert.deepEqual(reads, [...files, added].map(file => file.name));
});

test("unreadable companions block analysis, are retryable, and cancellation stops further reads", async () => {
  const annotation = file("cluster.gff3");
  let readable = false;
  const companion = { name: "cluster.fa", text: async () => {
    if (!readable) throw new Error("Disk unavailable");
    return ">chr\nATG";
  } };
  const cache = new WeakMap();
  const failed = await inspectInputFiles([annotation, companion], { cache });
  assert.equal(failed.ready, false);
  assert.match(failed.rows[0].reason, /cluster.fa: Disk unavailable/);
  assert.deepEqual(failed.rows[0].files, [annotation.name, companion.name]);
  readable = true;
  assert.equal((await inspectInputFiles([annotation, companion], { cache })).ready, true);
  let reads = 0;
  await assert.rejects(inspectInputFiles([annotation, companion], {
    readText: async () => { reads += 1; throw new DOMException("Cancelled", "AbortError"); },
  }), { name: "AbortError" });
  assert.equal(reads, 1);
});
