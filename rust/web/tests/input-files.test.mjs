import test from "node:test";
import assert from "node:assert/strict";
import { classifyInputFiles } from "../input-files.mjs";

const file = (name, text = "") => ({ name, text: async () => text });

test("GFF and companion FASTA share one removable input-file row, regardless of selection order", async () => {
  const validation = await classifyInputFiles([file("cluster.FNA"), file("cluster.gff3"), file("other.gbk")]);
  assert.equal(validation.ready, true);
  assert.equal(validation.rows.length, 2);
  assert.deepEqual(validation.rows[0].files, ["cluster.gff3", "cluster.FNA"]);
  assert.equal(validation.rows[0].description, "+ cluster.FNA");
});

test("missing and ambiguous companions remain attached to their input-file error", async () => {
  const missing = await classifyInputFiles([file("cluster.gff3")]);
  assert.equal(missing.ready, false);
  assert.equal(missing.rows[0].detail, "Missing FASTA");
  const ambiguous = await classifyInputFiles([file("cluster.fa"), file("cluster.gff"), file("cluster.fna")]);
  assert.equal(ambiguous.ready, false);
  assert.equal(ambiguous.rows.length, 1);
  assert.deepEqual(ambiguous.rows[0].files, ["cluster.gff", "cluster.fa", "cluster.fna"]);
});

test("embedded FASTA needs no companion, but orphan FASTA blocks analysis", async () => {
  const embedded = await classifyInputFiles([file("cluster.gff3", "##gff-version 3\n##FASTA\n>chr\nATG")]);
  assert.equal(embedded.ready, true);
  assert.equal(embedded.rows[0].description, "GFF3 with embedded FASTA");
  const orphan = await classifyInputFiles([file("cluster.fa")]);
  assert.equal(orphan.ready, false);
  assert.equal(orphan.rows[0].detail, "Missing GFF3");
});

test("empty selection and unsupported files cannot start a run", async () => {
  assert.equal((await classifyInputFiles([])).ready, false);
  assert.equal((await classifyInputFiles([file("notes.txt")])).ready, false);
});
