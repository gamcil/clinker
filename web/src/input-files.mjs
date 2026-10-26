export const isGenBankFile = file => /\.(?:gb|gbk|genbank|gbf|gbff)$/i.test(file.name);
export const isGffFile = file => /\.(?:gff|gff3)$/i.test(file.name);
export const isFastaFile = file => /\.(?:fa|fsa|fna|fasta|faa)$/i.test(file.name);
export const isSupportedInputFile = file => isGenBankFile(file) || isGffFile(file) || isFastaFile(file);
export const fileStem = name => name.replace(/\.[^.]+$/, "").toLowerCase();

// A GFF annotation and its sequence share one input-file row. Keep their
// filenames together so displaying and removing a row use the same pairing.
function classifyInputFiles(files, metadata) {
  const fastaByStem = new Map();
  const annotationStems = new Set(files.filter(isGffFile).map(file => fileStem(file.name)));
  for (const file of files.filter(isFastaFile)) {
    const stem = fileStem(file.name);
    if (!fastaByStem.has(stem)) fastaByStem.set(stem, []);
    fastaByStem.get(stem).push(file);
  }

  const rows = files.filter(file => !isFastaFile(file)).map(file => {
    const row = { name: file.name, files: [file.name], detail: "Ready", ready: true };
    const companions = isGffFile(file) ? fastaByStem.get(fileStem(file.name)) || [] : [];
    row.files.push(...companions.map(companion => companion.name));
    const unreadable = [file, ...companions].find(input => metadata.get(input)?.error);
    if (unreadable) {
      return { ...row, detail: "Unreadable", reason: `Could not read ${unreadable.name}: ${metadata.get(unreadable).error}`, ready: false };
    }
    if (isGenBankFile(file)) return { ...row, description: "GenBank" };
    if (!isGffFile(file)) {
      return { ...row, detail: "Unsupported", reason: "Choose GenBank or GFF3 with FASTA.", ready: false };
    }
    row.description = companions.length ? `+ ${companions.map(companion => companion.name).join(", ")}` : "GFF3";
    if (companions.length > 1) {
      return { ...row, detail: "Multiple FASTA files", reason: "Keep one matching FASTA file for this annotation.", ready: false };
    }
    if (companions.length === 1) return row;
    if (metadata.get(file).embeddedFasta) {
      return { ...row, description: "GFF3 with embedded FASTA" };
    }
    return { ...row, detail: "Missing FASTA", reason: "Add a same-name FASTA file, or use GFF3 with embedded FASTA.", ready: false };
  });

  for (const file of files.filter(isFastaFile)) {
    if (!annotationStems.has(fileStem(file.name))) {
      rows.push({ name: file.name, files: [file.name], detail: "Missing GFF3", reason: "Add a same-name GFF3 annotation file.", ready: false });
    }
  }
  return { rows, ready: rows.length > 0 && rows.every(row => row.ready) };
}

function genBankRecords(text) {
  const records = [];
  const locusLine = /^LOCUS\s+(\S+)\s+(\d+)\s+bp(?:\s|$)/gim;
  const matches = [...text.matchAll(locusLine)];
  for (let index = 0; index < matches.length; index++) {
    const match = matches[index];
    const length = Number(match[2]);
    if (Number.isSafeInteger(length) && length > 0) {
      const record = { name: match[1], start: 1, end: length, length };
      const header = text.slice(match.index, matches[index + 1]?.index ?? text.length);
      const organism = header.match(/^\s+ORGANISM\s+([^\r\n]+)/m)?.[1].trim()
        || header.match(/^SOURCE\s+([^\r\n]+)/m)?.[1].trim();
      if (organism) record.organism = organism;
      records.push(record);
    }
  }
  return records;
}

function fastaRecords(text) {
  const records = [];
  let current = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (current) records.push(current);
      current = { name: line.slice(1).trim().split(/\s+/, 1)[0], start: 1, end: 0, length: 0 };
    } else if (current) {
      current.length += line.replace(/\s/g, "").length;
    }
  }
  if (current) records.push(current);
  return records
    .filter(record => record.name && record.length > 0)
    .map(record => ({ ...record, end: record.length }));
}

function gffRecords(text) {
  const records = [...text.matchAll(/^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)\s*$/gim)]
    .map(match => ({
      name: match[1],
      start: Number(match[2]),
      end: Number(match[3]),
      length: Number(match[3]) - Number(match[2]) + 1,
    }))
    .filter(record => Number.isSafeInteger(record.start)
      && Number.isSafeInteger(record.end)
      && Number.isSafeInteger(record.length)
      && record.start > 0
      && record.end >= record.start);
  if (records.length > 0) return records;
  const embedded = text.match(/^##FASTA\s*$([\s\S]*)/im)?.[1] || "";
  return fastaRecords(embedded);
}

// Cache only sequence metadata, never complete file contents. Read sequentially
// so selecting many large files does not retain all their decoded text at once.
export async function inspectInputFiles(files, { readText = file => file.text(), cache = new WeakMap() } = {}) {
  const metadata = new Map();
  for (const file of files.filter(isSupportedInputFile)) {
    let info = cache.get(file);
    if (!info) {
      try {
        const text = await readText(file);
        info = {
          records: isGenBankFile(file) ? genBankRecords(text)
            : isGffFile(file) ? gffRecords(text) : fastaRecords(text),
          embeddedFasta: isGffFile(file) && /^##FASTA[\t ]*$/im.test(text),
        };
        cache.set(file, info);
      } catch (error) {
        if (error.name === "AbortError") throw error;
        info = { records: [], error: error.message || String(error) };
      }
    }
    metadata.set(file, info);
  }
  const fastaByStem = new Map(files.filter(isFastaFile).map(file => [fileStem(file.name), file]));
  const filesWithRecords = files.filter(file => isGenBankFile(file) || isGffFile(file)).map(file => {
    let records = metadata.get(file).records;
    if (isGffFile(file) && records.length === 0) {
      records = metadata.get(fastaByStem.get(fileStem(file.name)))?.records || [];
    }
    return { file, records };
  });
  return { ...classifyInputFiles(files, metadata), filesWithRecords };
}
