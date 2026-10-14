export const isGenBankFile = file => /\.(?:gb|gbk|genbank|gbf|gbff)$/i.test(file.name);
export const isGffFile = file => /\.(?:gff|gff3)$/i.test(file.name);
export const isFastaFile = file => /\.(?:fa|fsa|fna|fasta|faa)$/i.test(file.name);
export const isSupportedInputFile = file => isGenBankFile(file) || isGffFile(file) || isFastaFile(file);
export const fileStem = name => name.replace(/\.[^.]+$/, "").toLowerCase();

// A GFF annotation and its sequence share one input-file row. Keep their
// filenames together so displaying and removing a row use the same pairing.
export async function classifyInputFiles(files) {
  const fastaByStem = new Map();
  const annotationStems = new Set(files.filter(isGffFile).map(file => fileStem(file.name)));
  for (const file of files.filter(isFastaFile)) {
    const stem = fileStem(file.name);
    if (!fastaByStem.has(stem)) fastaByStem.set(stem, []);
    fastaByStem.get(stem).push(file);
  }

  const rows = await Promise.all(files.filter(file => !isFastaFile(file)).map(async file => {
    const row = { name: file.name, files: [file.name], detail: "Ready", ready: true };
    if (isGenBankFile(file)) return { ...row, description: "GenBank" };
    if (!isGffFile(file)) {
      return { ...row, detail: "Unsupported", reason: "Choose GenBank or GFF3 with FASTA.", ready: false };
    }
    const companions = fastaByStem.get(fileStem(file.name)) || [];
    row.files.push(...companions.map(companion => companion.name));
    row.description = companions.length ? `+ ${companions.map(companion => companion.name).join(", ")}` : "GFF3";
    if (companions.length > 1) {
      return { ...row, detail: "Multiple FASTA files", reason: "Keep one matching FASTA file for this annotation.", ready: false };
    }
    if (companions.length === 1) return row;
    try {
      if (/^##FASTA[\t ]*$/im.test(await file.text())) {
        return { ...row, description: "GFF3 with embedded FASTA" };
      }
      return { ...row, detail: "Missing FASTA", reason: "Add a same-name FASTA file, or use GFF3 with embedded FASTA.", ready: false };
    } catch (error) {
      return { ...row, detail: "Unreadable", reason: `Could not read file: ${error.message || String(error)}`, ready: false };
    }
  }));

  for (const file of files.filter(isFastaFile)) {
    if (!annotationStems.has(fileStem(file.name))) {
      rows.push({ name: file.name, files: [file.name], detail: "Missing GFF3", reason: "Add a same-name GFF3 annotation file.", ready: false });
    }
  }
  return { rows, ready: rows.length > 0 && rows.every(row => row.ready) };
}
