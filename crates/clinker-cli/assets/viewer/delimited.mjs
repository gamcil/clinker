export function parseDelimitedRows(text, delimiter) {
  const source = String(text).replace(/^\uFEFF/, "");
  if (delimiter !== "," && delimiter !== "\t") throw new Error("Choose comma or tab as the separator.");
  const rows = [];
  let row = [], field = "", quoted = false, closedQuote = false;
  for (let i = 0; i < source.length; i += 1) {
    const char = source[i];
    if (quoted) {
      if (char === '"' && source[i + 1] === '"') { field += '"'; i += 1; }
      else if (char === '"') { quoted = false; closedQuote = true; }
      else field += char;
      continue;
    }
    if (char === '"') {
      if (field.length || closedQuote) throw new Error(`Unexpected quote near character ${i + 1}.`);
      quoted = true;
    } else if (char === delimiter) {
      row.push(field); field = ""; closedQuote = false;
    } else if (char === "\n" || char === "\r") {
      if (char === "\r" && source[i + 1] === "\n") i += 1;
      row.push(field); field = ""; closedQuote = false;
      if (row.some(value => value.trim() !== "")) rows.push(row);
      row = [];
    } else {
      if (closedQuote && !/\s/.test(char)) throw new Error(`Unexpected text after a quoted field near character ${i + 1}.`);
      if (!closedQuote) field += char;
    }
  }
  if (quoted) throw new Error("The table contains an unclosed quoted field.");
  row.push(field);
  if (row.some(value => value.trim() !== "")) rows.push(row);
  return rows;
}

export function parseTwoColumnCsv(text, description) {
  return parseDelimitedRows(text, ",").map((fields, index) => {
    const row = fields.map(value => value.trim());
    if (row.length !== 2 || row.some(value => !value)) {
      throw new Error(`${description} row ${index + 1} must contain two non-empty columns.`);
    }
    return row;
  });
}
