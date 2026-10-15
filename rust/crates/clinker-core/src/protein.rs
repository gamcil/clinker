pub(crate) fn normalise_protein(protein: &[u8], output: &mut Vec<u8>) {
    output.clear();
    output.extend(protein.iter().map(|&amino_acid| {
        let amino_acid = amino_acid.to_ascii_uppercase();
        if b"ARNDCQEGHILKMFPSTWYVBZX".contains(&amino_acid) {
            amino_acid
        } else {
            b'X'
        }
    }));
}

pub(crate) fn normalise_translation(translation: &str) -> String {
    translation
        .split_whitespace()
        .collect::<String>()
        .trim_end_matches('*')
        .to_owned()
}

pub(crate) fn translate_standard(coding_sequence: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let index = |base: u8| match base.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    let translation = coding_sequence
        .as_chunks::<3>()
        .0
        .iter()
        .map(
            |codon| match (index(codon[0]), index(codon[1]), index(codon[2])) {
                (Some(a), Some(b), Some(c)) => TABLE[a * 16 + b * 4 + c] as char,
                _ => 'X',
            },
        )
        .collect::<String>();
    normalise_translation(&translation)
}
