def get_maximum_row_lengths(rows):
    """Finds the longest lengths of fields per column in a collection of rows."""
    lengths, total = [], len(rows[0])
    for index in range(total):
        largest = max(len(str(row[index])) for row in rows)
        lengths.append(largest)
    return lengths


def add_field_whitespace(rows, lengths):
    """Fills table fields with whitespace to specified lengths."""
    result = []
    for row in rows:
        fmt = [f"{row[index]:{length}}" for index, length in enumerate(lengths)]
        result.append(fmt)
    return result


def humanise(rows):
    """Formats a collection of fields as human-readable."""
    lengths = get_maximum_row_lengths(rows)
    table = add_field_whitespace(rows, lengths)
    return table


def get_link_values(link, decimals=4):
    return [
        link.query.label,
        link.target.label,
        f"{link.identity:.{decimals}f}",
        f"{link.similarity:.{decimals}f}",
    ]


def format_links(
    links,
    decimals=4,
    headers=False,
    delimiter=None,
):
    """Generates a summary table for a hit cluster.

    Args:
        links (list): collection of Link objects
        decimals (int): number of decimal points to show
        show_headers (bool): show column headers in output
        human (bool): use human-readable format
    Returns:
        summary table
    """
    rows = [
        get_link_values(link, decimals=decimals)
        for link in links
    ]
    if headers:
        hdrs = ["Query", "Target", "Identity", "Similarity"]
        rows.insert(0, hdrs)
    if not delimiter:
        delimiter = "  "
        rows = humanise(rows)
    return "\n".join(delimiter.join(row) for row in rows)


def format_alignment(
    alignment,
    decimals=4,
    alignment_headers=True,
    link_headers=True,
    delimiter=None,
):
    fmt = format_links(
        alignment.links,
        decimals=decimals,
        headers=link_headers,
        delimiter=delimiter,
    )
    if alignment_headers:
        header = f"{alignment.query.name} vs {alignment.target.name}"
        separator = "-" * len(header)
        fmt = f"{header}\n{separator}\n{fmt}"
    return fmt


def format_globaligner(
    aligner,
    decimals=4,
    alignment_headers=True,
    link_headers=True,
    delimiter=None,
):
    if not aligner.alignments:
        raise ValueError("No alignments are stored in aligner")
    fmts = []
    for alignment in aligner.alignments.values():
        fmt = format_alignment(
            alignment,
            alignment_headers=alignment_headers,
            link_headers=link_headers,
            delimiter=delimiter,
            decimals=decimals,
        )
        fmts.append(fmt)
    return "\n\n".join(fmts)
