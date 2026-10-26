"""Thin timing adapter for an explicit legacy Python source checkout."""
import io
import json
from pathlib import Path
import sys
from time import perf_counter



def main():
    request = json.load(sys.stdin)
    source = Path(request["legacy_source"])
    sys.path.insert(0, str(source))
    from Bio import SeqIO
    from clinker.classes import Cluster, Locus
    from clinker.align import align_clusters
    staged = [(Path(path).stem, Path(path).read_bytes()) for path in request["files"]]
    start = perf_counter()
    # Same construction as cluster_from_genbank, reading staged bytes rather
    # than reopening files: input acquisition is outside all three timers.
    clusters = [Cluster(name, [Locus.from_seqrecord(record, set_origin=True)
                              for record in SeqIO.parse(io.StringIO(raw.decode("utf-8")), "genbank")])
                for name, raw in staged]
    parse_ms = (perf_counter() - start) * 1000
    if request["inspect"]:
        print(json.dumps({"proteins": [[g.translation.strip() for locus in c.loci
                                       for g in locus.genes] for c in clusters]}))
        return
    align_start = perf_counter()
    result = align_clusters(*clusters, cutoff=request["identity"], jobs=request["workers"])
    alignment_ms = (perf_counter() - align_start) * 1000
    total_ms = (perf_counter() - start) * 1000
    counts = [sum(bool(g.translation.strip()) for locus in c.loci for g in locus.genes)
              for c in clusters]
    possible = sum(a * b for i, a in enumerate(counts) for b in counts[i + 1:])
    print(json.dumps({"parse_ms": parse_ms, "prepare_ms": 0, "prefilter_ms": None,
                     "alignment_ms": alignment_ms, "total_ms": total_ms,
                     "clusters": len(clusters), "proteins": sum(counts),
                     "amino_acids": sum(len(g.translation.strip()) for c in clusters
                                        for locus in c.loci for g in locus.genes),
                     "possible_pairs": possible, "eligible_pairs": None,
                     "candidate_pairs": possible, "aligned_pairs": possible,
                     "retained_links": sum(len(a.links) for a in result.alignments.values())}))


if __name__ == "__main__":
    main()
