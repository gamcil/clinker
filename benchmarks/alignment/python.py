import json
import sys
import time
from pathlib import Path

from clinker.align import Globaligner
from clinker.classes import parse_files


def main():
    dataset, identity, jobs, output = sys.argv[1:]
    paths = sorted(str(path) for path in Path(dataset).glob("*.gbk"))
    clusters = parse_files(paths)
    aligner = Globaligner()
    aligner.add_clusters(*clusters)
    started = time.perf_counter()
    aligner.align_stored_clusters(cutoff=float(identity), jobs=int(jobs))
    elapsed = time.perf_counter() - started
    links = sum(len(alignment.links) for alignment in aligner.alignments.values())
    Path(output).write_text(json.dumps({"elapsed_seconds": elapsed, "retained_links": links}) + "\n")


if __name__ == "__main__":
    main()
