"""Create benchmark copies without CDS fragments clipped outside a fetched region."""
import json
from pathlib import Path
import shutil
import sys

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition


def main():
    request = json.load(sys.stdin)
    output = Path(request["output"])
    datasets = []
    removed = []
    for source in request["datasets"]:
        source = Path(source)
        target = output / source.name
        target.mkdir(parents=True)
        datasets.append(str(target))
        for path in sorted(source.iterdir()):
            if path.suffix.lower() not in {".gbk", ".gb", ".gbff", ".genbank"}:
                continue
            records = list(SeqIO.parse(str(path), "genbank"))
            changed = False
            for record in records:
                keep = []
                for feature in record.features:
                    clipped = feature.type == "CDS" and any(
                        (isinstance(part.end, BeforePosition) and int(part.end) <= 1)
                        or (isinstance(part.start, AfterPosition) and int(part.start) >= len(record) - 1)
                        for part in feature.location.parts)
                    if clipped:
                        removed.append({"dataset": source.name, "file": path.name,
                                        "record": record.id, "location": str(feature.location),
                                        "protein_id": feature.qualifiers.get("protein_id"),
                                        "reason": "fuzzy CDS part extends outside fetched region at a boundary"})
                        changed = True
                    else:
                        keep.append(feature)
                record.features = keep
            if changed:
                SeqIO.write(records, str(target / path.name), "genbank")
            else:
                shutil.copyfile(path, target / path.name)
    print(json.dumps({"datasets": datasets, "removed": removed}))


if __name__ == "__main__":
    main()
