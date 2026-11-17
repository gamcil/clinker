import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('summary', Path(__file__).parents[1] / 'scripts' / 'summarize.py')
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)


class FigureComparisonTests(unittest.TestCase):
    def test_duplicate_labels_use_locus_tags(self):
        def figure(prefix):
            genes = [
                {'uid': prefix + '1', 'label': 'WP_057012404.1', 'names': {'locus_tag': 'TU73_RS12210'},
                 'start': 9208, 'end': 9679, 'strand': 1},
                {'uid': prefix + '2', 'label': 'WP_057012404.1', 'names': {'locus_tag': 'TU73_RS12240'},
                 'start': 13269, 'end': 13740, 'strand': 1},
            ]
            return {'clusters': [{'name': 'cluster', 'loci': [{'name': 'locus', 'genes': genes}]}],
                    'links': [{'query': {'uid': prefix + '1'}, 'target': {'uid': prefix + '2'}}],
                    'groups': [{'genes': [prefix + '1', prefix + '2']}]}

        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / 'python.json', Path(directory) / 'rust.json']
            for path, prefix in zip(paths, ('python', 'rust')):
                path.write_text(json.dumps(figure(prefix)))
            python, rust = (summary.figure_members(path) for path in paths)
            self.assertEqual(python, rust)
            self.assertEqual(tuple(map(len, rust)), (1, 1, 2))


if __name__ == '__main__':
    unittest.main()
