from .tree import Tree, Expression
from .enrichment import TreeTester
import csv

class TestHandler:
    def __init__(self, mapping_path, out_path, background, diff=None, up=None, down=None, alpha=0.05):
        self.mapping_path = mapping_path
        self.tree = Tree()
        self.tree.load_map(self.mapping_path)

        self.out_path = out_path

        # convert lists to sets for testing
        self.background = set(background)
        self.directional = any([up, down])

        self._diff = set() if diff is None else set(diff)
        self.up = set() if up is None else set(up)
        self.down = set() if down is None else set(down)

        if self.directional and self._diff:
            raise ValueError("Diff may only be provided if no up or down is provided")

        if not self.diff:
            raise ValueError("No DEGs provided")

        self.undetermined = self.background - self.diff

        self.alpha = alpha

    @property
    def diff(self):
        if self.directional:
            return self.up | self.down
        else:
            return self._diff

    def perform_test(self):
        tree = self.tree
        tree.add_degs(self.up, Expression.UP)
        tree.add_degs(self.down, Expression.DOWN)
        tree.add_degs(self._diff, Expression.DIFF)
        tree.add_degs(self.undetermined, Expression.UNDETERMINED)
        missing_targets = tree.unmapped_genes
        if missing_targets:
            print(f'Some targets of interest were not found in the mapping file: {",".join(missing_targets)}')
        results = TreeTester(tree).calculate_enrichment()
        results = TreeTester.fdr_correction(results, self.alpha)
        print('Write results to: ' + str(self.out_path))
        self.write_file(results)
        tree.clear_degs()

    def write_file(self, results):
        with open(self.out_path, 'w') as tsv_file:
            fieldnames = [
                'bin_code',
                'bin_name',
                'up',
                'down',
                'diff',
                'detected',
                'log2_enrichment',
                'qval_enrichment',
                'enriched',
                'log2_trend',
                'qval_trend',
                'trend_direction'
            ]
            writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerow({
                    'bin_code': results[0][0].code,
                    'bin_name': results[0][0].name,
                    'up': len(results[0][0].expression_map.up),
                    'down': len(results[0][0].expression_map.down),
                    'diff': len(results[0][0].expression_map.diff),
                    'detected': len(results[0][0].expression_map.detected)
            })
            for record in results[1:]:
                node = record[0]
                node_results = record[1]
                if (node_results.diff.enrichment > 0) & (node_results.diff.qval <= self.alpha):
                    enriched = "True"
                else:
                    enriched = "False"
                if node_results.trend.qval <= self.alpha:
                    if node_results.trend.enrichment <= 0:
                        trend_direction = 'Down'
                    else:
                        trend_direction = 'Up'
                else:
                    trend_direction = "None"
                writer.writerow({
                    'bin_code': node.code,
                    'bin_name': node.name,
                    'up': len(node.expression_map.up),
                    'down': len(node.expression_map.down),
                    'diff': len(node.expression_map.diff),
                    'detected': len(node.expression_map.detected),
                    'log2_enrichment': node_results.diff.enrichment,
                    'log2_trend': node_results.trend.enrichment,
                    'enriched': enriched,
                    'qval_enrichment': node_results.diff.qval,
                    'qval_trend': node_results.trend.qval,
                    'trend_direction': trend_direction,
                })
