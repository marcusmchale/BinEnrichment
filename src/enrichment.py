from numpy import log2, errstate
from scipy.stats import fisher_exact

from .tree import Node, Tree


class FisherResult:
	def __init__(self, result):
		self.pval = result[1]
		with errstate(divide='ignore'):
			self.enrichment = log2(result[0])
		self.qval = None


class EnrichmentResults:
	def __init__(
			self,
			diff: FisherResult,
			trend: FisherResult
	):
		self.diff = diff
		self.trend = trend
		self.map = {
			'diff': self.diff,
			'trend': self.trend,
		}


class NodeTester:
	def __init__(self, node: Node, root: Node):
		self.node = node
		self.root = root

	def get_tables(self):
		detected = len(self.node.expression_map.detected)
		up = len(self.node.expression_map.up)
		down = len(self.node.expression_map.down)
		diff = len(self.node.expression_map.diff)

		other_detected = len(self.root.expression_map.detected - self.node.expression_map.detected)
		other_up = len(self.root.expression_map.up - self.node.expression_map.detected)
		other_down = len(self.root.expression_map.down - self.node.expression_map.detected)
		other_diff = len(self.root.expression_map.diff - self.node.expression_map.detected)

		diff_table = [
			[diff, detected-diff],
			[other_diff, other_detected - other_diff]
		]
		trend_table = [
			[up, down],
			[other_up, other_down]
		]
		return {
			'diff': diff_table,
			'trend': trend_table
		}

	def get_results(self):
		tables = self.get_tables()
		diff = FisherResult(fisher_exact(tables['diff'], alternative='two-sided'))
		trend = FisherResult(fisher_exact(tables['trend'], alternative='two-sided'))
		return EnrichmentResults(diff, trend)


class TreeTester:
	def __init__(self, tree: Tree):
		self.tree = tree

	def calculate_enrichment_helper(self, node, results):
		if node is self.tree.root:
			results.append((node, None))
		else:
			tester = NodeTester(node, self.tree.root)
			node_results = tester.get_results()
			results.append(
				(node, node_results)
			)
		for child in node.children:
			self.calculate_enrichment_helper(child, results)

	def calculate_enrichment(self):
		print('Calculate enrichment for all bins')
		results = []
		self.calculate_enrichment_helper(self.tree.root, results)
		return results

	@staticmethod
	def fdr_correction(results, alpha):

		def bh_helper(results, key):
			# Apply Benjamini-Hochberg FDR correction of p-values
			len_results = len(results)
			results.sort(key=lambda x: x[1].map[key].pval)
			current_pval = None
			index = 0  # not using enumerate as we want to handle identical p-values without incrementing index
			# this is the method implemented in multtest package on bioconductor.
			# at least according to stack exchange https://stats.stackexchange.com/questions/18872
			for record in results:
				pval = record[1].map[key].pval
				if pval != current_pval:
					index += 1
				qval = pval * len_results / index
				if qval > 1:
					qval = 1
				record[1].map[key].qval = qval

		print('Apply FDR correction')
		root = results[0]
		results = results[1:]
		for key in results[0][1].map.keys():
			bh_helper(results, key)
		results.sort(key=lambda x: (not (x[1].diff.qval <= alpha and x[1].diff.enrichment >0), x[1].diff.qval))
		results.insert(0, root)
		return results
