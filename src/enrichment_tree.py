import scipy.stats as stats
import numpy


class EnrichmentTree:
	def __init__(self, bin_tree):
		self.bin_tree = bin_tree
		self.root_bin_node = bin_tree.root_node
		self.root_node = self.EnrichmentNode(self.root_bin_node)
		self.up_count = None
		self.down_count = None
		self.detected_count = None
		self.node_map = {self.root_bin_node.code: self.root_node}

	def find_or_add_corresponding_node(self, bin_node):
		if bin_node.code in self.node_map:
			return self.node_map[bin_node.code]
		parent_enrichment_node = self.find_or_add_corresponding_node(bin_node.parent)
		node = self.EnrichmentNode(bin_node, parent_enrichment_node)
		self.node_map[bin_node.code] = node
		return node

	def add_up_targets(self, up):
		self.up_count = len(up)
		missing_targets = set()
		for target in up:
			bin_nodes = self.bin_tree.get_nodes_for_target(target)
			if not bin_nodes:
				missing_targets.add(target)
				self.up_count -= 1
			else:
				for bn in bin_nodes:
					enrichment_node = self.find_or_add_corresponding_node(bn)
					enrichment_node.add_up_target(target)
					enrichment_node.add_detected_target(target)
		return missing_targets

	def add_down_targets(self, down):
		self.down_count = len(down)
		missing_targets = set()
		for target in down:
			bin_nodes = self.bin_tree.get_nodes_for_target(target)
			if not bin_nodes:
				missing_targets.add(target)
				self.down_count -= 1
			else:
				for bn in bin_nodes:
					enrichment_node = self.find_or_add_corresponding_node(bn)
					enrichment_node.add_down_target(target)
					enrichment_node.add_detected_target(target)
		return missing_targets

	def add_unresponsive_targets(self, unresponsive):
		missing_targets = set()
		for target in unresponsive:
			bin_nodes = self.bin_tree.get_nodes_for_target(target)
			if not bin_nodes:
				missing_targets.add(target)
				self.down_count -= 1
			else:
				for bn in bin_nodes:
					enrichment_node = self.find_or_add_corresponding_node(bn)
					enrichment_node.add_detected_target(target)
		return missing_targets

	def calculate_enrichment_helper(self, node, results):
		if node is not self.root_node:
			up_test_result = node.calculate_fisher_test(
				'up',
				self.up_count,
				self.detected_count
			)
			down_test_result = node.calculate_fisher_test(
				'down',
				self.down_count,
				self.detected_count
			)
			diff_test_result = node.calculate_fisher_test(
				'diff',
				self.up_count + self.down_count,
				self.detected_count
			)
			results.append(
				self.EnrichmentResult(
					node.bin_node,
					len(node.up_targets),
					len(node.down_targets),
					len(node.detected_targets),
					up_test_result[0],
					up_test_result[1],
					down_test_result[0],
					down_test_result[1],
					diff_test_result[0],
					diff_test_result[1]
				)
			)
		child_bins = node.bin_node.children
		for child_bin in child_bins:
			if child_bin.code in self.node_map:
				child_node = self.node_map[child_bin.code]
				self.calculate_enrichment_helper(child_node, results)

	def calculate_enrichment(self):
		print('Calculate enrichment for all bins')
		results = []
		self.calculate_enrichment_helper(self.root_node, results)
		return results

	class EnrichmentNode:
		def __init__(self, bin_node, parent=None):
			self.bin_node = bin_node
			self.parent = parent
			self.up_targets = set()
			self.down_targets = set()
			self.detected_targets = set()

		def add_up_target(self, target):
			self.up_targets.add(target)
			if self.parent:
				self.parent.add_up_target(target)

		def add_down_target(self, target):
			self.down_targets.add(target)
			if self.parent:
				self.parent.add_down_target(target)

		def add_detected_target(self, target):
			self.detected_targets.add(target)
			if self.parent:
				self.parent.add_detected_target(target)

		def calculate_fisher_test(
				self,
				differential_type,
				relevant_root_diff_count,  # i.e. up/down/diff depending on type
				root_detected_count
		):
			if differential_type == 'up':
				diff = len(self.up_targets)
			elif differential_type == 'down':
				diff = len(self.down_targets)
			else:
				diff = len(self.up_targets) + len(self.down_targets)
			detected = len(self.detected_targets)
			not_diff = detected - diff
			other_diff = relevant_root_diff_count - diff
			other_detected = root_detected_count - detected
			other_not_diff = other_detected - other_diff
			p_val = stats.fisher_exact(
				[
					[diff, not_diff],
					[other_diff, other_not_diff]
				],
				alternative='two-sided'
			)[1]
			try:
				sample_frequency = diff / detected
			except ZeroDivisionError:
				sample_frequency = 1
			try:
				background_frequency = relevant_root_diff_count / root_detected_count
			except ZeroDivisionError:
				background_frequency = 1
			enrichment = sample_frequency / background_frequency
			if enrichment == 0:
				log2_enrichment = float('-inf')
			else:
				log2_enrichment = numpy.log2(enrichment)
			return p_val, log2_enrichment

	class EnrichmentResult:
		def __init__(
				self,
				bin_node,
				up,
				down,
				detected,
				up_fisher,
				up_log2_enrichment,
				down_fisher,
				down_log2_enrichment,
				diff_fisher,
				diff_log2_enrichment
		):
			self.bin_node = bin_node
			self.up = up
			self.down = down
			self.detected = detected
			self.up_fisher = up_fisher
			self.up_log2_enrichment = up_log2_enrichment
			self.up_bh_fisher = None
			self.down_fisher = down_fisher
			self.down_log2_enrichment = down_log2_enrichment
			self.down_bh_fisher = None
			self.diff_fisher = diff_fisher
			self.diff_log2_enrichment = diff_log2_enrichment
			self.diff_bh_fisher = None

		def to_dict(self):
			d = {
				'bin_code': self.bin_node.code,
				'bin_name': self.bin_node.name,
				'up': self.up,
				'down': self.down,
				'detected': self.detected,
				'up_qval': self.up_bh_fisher,
				'down_qval': self.down_bh_fisher,
				'diff_qval': self.diff_bh_fisher,
				'up_log2_enrichment': self.up_log2_enrichment,
				'down_log2_enrichment': self.down_log2_enrichment,
				'diff_log2_enrichment': self.diff_log2_enrichment
			}
			return d
