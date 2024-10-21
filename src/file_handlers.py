from .tree import Tree, Expression
from .enrichment import TreeTester
import csv

class GeneFileHandler:
	def __init__(self, mapping_path, out_path, background_file, up_file = None, down_file = None, diff_file = None):
		self.mapping_path = mapping_path
		self.out_path = out_path
		self.directional = any([up_file, down_file])

		self.background_file_path = background_file
		self.up_file_path = up_file
		self.down_file_path = down_file
		self.diff_file_path = diff_file

		if diff_file and self.directional:
			raise ValueError("Diff file may only be provided if no up or down file is provided")

		self.up = None
		self.down = None
		self._diff = None
		self.undetermined = None


	@property
	def diff(self):
		if self.directional:
			return self.up | self.down
		else:
			return self._diff

	@staticmethod
	def read_file(file_path):
		if file_path is None:
			return set()

		with open(file_path, 'r') as genes_file:
			return set([line.rstrip() for line in genes_file])

	def perform_test(self):
		self.up = self.read_file(self.up_file_path)
		self.down = self.read_file(self.down_file_path)
		self._diff = self.read_file(self.diff_file_path)
		self.undetermined = self.read_file(self.background_file_path) - self.diff
		tree = Tree()
		tree.load_map(self.mapping_path)
		tree.add_degs(self.up, Expression.UP)
		tree.add_degs(self.down, Expression.DOWN)
		tree.add_degs(self._diff, Expression.DIFF)
		tree.add_degs(self.undetermined, Expression.UNDETERMINED)
		missing_targets = tree.unmapped_genes
		if missing_targets:
			print(f'Some targets of interest were not found in the mapping file: {",".join(missing_targets)}')
		results = TreeTester(tree).calculate_enrichment()
		results = TreeTester.fdr_correction(results)
		print('Write results to: ' + self.out_path)
		self.write_file(results)

	def write_file(self, results):
		with open(self.out_path, 'w') as tsv_file:
			fieldnames = [
				'bin_code',
				'bin_name',
				'up',
				'down',
				'diff',
				'detected',
				'diff_log2_enrichment',
				'diff_qval',
				'enriched',
				'bias_log2_enrichment',
				'bias_qval',
				'bias_direction'
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
				if (node_results.diff.enrichment > 0) & (node_results.diff.qval <= 0.05):
					enriched = "True"
				else:
					enriched = "False"
				if node_results.bias.qval <= 0.05:
					if node_results.bias.enrichment <= 0:
						bias_direction = 'Down'
					else:
						bias_direction = 'Up'
				else:
					bias_direction = "None"
				writer.writerow({
					'bin_code': node.code,
					'bin_name': node.name,
					'up': len(node.expression_map.up),
					'down': len(node.expression_map.down),
					'diff': len(node.expression_map.diff),
					'detected': len(node.expression_map.detected),
					'diff_log2_enrichment': node_results.diff.enrichment,
					'bias_log2_enrichment': node_results.bias.enrichment,
					'enriched': enriched,
					'diff_qval': node_results.diff.qval,
					'bias_qval': node_results.bias.qval,
					'bias_direction': bias_direction,
				})


class ResultFileHandler:
	def __init__(self, mapping_path, out_path):
		self.mapping_path = mapping_path
		self.out_path = out_path
		self.up = set()
		self.down = set()
		self.undetermined = set()

	def perform_test(
			self,
			file_paths,
			target_field='target_id',
			lrt_sig_field='qval.LRT',
			interaction_sig_field='interaction_qval',
			main_effect_sig_field='main_effect_qval',
			wt_sig_field='qval.WT',
			change_field='b',
			sep='\t',
			alpha=0.05,
			min_prop=1
	):
		tree = Tree()
		tree.load_map(self.mapping_path)
		if len(file_paths) == 1:
			self.load_data(
				file_paths[0],
				target_field,
				lrt_sig_field,
				wt_sig_field,
				interaction_sig_field,
				main_effect_sig_field,
				change_field,
				sep,
				alpha
			)
		else:
			try:
				min_prop = float(min_prop)
			except ValueError as e:
				raise Exception(
					'Last argument should be minimum proportion of samples for inclusion in overlap of DEGs'
				) from e
			self.load_multiple_data(
				file_paths,
				target_field,
				lrt_sig_field,
				wt_sig_field,
				interaction_sig_field,
				main_effect_sig_field,
				change_field,
				sep,
				alpha,
				min_prop
			)
		tree.add_degs(self.up, Expression.UP)
		tree.add_degs(self.down, Expression.DOWN)
		tree.add_degs(self.undetermined, Expression.UNDETERMINED)
		missing_targets = tree.unmapped_genes
		if missing_targets:
			print(
				'Some targets of interest were not found in the mapping file: '
				+ ','.join(missing_targets)
			)
		results = TreeTester(tree).calculate_enrichment()
		results = TreeTester.fdr_correction(results)
		print('Write results to: ' + self.out_path)
		self.write_file(results)


	@staticmethod
	def is_deg(alpha, lrt_p_val=None, wt_p_val=None, interaction_p_val=None, main_effect_p_val=None):
		if interaction_p_val is not None:
			return any([
				all([i <= alpha for i in [
					p for p in [interaction_p_val, main_effect_p_val, wt_p_val] if p is not None
				]]),
				all([interaction_p_val > alpha] + [i <= alpha for i in [
					p for p in [lrt_p_val, wt_p_val] if p is not None
				]])
			])
		else:
			return all([i <= alpha for i in [p for p in [lrt_p_val, wt_p_val] if p is not None]])

	def load_data(
			self,
			file_path,
			target_field,
			lrt_sig_field,
			wt_sig_field,
			interaction_sig_field,
			main_effect_sig_field,
			change_field,
			sep,
			alpha
	):
		up, down, unresponsive = self.read_results_file(
			file_path,
			target_field,
			lrt_sig_field,
			wt_sig_field,
			interaction_sig_field,
			main_effect_sig_field,
			change_field,
			sep,
			alpha
		)
		self.up = up
		self.down = down
		self.undetermined = unresponsive

	def load_multiple_data(
			self,
			file_paths,
			target_field,
			lrt_sig_field,
			wt_sig_field,
			interaction_sig_field,
			main_effect_sig_field,
			change_field,
			sep,
			alpha,
			min_prop
	):
		#
		to_add = [i for i in file_paths if not i.startswith("-")]
		to_remove = [i[1:] for i in file_paths if i.startswith("-")]
		print(
			f"Loading files: {to_add}.\n"
			f"Targets with significant change in the same direction in a minimum of {min_prop * len(to_add)} files"
			f" will be considered DEGs.\n"
		)
		if to_remove:
			print(
				f"The following input filenames are prepended with '-': {to_remove}. "
				f"Targets from these files will be removed from the DEG lists and considered as background. "
				f"This is useful to consider unique groupings as depicted in venn diagrams. "
			)
		up_list = list()
		down_list = list()
		unresponsive_list = list()
		to_remove_up = set()
		to_remove_down = set()
		for file_path in to_add:
			up, down, unresponsive = self.read_results_file(
				file_path,
				target_field,
				lrt_sig_field,
				wt_sig_field,
				interaction_sig_field,
				main_effect_sig_field,
				change_field,
				sep,
				alpha
			)
			print(f"Including targets identified in {file_path}")
			up_list.append(up)
			down_list.append(down)
			unresponsive_list.append(unresponsive)
		for file_path in to_remove:
			up, down, unresponsive = self.read_results_file(
				file_path,
				target_field,
				lrt_sig_field,
				wt_sig_field,
				interaction_sig_field,
				main_effect_sig_field,
				change_field,
				sep,
				alpha
			)
			print(f"Marking targets from {file_path} as unresponsive")
			to_remove_up |= up
			to_remove_down |= down
		detected = set.union(*up_list, *down_list, *unresponsive_list)
		up_frequency = {g: 0 for g in detected}
		down_frequency = {g: 0 for g in detected}
		for gene_list in up_list:
			for g in gene_list:
				up_frequency[g] += 1
		for gene_list in down_list:
			for g in gene_list:
				down_frequency[g] += 1
		self.up = set([g for g, v in up_frequency.items() if g not in to_remove_up and min_prop <= v / len(to_add)])
		self.down = set([g for g, v in down_frequency.items() if g not in to_remove_down and min_prop <= v / len(to_add)])
		self.undetermined = detected - self.up - self.down

	def read_results_file(
			self,
			file_path,
			target_field,
			lrt_sig_field,
			wt_sig_field,
			interaction_sig_field,
			main_effect_sig_field,
			change_field,
			sep,
			alpha
	):
		print('Start loading data from: ' + file_path)
		up = set()
		down = set()
		unresponsive = set()
		with open(file_path) as csv_file:
			reader = csv.DictReader(csv_file, delimiter=sep, quoting=csv.QUOTE_MINIMAL)
			for row in reader:
				target = row[target_field].replace('"', '')
				lrt_p_val = float(row[lrt_sig_field].replace('"', ''))
				if wt_sig_field in reader.fieldnames:
					wt_p_val = float(row[wt_sig_field].replace('"', ''))
				else:
					wt_p_val = None
				change = float(row[change_field].replace('"', ''))
				if interaction_sig_field in reader.fieldnames:
					interaction_p_val = float(row[interaction_sig_field].replace('"', ''))
					main_effect_p_val = float(row[main_effect_sig_field].replace('"', ''))
				else:
					interaction_p_val = None
					main_effect_p_val = None
				if self.is_deg(alpha, lrt_p_val, wt_p_val, interaction_p_val, main_effect_p_val):
					if change > 0:
						up.add(target)
					elif change < 0:
						down.add(target)
					else:
						print(
							'Somehow a significant fold-change is equal to 0'
							'Adding this gene (%s) to the background '
							'rather than those that exhibit a change' % target
						)
						unresponsive.add(target)
				else:
					unresponsive.add(target)
		return up, down, unresponsive

	def write_file(self, results):
		with open(self.out_path, 'w') as tsv_file:
			fieldnames = [
				'bin_code',
				'bin_name',
				'up',
				'down',
				'detected',
				'diff_log2_enrichment',
				'diff_qval',
				'enriched',
				#'up_log2_enrichment',
				#'up_qval',
				#'down_log2_enrichment',
				#'down_qval',
				'bias_log2_enrichment',
				'bias_qval',
				'bias_direction'
				#'diff_peers_log2_enrichment',
				#'diff_peers_qval',
				#'enriched_over_peers'
			]
			writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
			writer.writeheader()
			writer.writerow({
					'bin_code': results[0][0].code,
					'bin_name': results[0][0].name,
					'up': len(results[0][0].expression_map.up),
					'down': len(results[0][0].expression_map.down),
					'detected': len(results[0][0].expression_map.detected)
			})
			for record in results[1:]:
				node = record[0]
				node_results = record[1]
				if (node_results.diff.enrichment > 0) & (node_results.diff.qval <= 0.05):
					enriched = "True"
				else:
					enriched = "False"
				if node_results.bias.qval <= 0.05:
					if node_results.bias.enrichment <= 0:
						bias_direction = 'Down'
					else:
						bias_direction = 'Up'
				else:
					bias_direction = "None"
				#if node_results.diff_peers.qval <= 0.05:
				#	if node_results.diff_peers.enrichment >= 0:
				#		enriched_over_peers = "Enriched"
				#	else:
				#		enriched_over_peers = "Depleted"
				#else:
				#	enriched_over_peers = "None"
				writer.writerow({
					'bin_code': node.code,
					'bin_name': node.name,
					'up': len(node.expression_map.up),
					'down': len(node.expression_map.down),
					'detected': len(node.expression_map.detected),
					'diff_log2_enrichment': node_results.diff.enrichment,
					# 'up_log2_enrichment': node_results.up.enrichment,
					# 'down_log2_enrichment': node_results.down.enrichment,
					'bias_log2_enrichment': node_results.bias.enrichment,
					'enriched': enriched,
					# 'diff_peers_log2_enrichment': node_results.diff_peers.enrichment,
					'diff_qval': node_results.diff.qval,
					# 'up_qval': node_results.up.qval,
					# 'down_qval': node_results.down.qval,
					'bias_qval': node_results.bias.qval,
					# 'diff_peers_qval': node_results.diff_peers.qval,
					'bias_direction': bias_direction,
					# 'enriched_over_peers': enriched_over_peers
				})


