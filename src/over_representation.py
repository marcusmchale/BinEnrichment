from .bin_tree import BinTree
from .enrichment_tree import EnrichmentTree
import csv


class Handler:
	def __init__(self, mapping_path, out_path):
		self.mapping_path = mapping_path
		self.out_path = out_path
		self.up_targets = set()
		self.down_targets = set()
		self.unresponsive_targets = set()

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
		bin_tree = BinTree()
		bin_tree.load_mapping_file(self.mapping_path)
		en_tree = EnrichmentTree(bin_tree)
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
		en_tree.detected_count = (
				len(self.up_targets) +
				len(self.down_targets) +
				len(self.unresponsive_targets)
		)
		missing_targets = en_tree.add_up_targets(self.up_targets)
		missing_targets = missing_targets.union(en_tree.add_down_targets(self.down_targets))
		missing_targets = missing_targets.union(en_tree.add_unresponsive_targets(self.unresponsive_targets))
		if missing_targets:
			print(
				'Some targets of interest were not found in the mapping file: '
				+ ','.join(missing_targets)
			)
		# add handler for return from these (unrecognised genes)
		results = en_tree.calculate_enrichment()
		results = self.fdr_correction(results)
		print('Write results to: ' + self.out_path)
		self.write_file(
			results,
			en_tree.up_count,
			en_tree.down_count,
			en_tree.detected_count
		)
		print(
			f"Background DEG frequencies:\n"
			f"UP: {len(self.up_targets)}/{en_tree.detected_count}\n"
			f"DOWN: {len(self.down_targets)}/{en_tree.detected_count}"
		)

	@staticmethod
	def fdr_correction(results):
		# Apply Benjamini-Hochberg FDR correction of p-values
		print('Apply FDR correction')
		len_results = len(results)
		results.sort(key=lambda x: x.up_fisher)
		for i, record in enumerate(results):
			up_bh_fisher = record.up_fisher * len_results / (i + 1)
			if up_bh_fisher > 1:
				up_bh_fisher = 1
			record.up_bh_fisher = up_bh_fisher
		results.sort(key=lambda x: x.down_fisher)
		for i, record in enumerate(results):
			down_bh_fisher = record.down_fisher * len_results / (i + 1)
			if down_bh_fisher > 1:
				down_bh_fisher = 1
			record.down_bh_fisher = down_bh_fisher
		results.sort(key=lambda x: x.diff_fisher)
		for i, record in enumerate(results):
			diff_bh_fisher = record.diff_fisher * len_results / (i + 1)
			if diff_bh_fisher > 1:
				diff_bh_fisher = 1
			record.diff_bh_fisher = diff_bh_fisher
		return results

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
		up, down, unresponsive = self.read_file(
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
		self.up_targets = up
		self.down_targets = down
		self.unresponsive_targets = unresponsive

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
		up_list = list()
		down_list = list()
		unresponsive_list = list()
		for file_path in file_paths:
			up, down, unresponsive = self.read_file(
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
			up_list.append(up)
			down_list.append(down)
			unresponsive_list.append(unresponsive)
		detected = set.intersection(*up_list, *down_list, *unresponsive_list)
		up_frequency = {g: 0 for g in detected}
		down_frequency = {g: 0 for g in detected}
		for gene_list in up_list:
			for g in gene_list:
				up_frequency[g] += 1
		for gene_list in down_list:
			for g in gene_list:
				down_frequency[g] += 1
		self.up_targets = set([g for g, v in up_frequency.items() if min_prop <= v/len(file_paths)])
		self.down_targets = set([g for g, v in down_frequency.items() if min_prop <= v / len(file_paths)])
		self.unresponsive_targets = detected - self.up_targets - self.down_targets

	def read_file(
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
							'Adding this target (%s) to the background '
							'rather than those that exhibit a change' % target
						)
						unresponsive.add(target)
				else:
					unresponsive.add(target)
		return up, down, unresponsive

	def write_file(self, results, root_up, root_down, root_detected):
		with open(self.out_path, 'w') as tsv_file:
			fieldnames = [
				'bin_code',
				'bin_name',
				'up',
				'down',
				'detected',
				'diff_log2_enrichment',
				'diff_qval',
				'up_log2_enrichment',
				'up_qval',
				'down_log2_enrichment',
				'down_qval',
			]
			writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
			writer.writeheader()
			writer.writerow(
				{
					'bin_code': '0',
					'bin_name': 'root',
					'up': root_up,
					'down': root_down,
					'detected': root_detected
				}
			)
			for row in results:
				writer.writerow(row.to_dict())
