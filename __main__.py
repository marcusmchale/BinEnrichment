#!/usr/bin/env python


from src import FileHandler
import sys

result_format = sys.argv[1]
mapping_path = sys.argv[2]
out_path = sys.argv[3]
data_paths = [sys.argv[4]] if len(sys.argv) == 5 else sys.argv[4:-1]
min_prop = sys.argv[-1]
if result_format == 'sleuth':
	FileHandler(mapping_path, out_path).perform_test(
		data_paths,
		target_field='target_id',
		lrt_sig_field='qval.LRT',
		interaction_sig_field='interaction_qval',
		main_effect_sig_field='main_effect_qval',
		wt_sig_field='qval.WT',
		change_field='b',
		sep='\t',
		alpha=0.05,
		min_prop=min_prop
	)
elif result_format == 'limma':
	FileHandler(mapping_path, out_path).perform_test(
		data_paths,
		target_field='GeneID',
		lrt_sig_field='adj.P.Val',
		change_field='logFC',
		min_prop=min_prop
	)
elif result_format == 'wgcna':
	print('ok')
	FileHandler(mapping_path, out_path).perform_test(
		data_paths,
		target_field="gene",
		lrt_sig_field=sys.argv[5],
		change_field="sign",
		min_prop=1
	)
elif result_format == 'deseq2':
	print('ok')
	FileHandler(mapping_path, out_path).perform_test(
		data_paths,
		target_field="GeneID",
		lrt_sig_field="padj",
		change_field="log2FoldChange",
		min_prop=1
	)


