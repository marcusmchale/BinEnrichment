#!/usr/bin/env python


from src import Handler
import sys

mapping_path = sys.argv[2]
out_path = sys.argv[3]
data_paths = sys.argv[4:]

if sys.argv[1] == 'sleuth':
	Handler(mapping_path, out_path).perform_test(
		data_paths,
		target_field='target_id',
		lrt_sig_field='qval.LRT',
		interaction_sig_field='interaction_qval',
		main_effect_sig_field='main_effect_qval',
		wt_sig_field='qval.WT',
		change_field='b',
		sep='\t',
		alpha=0.05
	)
elif sys.argv[1] == 'limma':
	Handler(mapping_path, out_path).perform_test(
		data_paths,
		target_field='GeneID',
		lrt_sig_field='adj.P.Val',
		change_field='logFC',
	)
elif sys.argv[1] == 'custom':
	Handler(mapping_path, out_path).perform_test(
		data_paths,
		target_field=sys.argv[5],
		lrt_sig_field=sys.argv[6],
		change_field=sys.argv[7],
	)

