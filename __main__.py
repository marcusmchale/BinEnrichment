#!/usr/bin/env python

from src import Handler
import sys

mapping_path = sys.argv[1]
out_path = sys.argv[2]
data_path = sys.argv[3]

# Limma
#Handler(mapping_path, out_path).perform_test(
#	data_path,
#	target_field='GeneID',
#	lrt_sig_field='adj.P.Val',
#	change_field='logFC',
#)


#Sleuth
Handler(mapping_path, out_path).perform_test(
	data_path,
	target_field='target_id',
	lrt_sig_field='qval.LRT',
	interaction_sig_field='interaction_qval',
	main_effect_sig_field='main_effect_qval',
	wt_sig_field='qval.WT',
	change_field='b',
	sep='\t',
	alpha=0.05
)