#!/usr/bin/env python

from src import Handler
import sys

mapping_path = sys.argv[1]
out_path = sys.argv[2]
data_path = sys.argv[3]

Handler(mapping_path, out_path).perform_test(
	data_path,
	target_field='GeneID',
	lrt_sig_field='adj.P.Val',
	change_field='logFC',
)
