# MapMan Bin Enrichment analysis

Developed from an example provided by Marie Bolger at Forschungszentrum Juelich on 2019-14-01
That code is implemented on the [plabi website](https://plabipd.de/portal/bin-enrichment).

This package loads a Mercator4 mapping file and a DEG file containing matching target IDs, p-values and fold-change values  
in order to assess enrichment based on Fisher test (2 sided) for Up/Down/Differentially regulated targets.
It then performs Benjamini Hochberg correction for false discovery rate adjustment of p-values
 and writes all this to a tsv file.

"""
Usage:

import MapmanEnrichment
MapmanEnrichment.Handler(mapping_path, out_path).perform_test(
            file_path,
			target_field='target_id',
			sig_field='qval',
			change_field='b',
			sep='\t',
			alpha=0.05
)
"""