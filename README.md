# MapMan Bin Enrichment analysis
This tool loads a Mercator4 mapping file, a list of all genes to be considered
and either a single list of differentially expressed genes or lists describing up- and down-regulated genes. 
The tool provides a report of enrichment in differential expression 
and, if up- and down-regulated gene lists are provided, a report on 
directional trends within each bin.

Testing for enrichment in differential expression is assessed by 2-sided Fisher's exact test within each bin. 
Contingency tables for this test describe 
the ratio of differentially expressed genes to other genes within the bin
vs. a background of the same ratio for all genes outside the bin.
The directional trend is similarly assessed using a contingency table 
of the ratio of up-regulated to down-regulated genes within the bin
vs. the same ratio for all genes outside the bin.

Reported values for enrichment and trend are;
  - Benjamini Hochberg corrected p-values from each test (qval_enrichment and qval_trend) 
  - the log2 transformed odds ratio from each test (log2_enrichment and log2_trend)
  - and a boolean value for significance of the enrichment (enriched),
  - and either Up, Down or None to describe any significant trend (trend_direction).

The basic enrichment functionality is implemented on the [plabi website](https://plabipd.de/portal/bin-enrichment).

Input is normally best taken from simple files with each line containing the corresponding ID from the mapping file
though may also be taken from bash arrays. If both are used they will be combined prior to analysis. 
If a simple list of differentially expressed genes is provided, 
the up- and down-regulated gene lists may not be provided.


"""
Usage for enrichment only:

python BinEnrichment/ \
  --mapping mapping.txt \ 
  --output "out_bins.tsv" \
  --background_file "all_genes.txt" \'
  --target_file "deg.txt"

OR for directional trends:

python BinEnrichment/ \
  --mapping mapping.txt \ 
  --output "out_bins.tsv" \
  --background_file "all_genes.txt" \
  --up_file "up.txt" \
  --down_file "down.txt"

"""

