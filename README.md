# MapMan Bin Enrichment analysis

Developed from an example provided by Marie Bolger at Forschungszentrum Juelich on 2019-14-01
That code is implemented on the [plabi website](https://plabipd.de/portal/bin-enrichment).

There may be an issue in the way that code handled gene lists that resulted in inflated gene counts in a parent bin
where a gene is present in multiple child bins. I fixed this by collecting genes into a set instead of a list. 

I expanded the functionality by separating the up/down regulated targets for independent testing
and reporting the raw counts, fold-changes and significant directional biases. 
These are all useful in interpreting and presenting the data.

This tool loads a Mercator4 mapping file and a DEG file containing matching target IDs, p-values and fold-change values  
in order to assess enrichment based on Fisher test (2 sided) for Up/Down/Differentially regulated targets.
It then performs Benjamini Hochberg correction for false discovery rate adjustment of p-values
 and writes all this to a tsv file.

"""
Usage:

Single file: 
python BinEnrichment/ sleuth mapping.txt "out_bins.tsv" "in.tsv"

Recursive on directories:

for i in `find . -name *result.txt`; 
	do python BinEnrichment/ sleuth mapping.txt "${i%_result.txt}_bins.tsv" "${i}";
done

Note for this if filenames have spaces then change the IFS
IFS=$(echo -en "\n\b")

"""