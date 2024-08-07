# General Information

If you use this code, please cite our [work](https://doi.org/10.3389/fbinf.2022.834034):
Morlion, A.; Hulstaert, E.; Vromman, M.; Anckaert, J.; Everaert, C.; Vandesompele, J.; Mestdagh, P. CiLiQuant: Quantification of RNA Junction Reads Based on Their Circular or Linear Transcript Origin. Frontiers in Bioinformatics, 2022. doi: 10.3389/fbinf.2022.834034

Feel free to ask questions/inform about issues in the [issues section](https://github.com/OncoRNALab/CiLiQuant/issues)

# CiLiQuant.py

This script separates junction reads based on their linear or circular origin (if a forward splice junction falls within another backsplice junction, the circular or linear origin cannot be determined -> ambiguous). Only the non-ambiguous junction reads are used to compare the relative linear and circular transcript abundance. 

Together with the circ fraction (circ/(circ+lin)), a 95% confidence interval is provided.

## Installation & dependencies
The following libraries are required:
* [Python](https://www.python.org) 3.6
* [Pandas](https://pandas.pydata.org) 1.0

To use this tool, you can simply download the CiLiQuant.py script or clone the repository using the command line:
```
git clone https://github.com/OncoRNALab/CiLiQuant.git
```

## Input
Each time BED-format for first 6 columns (chromosome, start, stop, name, score, strand), but score column may be filled with dots instead of numbers or may already contain the nr of reads
- a (forward-splice) junction file that contains coordinates of junctions and number of junction spanning reads (e.g. TopHat's junctions.bed or STAR's SJ.out.tab)
  - Note STAR's SJ.out.tab is not in BED-format and strand info formatting needs to be consistent with that of the back-splice and exon or gene file (+/-). This can be changed with:
    ```awk 'BEGIN {OFS="\t"} { if ($4 == 1) { $4 = "+" } else if ($4 == 2) { $4 = "-" } print $1,$2,$3,".",".",$4,$7 }' SJ.out.tab > SJ.out.tab.modified ```
- a corresponding back-splice junction file that contains coordinates of backsplice junctions and number of junction spanning reads (e.g. CIRCexplorer or find_circ output)
- an exon or gene file that contains start and stop positions of the (exons of the) genes of interest. Note that fractions will only be calculated for circRNAs that can be assigned to these genes of interest.

Note that if you want to exclude junction reads with low counts for the circular fraction calculation, you need to remove these junctions in the input files before running the script.

## Run script
```
Python3 CiLiQuant.py --help                                                                       
usage: CiLiQuant.py [-h] -j JUNCTIONS -b BACKSPLICE -e EXONS -bc BSREADS_COLUMN -fc FSREADS_COLUMN [-v OVERLAP_COLUMN] [-o OUTPUT]
                    [-n NAME_PREFIX] [-ff FSFILTER] [-bf BSFILTER] [-s {yes,no}]

optional arguments:
  -h, --help            show this help message and exit
  -j JUNCTIONS, --junctions JUNCTIONS
                        Normal (forward) junction file (BED format)
  -b BACKSPLICE, --backsplice BACKSPLICE
                        Backsplice junction file (BED format)
  -e EXONS, --exons EXONS
                        Exon or gene file (BED format)
  -bc BSREADS_COLUMN, --bsreads_column BSREADS_COLUMN
                        Column in backsplice junction file that contains the number of junction reads
  -fc FSREADS_COLUMN, --fsreads_column FSREADS_COLUMN
                        Column in forward junction file that contains the number of junction reads
  -v OVERLAP_COLUMN, --overlap_column OVERLAP_COLUMN
                        In case the start and stop columns contain the max spanning read positions instead of exact sites, indicate which
                        column contains overlap at left and right side of junction (e.g. column 11 in TopHat junction file, not needed in
                        STAR junction file)
  -o OUTPUT, --output OUTPUT
                        Optional output directory
  -n NAME_PREFIX, --name_prefix NAME_PREFIX
                        Optional pefix (e.g. sample name) for output files
  -ff FSFILTER, --fsfilter FSFILTER
                        Filter out forward junctions with less reads than this filter (default=1)
  -bf BSFILTER, --bsfilter BSFILTER
                        Filter out backsplice junctions with less reads than this filter (default=1)
  -s {yes,no}, --strand {yes,no}
                        Run in strand-specific mode or not? (default=yes)
```
Remark: A long exon file will require more time to process.
If it really takes too long, you can split the exon file per chromosome and run the script on the exon file subsets (no need to split the forward splice and backsplice junction files - the script will only use the junctions that are between the exon file boundaries).

## Output
- **CiLiQuant_circ.txt**: a file that for every unique backsplice junctions (within boundaries of exon file), gives the nr of left and right flanking junction reads together with the estimated circRNA fraction for this site and CI	
  - *bs_reads*: total number of backsplice reads for each circ_id
  - (*lfl_reads, rfl_reads*): sum of linear only reads that partially overlap with backsplice junction on the left and right resp.
  - (*lfl_junctions, rfl_junctions*): number of unique flanking junctions on the left and right of circ_id resp. (sometimes 2 or more unique junctions partially overlap with the backsplice on one side)
  - (*lfl_junctions_ambi, rfl_junctions_ambi*): number of ambiguous (inside other backsplice) flanking junctions on left and right of circ_id resp.
  - (*linfl_reads, linfl_reads_av*): sum of linear only flanking reads and average (sum divided by nr of linear only flanking junctions) resp.
  - *circ_fraction_fl*: circRNA fraction for circ_id based on its backsplice reads and linear flanking junction reads: bs_reads/(bs_reads+linfl_reads_av)
  - (*p_AC_fl, ci_lower_AC_fl, ci_upper_AC_fl*): cf. circ_fraction_fl, here estimated proportion and CI based on modified Wald (Agresti-Coull, 1998; Brown et al, 2001)
  - *lin_reads_gene_av*: average nr of linear only junction reads in gene; see lin_reads_av in CiLiQuant_gene.txt below
  - *circ_fraction_all*: circRNA fraction for circ_id based on its backsplice reads and all linear only junction reads in gene: bs_reads/(bs_reads+lin_reads_gene_av)
  - (*p_AC_all, ci_lower_AC_all, ci_upper_AC_all*): estimated proportion and CI of circ_fraction_all based on modified Wald (Agresti-Coull, 1998; Brown et al, 2001)
- **CiLiQuant_gene.txt**: a file that for every gene in the exon file shows how much backsplice, linear only and ambiguous reads and junctions there are together with the estimated circRNA fraction for this gene and CI
  - (*linear_reads, linear_junctions*): sum of linear only reads for each gene_id & corresponding nr of unique linear junctions resp.
  - (*ambiguous_reads, ambiguous_junctions*): sum of ambiguous junction reads & corresponding nr of unique ambiguous junctions resp.
  - (*circ_reads, circ_junctions*): sum of circular (backsplice) reads for each gene_id & corresponding nr of unique backsplice junctions resp.
  - *lin_reads_av*: linear only junction reads corrected for nr of linear only junctions
  - *circ_reads_av*: circular (backsplice) only junction reads corrected for nr of circular only junctions
  - *circ_fraction*: circRNA fraction in gene based on average of all backsplice and linear only junction reads in gene: circ_reads_av/(circ_reads_av+lin_reads_av)
  - (*p_AC, ci_lower_AC, ci_upper_AC*): estimated proportion and CI of circ_fraction based on modified Wald (Agresti-Coull, 1998; Brown et al, 2001)
