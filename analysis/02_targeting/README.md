Stage 02: Targeted Gene Selection
Overview

To prevent the "Curse of Dimensionality" in Machine Learning, we restrict our analysis to biologically relevant regions.
What we have (Input)

    genes_uc.txt: A list of 489 genes associated with Ulcerative Colitis (Open Targets).

What we do (Process)

    Coordinate Mapping: Convert gene symbols to genomic coordinates (GRCh37/hg19).

    BED Generation: Create a Browser Extensible Data (BED) file defining our regions of interest.

What we get (Output)

    uc_genes.bed: A file defining the "search space" for our analysis.

Theoretical Background

Machine Learning models struggle with "p > n" problems (more features than samples). By using a Hypothesis-Driven Approach, we focus on genes with known clinical precedence or GWAS associations, significantly increasing the Signal-to-Noise ratio of our future model.