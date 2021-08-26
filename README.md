# Inference_of_Gene_regulatory_networks

The aim of this project was to infer a gene regulatory network from single cell data using Coulas to find differentially coexpressed genes.<br />
The data used a example data consists of 3k PBMCs from a Healthy Donor and is  freely available from 10x Genomics.
It can be downloaded from this [web site](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k) or using 
>wget http://<i></i>cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
>
on linux.

# Method

After doing the preprocessing on should end up with 2 clusters of expression data.
In order to infer a gene regulatory network from the single cell data a copula is computed for a pair of transcription factor and gene. The distributions of each genes expression data are used as the marginal distributions. This process is done twice, once in each cluster. Afterwards the KS-Distance between the 2 Copulas is computed.


# Documentation 
Script name | functionality
------------ | -------------
Preprocessing_singlecell.ipynb | Does the preprocessing of the single cell data as described in the [scanpy tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
Extract_clusters.ipynb | Extracts clusters from the preprocessed data
Copula_Networkinference.py | <ul><li>Computes KS-Distances between a Tf and a gene</li><li>Can generate random KS-Distances</li></ul>
Merge_bootstrapps.py | Merges multiple Bootstrapps into one file and computes the mean squared error
Filtering.py | Filters the result based on error and assigns p-values to the KS-Distances (if a file with random KS-Distances is provided)
shufflegeneexpressiondata.py | <ul><li>shuffles the entire data</li><li>shuffles the entrie data except for the transcription factors</li></ul>
PlottingFunctions.py | <ul><li>plot all KS-Distances for edges starting from a Tf</li><li>plot the gene expression data of a gene</li><li>plot the distribution of errors</li><li>plot the distribution of p-values</li><li>plot the KS-Distances for edges starting from a Tf against the same edges from a differnt file</li><li>plot the copula for a Tf/gene pair in each cluster</li></ul>
PlottingVenn.py | Compare the edges most likely to be in the network against a gold standard to validate the result
