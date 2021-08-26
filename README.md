# Inference_of_Gene_regulatory_networks

The aim of this project was to infer a gene regulatory network from single cell data using Coulas.<br />
The data used a example data consists of 3k PBMCs from a Healthy Donor and is  freely available from 10x Genomics.
It can be downloaded from this [web site](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k) or using 
>wget http://<i></i>cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
>
on linux.

# Documentation 
Script name | functionality
------------ | -------------
Preprocessing_singlecell.ipynb | Does the preprocessing of the single cell data as described in the [scanpy tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
Extract_clusters.ipynb | Extracts clusters from the preprocessed data
Copula_Networkinference.py |
