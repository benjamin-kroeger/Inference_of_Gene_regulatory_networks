# Example Data explained

File name | contents
------------ | -------------
10000randomKsdists_onrandclusters.npy | 10000 random KS-Distances computed on completly shuffled expression data
10000randomKsdists_tfnormalrandclusters.npy | 10000 random KS-Distances computed on completly shuffled expression data, except for the Tfs
10000tfnormal_varygene_randomKsdists.npy | 10000 random KS-Distances computed between the Copulas for TF/gene1 and for TF/gene2
CopulaKSdist_beta_average_with70percent....| KS-Distances from merged Bootstrapps with an Error column and a p-value colum (Dataframe)
GENIE3_output_all_tfs.csv | Edges and scores computed by genie3
cluster0.h5 | cluster0 as a Dataframe (see. Extract_clusters.py)
cluster0tfnormal_random.h5  | shuffled expression data of cluster0 except for the Transcription factors 
