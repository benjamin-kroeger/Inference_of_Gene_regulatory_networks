import numpy as np
import pandas as pd
import random
import ntpath
import concurrent.futures


class shuffleGeneExpressions:
    """This class contains function to shuffle the gene expression data of a data frame, to further randomize the approach to get
    random ks_distances."""

    def __init__(self, pathtoadata):

        self.filename = ntpath.basename(pathtoadata)

        self.adata_df = pd.read_hdf(pathtoadata, 'df')

    def shuffleverything(self):

        # extract shape
        columnnames = self.adata_df.columns.tolist()
        numberrows, _ = self.adata_df.shape
        _, numberrcols = self.adata_df.shape

        # put all values into a list
        allvalues = []
        for row in self.adata_df.values.tolist():
            allvalues.extend(row)

        # shuffle the list
        random.shuffle(allvalues)

        # rebuild the structure of the df
        newlist = []
        counter = 0
        for i in range(0, numberrows):
            newlist.append(allvalues[counter: (counter + numberrcols)])
            counter = counter + numberrcols

        # recreate the dataframe
        newlist = np.array(newlist)

        self.adata_df = pd.DataFrame(newlist)
        self.adata_df.columns = columnnames

        print(self.adata_df)

    def shufflewithouttfs(self):

        genes = self.adata_df.columns.tolist()
        # get all tfs
        with open('../Human_allTFs.txt', 'r') as f:
            tfs = [gene.strip() for gene in f.readlines()]
            tfs = [tf for tf in tfs if tf in genes]

        # safe protion of df with tfs in it
        tf_df = self.adata_df.loc[:, tfs]
        tf_df.index = list(range(len(tf_df)))
        print(tf_df.index)
        tf_df.reset_index(drop=True)

        # drop the tfs
        self.adata_df.drop(tfs, axis=1, inplace=True)

        self.shuffleverything()

        self.adata_df = pd.concat([self.adata_df, tf_df], axis=1)

        print(self.adata_df)

    def savedf(self):

        self.adata_df.to_hdf(self.filename.rstrip('.h5') + 'tfnormal_random.h5', 'df')


test = shuffleGeneExpressions('/home/benjaminkroeger/PycharmProjects/Gene_regulatory_networks/Data_h5_Scellnetor/Scellnetor_cluster1.h5')
test.shufflewithouttfs()
test.savedf()
