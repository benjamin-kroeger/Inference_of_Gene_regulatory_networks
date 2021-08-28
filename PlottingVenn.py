import matplotlib_venn as vplot
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import random


class vennplots:

    """This class contains functions to validate the results.
    In this Code the Genie3 data is used as the gold standard against which the data is compared"""

    def __init__(self, pathtoCopulainf, pathtoGenie3):
        self.copula_df = pd.read_hdf(pathtoCopulainf, 'df')

        self.genie3_df = pd.read_csv(pathtoGenie3, names=['Tf', 'gene', 'score'])

    def plotvenndiagramnaive(self, Copulapvaluecutoff, Genie3cutoff):

        # filter the Copula df based on p-value of the ks_distance

        filteredcopula_df = self.copula_df[self.copula_df['pvalue'] >= Copulapvaluecutoff]

        # filter the Genie3 df

        filteredGenie3_df = self.genie3_df[self.genie3_df['score'] >= Genie3cutoff]

        # only keep the genenames

        filteredGenie3_df = filteredGenie3_df[['Tf', 'gene']]
        filteredcopula_df = filteredcopula_df[['Tf', 'gene']]

        # get set of edges form both datasets

        copula_tuples = [tuple(x) for x in filteredcopula_df.to_numpy()]

        genie3_tuples = [tuple(x) for x in filteredGenie3_df.to_numpy()]

        total = len(set(copula_tuples).union(set(genie3_tuples)))

        # create a venn diagramm with both sets to check if there is a significant overlap

        test = vplot.venn2([set(copula_tuples), set(genie3_tuples)], set_labels=('Copula_edges', 'Genie3_edges'),
                           subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x / total):1.0%}" + ")")
        plt.title('Copula pvalue cutoff: {0:4.1f} // Genie3 cutoff: {1:4.1f}'.format(Copulapvaluecutoff, Genie3cutoff))
        plt.savefig('Venndiagramms/Venndiagramm_pvaluecut{0:4.1f}_Genie3cut{1:4.1f}.png'.format(Copulapvaluecutoff, Genie3cutoff))
        plt.close()

    def plotvenndiagramms(self, numTopTf,randomcopulaedges):

        """This function also creats venndiagramms but has different filtering.
        Only the edges with the lowest pvlaues are kept, if randomcopulaedges is disabled.
        If it' enabled no filtering takes place."""

        # sort dfs by tf and only take the top n edges

        self.copula_df.sort_values(by=['Tf', 'pvalue'], inplace=True,ascending=True)
        self.genie3_df.sort_values(by=['Tf', 'score'], inplace=True,ascending=False)

        # get length of dfs

        orderedcopula_df = self.copula_df[['Tf', 'gene']]
        orderedgenie3_df = self.genie3_df[['Tf', 'gene']]

        orderedcopula_df.reset_index(drop=True, inplace=True)
        orderedgenie3_df.reset_index(drop=True, inplace=True)

        # get num of tf

        if self.copula_df['Tf'].nunique() == self.genie3_df['Tf'].nunique():
            numtfs = self.genie3_df['Tf'].nunique()
        else:
            raise Exception('The 2 dfs have different tfs')

        tfs = orderedcopula_df['Tf'].unique()
        edgescopula = []
        edgesgenie3 = []

        # create seperate df for each tf and take the head and put it into tuples
        for tf in tfs:

            tf_df_copula = orderedcopula_df[orderedcopula_df['Tf'] == tf]
            tf_df_copula.reset_index(drop=True,inplace=True)

            # either pick a certain amount of random genes or only keep a certain number of genes from
            # the top
            if randomcopulaedges:
                randlist = list(np.random.permutation(np.arange(0,len(tf_df_copula.index)))[:numTopTf])
                edgescopula.extend([tuple(x) for x in tf_df_copula.iloc[randlist].to_numpy()])
            else:

                edgescopula.extend([tuple(x) for x in tf_df_copula.head(numTopTf).to_numpy()])

            tf_df_genie3 = orderedgenie3_df[orderedgenie3_df['Tf'] == tf]

            edgesgenie3.extend([tuple(x) for x in tf_df_genie3.head(numTopTf).to_numpy()])

        total = len(set(edgescopula).union(set(edgesgenie3)))

        venn_plot = vplot.venn2([set(edgescopula),set(edgesgenie3)],set_labels=('Copula_edges', 'Genie3_edges'),
                                subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x / total):1.0%}" + ")")
        plt.title('Using the top {0} edeges of each Tf // randomly drawn copuaedges {1}'.format(numTopTf,randomcopulaedges))
        plt.savefig('Venndiagramms3/Venndiagramm_topgenes{0}_random_{1}.png'.format(numTopTf,randomcopulaedges))
        plt.close()


    def iterateovercutoffs(self, pstepsize):

        """This function is used to create multiple venn diagramms and iterate over
         multiple cutoffs to dertermine the best one."""

        for i in np.arange(0, 1, pstepsize):
            self.plotvenndiagramnaive(i, 0)

    def iterateovertopnum(self,stepsize,topnum):
        # starts at 10
        """This function tests different filtering parameters"""
        for i in range(1,topnum,stepsize):
            self.plotvenndiagramms(i,False)
            self.plotvenndiagramms(i,True)


plotting = vennplots('../Gene_correlation_data/CopulaKSdist_beta_average_with70percent_randomsampling_MAE_pvalues.h5',
                     '../Gene_correlation_data/GENIE3_output_all_tfs.csv')
plotting.iterateovertopnum(50,600)
