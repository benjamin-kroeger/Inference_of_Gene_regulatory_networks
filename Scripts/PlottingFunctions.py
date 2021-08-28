import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import random
import os


class plot_distributions:

    """This class is used to create plots for the Ks_distances"""

    def __init__(self, pathtoinferedscores):

        # check how toread the file with the scores
        if pathtoinferedscores[-4:] == '.csv':
            self.tf_gene_score_df = pd.read_csv(pathtoinferedscores)
            # read the file
            if len(self.tf_gene_score_df.columns) == 3:
                self.tf_gene_score_df.columns = ['TF', 'gene', 'score']
            elif len(self.tf_gene_score_df.columns) == 4:
                self.tf_gene_score_df.columns = ['TF', 'gene', 'score', 'error']

            self.filename = pathtoinferedscores.split('/')[-1].rstrip(".csv")

        if pathtoinferedscores[-3:] == '.h5':
            self.tf_gene_score_df = pd.read_hdf(pathtoinferedscores,'df')
            self.filename = os.path.basename(pathtoinferedscores).rstrip('.h5')

    def plotvaluedistribution(self, tfs):

        """This function plots all values for the edges starting from given transcription factors"""

        # dataframes with only one transrcription factor
        tf_dfs = []

        # create multiple dataframes for each transcription factor
        for tf in tfs:
            tf_df = self.tf_gene_score_df.loc[self.tf_gene_score_df['TF'] == (' ' + tf)]
            tf_dfs.append(tf_df)

        # plot the ks_distances as distrinbution and kde
        for tf_df in tf_dfs:
            values = tf_df['score'].to_list()

            ax = sns.displot(values, kde=True)
            ax.set_xlabels(tf_df['TF'].iloc[0])
            plt.tight_layout()
            plt.savefig(tf_df['TF'].iloc[0] + '_distribution.png')

            ax = sns.displot(values, kind='kde')
            ax.set_xlabels(tf_df['TF'].iloc[0])
            plt.tight_layout()
            plt.savefig(tf_df['TF'].iloc[0] + '_density.png')


    def plotdistofrandomGenes(self,amount):

        """This function plots the expression data of a random gene"""

        allgenes = self.cluster.columns.tolist()

        for i in range(amount):
            gene = random.choice(allgenes)

            values = self.cluster.filter(items=[gene])

            ax = sns.displot(values, kde=True)
            ax.set_xlabels(gene)
            plt.tight_layout()
            plt.savefig(self.filename+ '_' + gene + '_Distribution.png')




    def ploterrordistribution(self):

        """This function plots the distribution of errors"""

        if len(self.tf_gene_score_df.columns) == 3:
            print("There is no error column")
            return

        errors = self.tf_gene_score_df[self.tf_gene_score_df.columns[3]]

        ax = sns.displot(errors, kde=True)
        ax.set_xlabels("Mean absolute Error from the Mean")
        plt.tight_layout()
        plt.savefig(self.filename + '_standarderror_dist.png')

    def plotpvaluedistribution(self):

        """"This function plots the distribution of p-values"""
        if len(self.tf_gene_score_df.columns) <=4:
            print("There is no pvalue column")
            return

        pvalues = self .tf_gene_score_df[self.tf_gene_score_df.columns[4]]

        ax = sns.displot(pvalues, kde=True)
        ax.set_xlabels("pvalues")
        plt.tight_layout()
        plt.savefig(self.filename + 'pvalue_dist.png')


    def plotcorrelationksdiststderror(self):

        """This function is used to plot the correlation between the ks_distance and error"""

        if len(self.tf_gene_score_df.columns) == 3:
            print("There is no std error column")
            return

        corrmat = np.corrcoef(self.tf_gene_score_df['score'],self.tf_gene_score_df['error'])
        print(corrmat)
        print(pearsonr(self.tf_gene_score_df['score'],self.tf_gene_score_df['error']))

        plt.scatter(self.tf_gene_score_df['score'],self.tf_gene_score_df['error'])
        plt.xlabel('Ksdist')
        plt.ylabel('stderror')
        plt.show()

class plot_distributions_against:

    """This class contains functions to plot 2 distributions against each other"""

    def __init__(self, path0, path1):

        # path0 and path1 should lead to two different files with computed ks_distances


        if path0.endswith('.csv'):
            self.tf_gene_score_df0 = pd.read_csv(path0, header=None)
            self.assign_column_names(self.tf_gene_score_df0)
        elif path0.endswith('.h5'):
            self.genecluster0 = pd.read_hdf(path0, 'df')

        if path1.endswith('.csv'):
            self.tf_gene_score_df1 = pd.read_csv(path1, header=None)
            self.assign_column_names(self.tf_gene_score_df1)
        elif path1.endswith('.h5'):
            self.genecluster1 = pd.read_hdf(path1, 'df')


    def assign_column_names(self, df):

        if len(df.columns) == 3:
            df.columns = ['TF', 'gene', 'score']
        elif len(df.columns) == 4:
            df.columns = ['TF', 'gene', 'score', 'error']


    def computerandomCopula(self):

        allgenes = self.genecluster0.columns

        gene0=random.choice(allgenes)
        gene1=random.choice(allgenes)

        df_cluster0 = pd.concat([self.genecluster0[gene0], self.genecluster0[gene1]], axis=1)
        df_cluster1 = pd.concat([self.genecluster1[gene0], self.genecluster1[gene1]], axis=1)

        self.computecopulaplots(df_cluster0, df_cluster1)


    def computespecificCopula(self,gene0,gene1):

        df_cluster0 = pd.concat([self.genecluster0[gene0],self.genecluster0[gene1]], axis=1)
        df_cluster1 = pd.concat([self.genecluster1[gene0], self.genecluster1[gene1]], axis=1)

        self.computecopulaplots(df_cluster0,df_cluster1)


    def computecopulaplots(self, cluster0, cluster1):

        """This function creates 2 plots with the 2 different marginal distributions.
        The first plot is for the copula for the 2 genes in cluster0.
        The first plot is for the copula for the 2 genes in cluster1.
        """

        dist0 = np.array(cluster0.iloc[:,0].tolist())
        dist1 = np.array(cluster0.iloc[:,1].tolist())

        plot = sns.jointplot(x=dist0, y=dist1, kind='kde')
        plot.set_axis_labels(cluster0.columns[0], cluster0.columns[1], fontsize=16)
        plt.tight_layout()
        plt.title('cluster0')
        plt.show()

        dist0 = np.array(cluster1.iloc[:, 0].tolist())
        dist1 = np.array(cluster1.iloc[:, 1].tolist())

        plot = sns.jointplot(x=dist0, y=dist1, kind='kde')
        plot.set_axis_labels(cluster1.columns[0], cluster1.columns[1], fontsize=16)
        plt.tight_layout()
        plt.title('cluster1')
        plt.show()

    def plot_dists_against(self, tfs):

        """This function plots the ks_distances from cluster 0 against the ks_distances in cluster 1
        for a specific transcription factor"""


        tf_dfs0 = []
        tf_dfs1 = []

        for tf in tfs:
            tf_df = self.tf_gene_score_df0.loc[self.tf_gene_score_df0['TF'] == (' ' + tf)]
            tf_dfs0.append(tf_df)
            tf_df = self.tf_gene_score_df1.loc[self.tf_gene_score_df1['TF'] == (' ' + tf)]
            tf_dfs1.append(tf_df)

        for df0, df1 in zip(tf_dfs0, tf_dfs1):
            values0 = df0['score'].to_list()
            values1 = df1['score'].to_list()

            # zip the 2 scores together to plot them
            score_df = pd.DataFrame(zip(values0, values1))
            score_df.columns = ['without', 'with bootstrap']

            ax=sns.displot(data=score_df, kind= 'kde')
            ax.set_xlabels(df0['TF'].iloc[0])
            plt.tight_layout()
            plt.savefig(df0['TF'].iloc[0] + '_density_comparisson.png')

            ax =sns.displot(data=score_df)
            ax.set_xlabels(df0['TF'].iloc[0])
            plt.tight_layout()
            plt.savefig(df0['TF'].iloc[0] + '_distribution_comparisson.png')

def plotvalues(pathtovalues):

    """This function plots the randomly generated Ks_distances"""
    with open(pathtovalues,'rb') as f:
        randks  =np.load(f)

        ax = sns.displot(randks, kde=True)
        ax.set_xlabels("Random KS-dists")
        plt.tight_layout()
        plt.savefig('RandomKsDistances_10000tfnormalclusters_varied.png')


def plotGenie3scores(pathtoGenie3):

    Genie3_df =  pd.read_csv(pathtoGenie3,names=['Tf', 'gene', 'score'])
    plt.figure(figsize=(9.0, 4.8))
    scores = Genie3_df['score'].to_list()
    plt.title('Distribution of Genie3 Scores (roughly 19_000 are bigger than 1')

    scoredist = sns.histplot(scores,kde=True)

    scoredist.set(xlim=(0, 5))
    plt.show()


plotGenie3scores( '../Gene_correlation_data/GENIE3_output_all_tfs.csv')


# test = plot_distributions('../Gene_correlation_data/CopulaKSdist_beta_average_with70percent_randomsampling_MAE_pvalues.h5')
# test.plotpvaluedistribution()
