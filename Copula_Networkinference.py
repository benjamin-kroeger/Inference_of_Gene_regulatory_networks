import pandas as pd
import concurrent.futures
from scipy import stats
import numpy as np
import ntpath

from copulae import EmpiricalCopula
import csv
import os
import random
from tqdm import tqdm


class Edge():
    """This class is used to store the computed edges between two genes"""

    def __init__(self, start, end, ks_dist):
        # storing directed edges
        self.start = start
        self.end = end
        # this variable is stored to store the ks Distance between 2 Copulas
        self.ks_dist = ks_dist

    def __str__(self):
        return str(str(self.start) + ' ==> ' + self.end + ':\t\t' + str(self.ks_dist))


class Network_Inference():
    """This class infers edges from the single cell gene expression data.

       The class needs 2 hdf files with Gene expression data form 2 different single cell clusters.
       The transcription factors have to be given to the class.
       Randomsampling can be enabled, if multiple runs are performed to obtain bootstrapped data."""

    def __init__(self, path_to_cluster0, path_to_cluster1, transcription_factors, randomsampling):

        self.random = randomsampling
        samplesizepercentagne = 0.7

        # read the cluster data
        self.cluster0 = pd.read_hdf(path_to_cluster0, 'df')

        self.cluster1 = pd.read_hdf(path_to_cluster1, 'df')

        # if randomsampling is enabled draw 70 of each cluster as a sample
        if randomsampling:
            self.cluster0 = self.cluster0.sample(int(len(self.cluster0.index) * samplesizepercentagne))
            self.cluster1 = self.cluster1.sample(int(len(self.cluster1.index) * samplesizepercentagne))

        # check if the clusters have the same columns/genes
        if set(self.cluster0.columns.tolist()) == set(self.cluster1.columns.tolist()):
            self.genes = list(self.cluster0.columns)
        else:
            raise ValueError('The genes of the 2 clusters are not the same')

        # initialize transcriptionfactors if the df contains it
        self.transcription_factors = [factor for factor in transcription_factors if factor in self.genes]
        self.edges = []

    def compute_differential_coexpression(self, df_cluster0, df_cluster1):

        """This function computes the ks_distance between 2 copulas.
           The first copula is infered from a transcription factor gene pair from the first cluster
           The second copula is infered from the same pair but with the data from the second cluster"""

        # create copula for each pair of genes

        emp_cop_cluster0 = EmpiricalCopula(2, data=df_cluster0, smoothing='beta')
        emp_cop_cluster1 = EmpiricalCopula(2, data=df_cluster1, smoothing='beta')

        # get the cdf from each copula

        cdf_cop_0 = emp_cop_cluster0.cdf(emp_cop_cluster0.data)
        cdf_cop_1 = emp_cop_cluster1.cdf(emp_cop_cluster1.data)

        # perform the ks test distance and discard p-value

        differential_expression, _ = stats.kstest(cdf_cop_0, cdf_cop_1)

        # get the names of the Tf and the gene from the 1st and 2nd colum of the dataframe
        genenames = df_cluster1.columns.tolist()

        # retrun the gene pair along with the ks_distance between the 2 copulas
        return genenames[0], genenames[1], differential_expression

    def compute_edges(self):
        # add progressbar with the tqdm package

        with tqdm(total=(len(self.transcription_factors)) * (len(self.genes) - 1)) as pbar:
            # compare each transcription factor against each gene

            # multiprocess for better speed
            with concurrent.futures.ProcessPoolExecutor() as executor:

                # list to store factor gene and future object
                results = []

                for factor in self.transcription_factors:
                    # remove factor so that it isn't compared with itsself
                    self.genes.remove(factor)


                    for gene in self.genes:
                        # compute a dataframe with factor/gene in cluster0 and cluster1
                        df_factor_gene_0 = pd.concat([self.cluster0[factor], self.cluster0[gene]], axis=1)
                        df_factor_gene_1 = pd.concat([self.cluster1[factor], self.cluster1[gene]], axis=1)

                        # create process that executes compute diff exp with the 2 data frames
                        f1 = executor.submit(self.compute_differential_coexpression, df_factor_gene_0, df_factor_gene_1)
                        # store respective factor gene and future object in results
                        results.append(f1)

                    # reappend previously removed factor
                    self.genes.append(factor)

                    # use results of finished processes
                for f in concurrent.futures.as_completed(results):
                    transcriptionfactor, gene, ks_dist = f.result()

                    edge = Edge(transcriptionfactor, gene, ks_dist)
                    # append edge object
                    self.edges.append(edge)
                    pbar.update(1)

    def profilingfunction(self):
        """This function is identical to compute edges except for the multiprocessing, in order to find bottlenecks"""
        with tqdm(total=(len(self.transcription_factors)) * (len(self.genes) - 1)) as pbar:
            # compare each transcription factor against each gene

            # list to store factor gene and future object
            results = []

            for factor in self.transcription_factors:
                # remove factor so that it isn't compared with itsself
                self.genes.remove(factor)
                # multiprocess for better speed 16x

                for gene in self.genes:
                    # compute a dataframe with factor/gene in cluster0 and cluster1
                    df_factor_gene_0 = pd.concat([self.cluster0[factor], self.cluster0[gene]], axis=1)
                    df_factor_gene_1 = pd.concat([self.cluster1[factor], self.cluster1[gene]], axis=1)

                    # create process that executes compute diff exp with the 2 data frames
                    results.append(tuple(self.compute_differential_coexpression(df_factor_gene_0, df_factor_gene_1)))
                    pbar.update(1)

                # reappend previously removed factor
                self.genes.append(factor)

                # use results of finished processes

    def generaterandomKsdists(self, randomsamples):

        """This function computes random Ks_distances between 2 arbitrary genes"""

        allgenes = self.cluster0.columns.tolist()
        results = []
        futureobjs = []

        # add a progressbar
        with tqdm(total=randomsamples) as pbar:

            with concurrent.futures.ProcessPoolExecutor() as executor:

                # draw 10 000 random genepairs
                for i in (range(randomsamples)):
                    # chose 2 random genes
                    gene1 = random.choice(allgenes)
                    gene2 = random.choice(allgenes)

                    # create df with only the 2 genes in it
                    df_factor_gene_0 = pd.concat([self.cluster0[gene1], self.cluster0[gene2]], axis=1)
                    df_factor_gene_1 = pd.concat([self.cluster1[gene1], self.cluster1[gene2]], axis=1)

                    f1 = executor.submit(self.compute_differential_coexpression, df_factor_gene_0, df_factor_gene_1)
                    futureobjs.append(f1)

                for res in concurrent.futures.as_completed(futureobjs):
                    _, _, ksdist = res.result()
                    results.append(ksdist)
                    # update the progressbar
                    pbar.update(1)


        outputfilename = (str(randomsamples) + 'randomKsdists.npy')
        with open(outputfilename, 'wb') as f:
            np.save(f, np.array(results))

    def generaterandomKsdistanceswithtfs(self, randomsamples, varyinggenes):

        """This function computes random Ks_distances between a random transcription factor and a gene.
           If varyinggenes is enabled then 2 different genes are used to compute the copulas.
           It's recommended to use with a df with shuffled gene expression data where the tfs are unchanged."""

        # randomsamples = how many random samples shall be drwan
        # varyinggenes = shall the tf gene pair be the same in both clusters or shall 2 different genes be used

        allgenes = self.cluster0.columns.tolist()
        results = []
        futureobjs = []

        with tqdm(total=randomsamples) as pbar:

            with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:

                # draw 10 000 random genepairs
                for i in (range(randomsamples)):
                    # chose 2 random genes
                    transcriptionfactor = random.choice(self.transcription_factors)
                    gene2 = random.choice(allgenes)

                    if varyinggenes:
                        # use 2 different genes
                        gene3 = random.choice(allgenes)
                        df_factor_gene_0 = pd.concat([self.cluster0[transcriptionfactor], self.cluster0[gene2]], axis=1)
                        df_factor_gene_1 = pd.concat([self.cluster1[transcriptionfactor], self.cluster1[gene3]], axis=1)
                    else:
                        # use the same gene in each cluster
                        # create df with only the 2 genes in it
                        df_factor_gene_0 = pd.concat([self.cluster0[transcriptionfactor], self.cluster0[gene2]], axis=1)
                        df_factor_gene_1 = pd.concat([self.cluster1[transcriptionfactor], self.cluster1[gene2]], axis=1)

                    f1 = executor.submit(self.compute_differential_coexpression, df_factor_gene_0, df_factor_gene_1)
                    futureobjs.append(f1)

                for res in concurrent.futures.as_completed(futureobjs):
                    _, _, ksdist = res.result()
                    results.append(ksdist)
                    # update the progressbar
                    pbar.update(1)

        print(results)
        if varyinggenes:
            outputfilename = (str(randomsamples) + 'tfnormal_varygene_randomKsdists.npy')
        else:
            outputfilename = (str(randomsamples) + 'tfnormal_randomKsdists.npy')

        with open(outputfilename, 'wb') as f:
            np.save(f, np.array(results))

    def edges_to_hdf(self, outputfilename):

        name = ntpath.basename(outputfilename)

        if self.random:
            name = name + '_with_random_sampling'

        name = name + '_beta_average'
        name = 'Ordered_' + name

        # create the dataframe
        df_list = []
        for edge in self.edges:
            df_list.append([edge.start, edge.end, edge.probability])

        edges_df = pd.DataFrame(df_list, columns=['Tf', 'Gene', 'Ks_dist'])
        edges_df.sort_values(by=['Ks_dist'], inplace=True, ascending=False)
        edges_df.to_hdf(name + '.h5', 'df', index=False)

    def edges_to_csv(self, outputfilename):
        name = outputfilename;

        if self.random:
            name = name.split('.')
            name.insert(-2, '_with_random_sampling')
            name = ''.join(name)

        with open(name, 'w+') as f:
            writer = csv.writer(f, delimiter=',')
            for edge in self.edges:
                writer.writerow([edge.start, edge.end, edge.probability])

        copula_scores = pd.read_csv(name, names=['Tf', 'target', 'score'])
        copula_scores.sort_values(by=['score'], inplace=True, ascending=False)
        os.remove(name)
        copula_scores.to_csv('Ordered_' + name, index=False, header=False)


tfs = []

with open('../Human_allTFs.txt', 'r') as f:
    tfs = f.readlines()

tfs = [tf.strip() for tf in tfs]
# test = Network_Inference('Data_h5/cluster0.h5', 'Data_h5/cluster1.h5', tfs, False)
# test.compute_edges()
# test.edges_to_csv('Coupulapred.csv')
tfs = ['Ugp2', 'Ddit3', 'Nedd8', 'Hes6', 'Sp3', 'Sertad2', 'Cdkn2c', 'Ostf1', 'Smarca4', 'Rab18', 'Shprh', 'Baz2a', 'Brd8',
       'Ing1', 'Bbx', 'Rab8a', 'Stat3', 'Zfp292', 'Pqbp1', 'Bclaf1', 'Irf1', 'Chd4', 'Nrf1', 'Thrap3', 'Irf3', 'Mef2a', 'Ppp1r10',
       'Snrpb', 'Znrd1', 'Hcfc1', 'Ikbkb', 'Ssrp1', 'Klf2', 'Ash1l', 'Bdp1', 'Stag1', 'Fli1', 'Nfx1', 'Ankhd1', 'Nfatc3', 'Tbx21', 'Mapk1',
       'Cnot8', 'Drg1', 'Nsd1', 'Irf7', 'Zfp148', 'Zfx', 'Irf2', 'Taf6l', 'Nr3c1', 'Pias1', 'Dek', 'Rbbp7', 'Taf7', 'Mtpn', 'Trp53',
       'Ankrd10', 'Lsm4', 'Rest', 'Ctcf', 'Pbxip1', 'Dazap2', 'Taf12', 'Gtf2e2', 'Prkar1a', 'Nfyb', 'Papola', 'Nmi', 'Uhrf2', 'Nr4a2',
       'Arid5b', 'Ccnt2', 'Mnat1', 'Zfp207', 'Ifnar2', 'Baz1b', 'Plrg1', 'Elf2', 'Klf3', 'Zfp110', 'Khdrbs1', 'C1d', 'Ubtf', 'Elf1', 'Hif1a',
       'Setbp1', 'Rora', 'Arid4a', 'Hsbp1', 'Sra1', 'Stat4', 'Strap', 'Dnmt1', 'Creb1', 'Rela', 'Tcf7', 'Fus', 'Taf3', 'Cdkn2d', 'Ash2l', 'Zbtb1',
       'Ccnh', 'Zmynd11', 'Psmc5', 'Atf7ip', 'Smarca2', 'Psmc3', 'Tcerg1', 'Zfp68', 'Rbl2', 'Sin3b', 'Elk4', 'Sirt7', 'Mllt10', 'Mtf2', 'Dido1',
       'Yaf2', 'Polr2b', 'Helz', 'Ezh2', 'Baz1a', 'Xbp1', 'Lef1']

Inferksdists = Network_Inference('../Data_h5_TCellPaper/cluster2Tcell.h5',
                                 '../Data_h5_TCellPaper/cluster5Tcell.h5', tfs, True)

# Inferksdists.generaterandomKsdistanceswithtfs(10_000,True)


Inferksdists.compute_edges()

Inferksdists.edges_to_hdf('CopluaKSdist_Tcell')
