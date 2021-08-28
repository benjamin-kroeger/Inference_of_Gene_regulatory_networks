import numpy as np
import pandas as pd

class filteringresults:

    """This class contains the function to filter coputed edges based on the previously computed error and assigns p-vlaues to the
    edges."""
    def __init__(self, pathtofile,fileending):

        # reading the data and checking if the error column has already been computed
        if fileending == '.csv':

            self.tf_gene_df = pd.read_csv(pathtofile, header=None)
            if len(self.tf_gene_df.columns) < 4:
                raise ValueError('There are not enough columns in the dataframe')

            self.tf_gene_df.columns = ['Tf','gene','score','MAE_error']

            self.filename = pathtofile.split('/')[-1].rstrip(".csv")

        elif fileending == '.h5':
            self.tf_gene_df = pd.read_hdf(pathtofile,'df')

            if len(self.tf_gene_df.columns) < 4:
                raise ValueError('There are not enough columns in the dataframe')

            self.tf_gene_df.columns = ['Tf', 'gene', 'score', 'MAE_error']

            self.filename = pathtofile.split('/')[-1].rstrip(".h5")

        self.errorfiltered = False
        self.pvaluefiltered = False


    def filterbasedonstderror(self,cutoffvalue):

        if cutoffvalue > self.tf_gene_df['Mean_abserror'].max():
            raise ValueError('Cutoff value is bigger than the largest std error')

        indextodelete = self.tf_gene_df[self.tf_gene_df['Mean_abserror'] > cutoffvalue].index

        self.tf_gene_df.drop(indextodelete,inplace=True)

        self.errorfiltered = True

    def assignpvalue(self):

        # given a file with random ks_distances this function computes the p-value for each edge

        with open('../Data_h5/10000tfnormal_varygene_randomKsdists.npy', 'rb') as f:
            randks = np.load(f)

        randks.sort()

        pvalues = []


        for ksdist in self.tf_gene_df['score']:
            # searchsorted returns index at which to insert
            pvalues.append(1-np.searchsorted(randks,ksdist)/10_000)

        print(pvalues)
        print(len(pvalues))
        self.tf_gene_df['pvalue'] = pvalues

        print(self.tf_gene_df)

        self.tf_gene_df.to_hdf(self.filename + '_pvalues.h5','df')


        self.pvaluefiltered = True



    def saveresult(self):

        outputname = ''
        if self.errorfiltered:
            outputname = outputname + '_stderror'

        if self.pvaluefiltered:
            outputname = outputname + '_pvalue'



        self.tf_gene_df.to_csv(self.filename + (outputname + '_filtered.csv'), index=None)





filtering = filteringresults('../Gene_correlation_data/CopulaKSdist_beta_average_with70percent_randomsampling_MAE.h5','.h5')
filtering.assignpvalue()
