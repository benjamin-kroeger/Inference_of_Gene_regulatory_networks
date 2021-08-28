import numpy as np
import pandas as pd
import csv
import scipy.stats as stats
import glob, os



class merge_bootstraps:

    """This class contains functions to merge multiple bootstrapps into one file"""

    def __init__(self, pathtofolder,filetype):

        # stores the different bootstrapped dataframes
        self.bootstarp_dfs = []
        counter = 0

        if filetype == 'csv':
            # read all the csv files in the specified Folder
            for file in os.listdir(pathtofolder):
                if file.endswith('.csv'):
                    # make sure each score columm has a unique name using counter
                    df = pd.read_csv(os.path.join(pathtofolder, file), names=['Tf', 'target', 'score' + str(counter)])
                    # safe created df
                    self.bootstarp_dfs.append(df)
                    counter += 1
        elif filetype == '.h5':

            # read all the h5 files in the specified Folder
            for file in os.listdir(pathtofolder):
                if file.endswith('.h5'):
                    # make sure each score columm has a unique name using counter
                    df = pd.read_hdf(os.path.join(pathtofolder, file), 'df')
                    df.columns = ['Tf', 'target', 'score' + str(counter)]
                    # safe created df
                    self.bootstarp_dfs.append(df)
                    counter += 1


    def meanabserror(self, ks_scores):

        mean = np.mean(ks_scores)
        errors=[]

        for dist in ks_scores:

            errors.append(abs(mean-dist))

        return sum(errors)/len(ks_scores)

    def merge_bootstarps(self):

        # combine the first 2 df into one
        merged_df = pd.merge(self.bootstarp_dfs[0], self.bootstarp_dfs[1], on=['Tf', 'target'])

        # merge the rest of the data frames into one each time adding another score collum
        for i in range(2, len(self.bootstarp_dfs)):
            merged_df = pd.merge(merged_df, self.bootstarp_dfs[i], on=['Tf', 'target'])

        # after merging all df one gets a merged_df with tf // gene // score1 // score2 // ...
        meanscores = []
        standard_errors = []

        for i in range(0, len(merged_df)):
            # get i-th row of the merged_df to numpy converts it into a 1dim array with gene and scores
            # only take the values from the score columns
            ks_scores = merged_df.iloc[i].to_numpy()[2:]

            # calculate the means and absolute errors and store them
            mean = (np.mean(ks_scores))
            meanabserror = self.meanabserror(ks_scores)
            meanscores.append(float(mean))
            standard_errors.append(float(meanabserror))

        # delete the score collumns from the merged_df
        merged_df.drop(merged_df.iloc[:, 2:], inplace=True, axis=1)

        # add the mean and error as 2 new columns
        merged_df.insert(2, 'Ks_dist_score', meanscores, True)
        merged_df.insert(3, 'Mean_abserror', standard_errors, True)

        # sort by ks dist axis=0 sind zeilen
        merged_df.sort_values(by=['Ks_dist_score'], ascending=False)

        merged_df.to_hdf('CopulaKSdist_beta_average_with70percent_randomsampling_MAE.h5','df')


test = merge_bootstraps('../Gene_correlation_data/Bootstraps_with_random_sampling/Samplesize70percent','.h5')
test.merge_bootstarps()
