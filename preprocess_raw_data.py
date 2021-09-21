# import matplotlib
# matplotlib.use('agg')
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 14:13:51 2020

This module:
    * Cleans a raw dataset
    * Preprocessesthe data 
        - filter sparse genes with 1/3 missing value in its gene expression, 
        i.e., keep genes whose time series have >= 2/3 of non-zeros.
    * Plots correlation distribution using 10% data sample
    * Finds the correlation value at 95% Quantile 
        
@author: Shuyue Xue
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def configure_SLV_dataframe():
    """
    Configure data files obtained from the published datasets online. 
    
    parameters
    -------------   
    data_file : string
        name of the data file.
    time_file : string
        name of the time information file.

    Returns
    -------------
    data_df : Dataframe
    """
    h1_data = 'Hourly1TimeSeries_SLV_KallistoNormedGeneGencodeGC.csv'
    h2_data = 'Hourly2TimeSeries_SLV_KallistoNormedGeneGencodeGC.csv'
    time_file = 'TimesHourly.csv'
    data_files = [h1_data, h2_data]
    
    h1, h2 = [pd.read_csv(file, index_col=0, header=None) for file in data_files]
    time_points = pd.read_csv(time_file, header=None)
    
    # Add time points as the columns
    h1.columns = round( time_points.iloc[:, 0]*len(time_points) )
    h1.columns.name = 'Hour'
    h1.index.name = 'genes'
    
    h2.columns = round( time_points.iloc[:, 0]*len(time_points) )
    h2.columns.name = 'Hour'
    h2.index.name = 'genes'
    
    return h1, h2


def configure_Bcell_datasets():
    """
    Create a tupe of dataframes from (B cell) KallistoNormed datafile. 

    Returns
    -------
    drt : data frame
        Dataframe for the Reapted Treated data.
    dru : data frame
        Dataframe for the Reapted Untreated data.
    dt : dataframe
        Dataframe for the Treated data.
    du : Dataframe
        Dataframe for the Untreated data.

    """
    data_file = 'RTX_KallistoNormedGeneGencodeGC.csv'

    # Configure the dataframes
    data_df = pd.read_csv(data_file, index_col=0, header=0)
    data_df.columns.name = 'time'
    data_df.index.name = 'genes'
    drt, dru, dt, du = [ data_df.iloc[:, i:i+6] \
                        for i in range(data_df.shape[1]) if 0 == (i%6) ]
    for df in [drt, dru, dt, du]:
        df.columns = [0, 1, 2, 4, 7, 15]

    return drt, dru, dt, du


def clean_raw_data(df):
    """
    Quality Control of the raw dataest.
    
    parameters
    -------------   
    df : data frame
        The dataframe of the raw data. 
        Must have the index of genes & column of time points.
    
    Returns
    -------------
    df : data frame
        A data frame wihout constant genes or duplicates. 
    """

    # # Check missing values
    # print('Check the missing values in the dataframe .... ')
    # number_of_missing_values = df.isnull().sum().sum()
    # print('# of missing values : ', number_of_missing_values)

    # Delete genes with 0 expressions at all time points (all-0 genes)
    df_clean = df.loc[ (df!=0).any(axis=1) ]
    # num_of_all_0_genes = df.shape[0] - df_clean.shape[0]
    # print('# of all-0 genes: ', num_of_all_0_genes)

    # Check the number of constant genes:
    # num_of_constant_genes = (df_clean.var(axis=1) == 0).sum()
    # print('After removing all-0-genes -------------')
    # print('num of constant genes = ', num_of_constant_genes)
    
    # Delete constant genes:
    df_clean = df_clean.loc[ df_clean.var(axis=1) != 0 ]
    # print('shape of clean data frame', df_clean.shape)

    # # Check the duplicate genes: True if all value counts == 1
    # print('no duplicates? ', (df_clean.index.value_counts() == 1 ).all() )

    # print('Finished Data QC ----')
    return df_clean


def preprocess_data(df, missing_value_threshold=1/3,):
    """
    Prof. George's preprocessing algorithm.
    
    parameters
    -------------   
    df : data frame
        The clean data frame results from the clean_raw_data() method.
    missing_value_threshold : float, default 1/3 
        The fraction of the missing values allowed in the time series 
        in defining the sparse genes.
    
    Returns
    -------------
    df : data frame
        Clean data frame
    """
    # print('\nStart preprocessing ....')
    # Step 1. set 0s to be NaN:
    df.mask( df==0.0, inplace=True )

    # Step 2. set values <1 to be 1
    df.mask( df<1, other=1, inplace=True )

    # # Check the number of constant genes with all-1s resulted from step 2
    # num_of_constant_genes = (df.var(axis=1) == 0).sum()
    # print('# of constant genes (all-1 genes) = ', num_of_constant_genes)

    # Step 3. remove the constant genes again
    df = df.loc[ df.var(axis=1) != 0 ]
    # num_of_constant_genes = (df.var(axis=1) == 0).sum()
    # print('No more constant genes? ', num_of_constant_genes == 0)

    # Step 4. Filter sparse genes
    number_of_non_NaNs_threshold = \
    df.shape[1] * (1 - missing_value_threshold)

    df.dropna( thresh=number_of_non_NaNs_threshold, inplace=True )
    # print('Shape of the dense dataframe : ', df.shape)
    # print('Finished preprocessing ----')

    return df


def get_pairwise_fold_change_genes(T, U,
         q_fold=1/3, distribution=True, title=''):
    """
    Select genes when the mean of the pairwise diff/sum ratio, 
    < (x_t - x_u) / (x_t + x_u) > is above/below the +/- quantile cutoff.
    
    parameters
    -------------   
    df1, df2 : data frame
        The gene data frame results from the preprocess_data() method.
    
    Returns
    -------------
    high_change_genes : set
        A set of gene ids with high fold changes. 
    """
    # print('\nStart selecting high changed genes....')
    # 1. Calculate the (x_t - x_u) / (x_t + x_u)
    ratio = (T - U) / (T + U)

    # check the all-NaNs:
    # print( 'All-NaN genes?: \n', ratio.isnull().all(axis=1) )
    # print( '# All-NaN genes: ', ratio.isnull().all(axis=1).sum() )
    # remove the all-NaNs:
    ratio.dropna(axis=0, how='all', inplace=True)
    # print( '# All-NaN genes after clean up: ',
    #       ratio.isnull().all(axis=1).sum() )

    # 2. Find the mean in each row in the ratio :
    ratio_mean = ratio.mean(axis=1)

    # Plot the distribution
    if distribution: # Plot the distribution of the ratio
        print('Plotting Ratio distribution....')
        plt.clf()
        bin_num = int( round( np.sqrt (ratio_mean.shape[0])) )

        ratio_mean.plot(kind='hist', bins=bin_num,
                   title=title+' Fold Change Dist',
                   label=title)
        plt.legend()
        plt.savefig(title+'fold change dist.pdf')

    # Select genes using quantile cutoff
    cutoff_high = ratio_mean.quantile(q = q_fold)
    cutoff_low = ratio_mean.quantile(q = round(1-q_fold, 3) )
    # print('right at %s quantile = '%round(q_fold, 3),  cutoff_high)
    # print('left at  %s quantile = '% round(1-q_fold, 3), cutoff_low)

    # high_change_genes = ratio_mean[ (ratio_mean >= cutoff_high) ^
    #                           (ratio_mean <= cutoff_low) ]

    high_change_genes = ratio_mean[ (ratio_mean >= cutoff_high) |
                              (ratio_mean <= cutoff_low) ]
    # print('# of selected genes = ', high_change_genes.shape[0])
    # print('End selecting genes----')
    return set(high_change_genes.index)


if __name__ == '__main__':
    import os
    # ------------- Parameters ---------------------------------------
    NaN_in_sparse_genes = 1/3   # Missing Values in sparse genes
    q = 0.95                    # Quantile for the cut-off correlation
    #-----------------------------------------------------------------

    # Load the dataset: H1 doesn't have vaccine, H2 HAS vaccine
    os.chdir('./SLV_data')
    h1_data = 'Hourly1TimeSeries_SLV_KallistoNormedGeneGencodeGC.csv'
    time = 'TimesHourly.csv'
    df = configure_SLV_dataframe(h1_data, time).iloc[:1313]

    # Output results into Results Directory
    os.chdir('../')
    if not os.path.exists('SLV_results'):
            os.makedirs('SLV_results')
    os.chdir('./SLV_results')

    df = clean_raw_data(df)
    df = preprocess_data(df)