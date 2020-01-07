#!/usr/bin/env python3

import numpy as np
import pandas as pd
from functools import reduce
import re

columns = ['INDEX', 'X_AXIS', 'Y_AXIS', 'DX', 'DY', 'X_PPM', 'Y_PPM', 'X_HZ', 'Y_HZ', 'XW', 'YW', 'XW_HZ',
           'Y_HZ', 'X1', 'X3', 'Y1', 'Y3', 'HEIGHT', 'DHEIGHT', 'VOL', 'PCHI2', 'TYPE', 'ASS', 'CLUSTID', 'MEMCNT'] #Column header of input data

t1=pd.read_csv("test_t1.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t2=pd.read_csv("test_t2.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t3=pd.read_csv("test_t3.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t4=pd.read_csv("test_t4.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t5=pd.read_csv("test_t5.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t6=pd.read_csv("test_t6.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t7=pd.read_csv("test_t7.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t8=pd.read_csv("test_t8.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t9=pd.read_csv("test_t9.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t10=pd.read_csv("test_t10.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t11=pd.read_csv("test_t13.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t12=pd.read_csv("test_t14.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])
t13=pd.read_csv("test_t15.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])


df_list = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13] #Collect dataframes into a list

def extract(df, y):  # Process each dataframe to extract assignment and desired column (e.g., extract (df, 'X_PPM'))
    df = df[~df.ASS.str.contains("None")]  # drop rows with 'None' in assignment 'ASS' column of an nmrDraw peak.tab
    df = df[['ASS', y]]  # y = 'X_PPM', 'Y_PPM', 'Z_PPM', etc.
    df['sort'] = df['ASS'].str.extract('(\d+)', expand=False).astype(int)  # for sorting ASS column (resID-Num: M4)
    df.sort_values(by=['sort'], inplace=True, ascending=True)
    df = df.drop('sort', axis=1)  # drop 'sort' column, used for sorting the dataframe by ASS
    return df

def merge(dataframes):  # merge list of dataframes
    df = reduce(lambda left, right: pd.merge(left, right, on=['ASS'], how='outer'),
                dataframes)  # merge all based on 'ASS'
    res = df.iloc[:, 0]  # separate residue name-num into an array
    data = df.iloc[:, 1:14]  # separate data into an array
    resn_split = [re.split(r'(\d+)', res)[0:2] for res in
                  res]  # split resID-Num (e.g., E4) col into two cols (str, num)
    df = pd.DataFrame(np.concatenate([resn_split, data], axis=1))  # combine arrays for res name, num, data into one df
    df.columns = (['ID', 'Num', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12','t13'])  # name cols
    return df


df1H = merge([extract(df, 'X_PPM') for df in df_list])  # extract desired columns from each dataframe and merge them according to ASS

df15N = merge([extract(df, 'Y_PPM') for df in df_list])  

#Save Extracted, Merged, Aligned and Sorted Chemical Shift Table to CSV
df1H.to_csv('1H_CS_Table.csv', index=False, float_format='%.2f')

df15N.to_csv('15N_CS_Table.csv', index=False, float_format='%.2f')
