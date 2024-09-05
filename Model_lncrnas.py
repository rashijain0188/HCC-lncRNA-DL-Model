# -*- coding: utf-8 -*-


import rpy2
import numpy as np
import pandas as pd
import scipy.stats as stats
import warnings
import keras
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from tkinter import Tk
from tkinter.filedialog import askopenfilename
limma = importr('limma')
edgeR = importr('edgeR')

# Function to prompt the user to select a file
print("\n Upload the lncRNA counts data for the samples with lncRNAs as rows and patient samples as columns : \n")
def upload_data():
    # Hide the root window of tkinter
    Tk().withdraw()
    
    # Open the file dialog and ask the user to select a file
    file_path = askopenfilename(title="Select a data file", filetypes=[("CSV files", "*.csv")])
    
    # Check if a file was selected
    if file_path:
        # Load the selected file into a DataFrame
        df = pd.read_csv(file_path)
        print("Data loaded successfully!")
        return df
    else:
        print("No file selected.")
        return None

# Call the function and get the DataFrame
df = upload_data()
#df = df.set_index(['IDs'])
df= df.set_index(['Unnamed: 0'])

# Display the first few rows of the DataFrame if data was loaded
if df is not None:
    print(df.head())

# Activate the automatic conversion of pandas DataFrames to R data.frames
pandas2ri.activate()

def edger_calcnormfactors(counts_df, ref=None, logratio_trim=0.3,
                          sum_trim=0.05, acutoff=-1e10, verbose=False):

    # discard genes with all-zero counts
    Y = counts_df.values.copy()
    allzero = np.sum(Y>0,axis=1)==0
    if np.any(allzero):
        Y = Y[~allzero,:]

    # select reference sample
    if ref is None:  # reference sample index
        f75 = np.percentile(Y/np.sum(Y,axis=0), 75, axis=0)
        ref = np.argmin(np.abs(f75-np.mean(f75)))
        if verbose:
            print('Reference sample index: '+str(ref))

    N = np.sum(Y, axis=0)  # total reads in each library

    # with np.errstate(divide='ignore'):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # log fold change; Mg in [1]
        logR = np.log2((Y/N).T / (Y[:,ref]/N[ref])).T
        # average log relative expression; Ag in [1]
        absE = 0.5*(np.log2(Y/N).T + np.log2(Y[:,ref]/N[ref])).T
        v = (N-Y)/N/Y
        v = (v.T + v[:,ref]).T  # w in [1]

    ns = Y.shape[1]
    tmm = np.zeros(ns)
    for i in range(ns):
        fin = np.isfinite(logR[:,i]) & np.isfinite(absE[:,i]) & (absE[:,i] > acutoff)
        n = np.sum(fin)

        loL = np.floor(n*logratio_trim)+1
        hiL = n + 1 - loL
        loS = np.floor(n*sum_trim)+1
        hiS = n + 1 - loS
        rankR = stats.rankdata(logR[fin,i])
        rankE = stats.rankdata(absE[fin,i])
        keep = (rankR >= loL) & (rankR <= hiL) & (rankE >= loS) & (rankE <= hiS)
        # in [1], w erroneously defined as 1/v ?
        tmm[i] = 2**(np.nansum(logR[fin,i][keep]/v[fin,i][keep]) / np.nansum(1/v[fin,i][keep]))

    tmm = tmm / np.exp(np.mean(np.log(tmm)))
    return tmm

def voom_transform(counts_df):
    """Apply counts transformation from limma-voom"""
    lib_size = counts_df.sum(0)
    norm_factors = edger_calcnormfactors(counts_df)
    return np.log2((counts_df + 0.5) / (lib_size*norm_factors + 1) * 1e6)

df2= voom_transform(df)
df3=df2.T

NT_model = keras.models.load_model("D:/Rashi/lncRNAs 24.07.24/Normal vs Tumor/NT3.keras")
NT_model.summary()
df4 = NT_model.predict(df3)
df4 =np.argmax(df4, axis = 1)
mapping = {0: "Normal", 1: "Tumor"}
text_labels = np.vectorize(mapping.get)(df4)
print(text_labels)

# Prepare a list to hold results
results = []

# Ensure the length of text_labels matches the number of indices in df3
if len(text_labels) != len(df3.index):
    raise ValueError("The length of text_labels must match the number of indices in df3.")

# Iterate over the index of df3
for idx, label in zip(df3.index, text_labels):
    results.append({
        'Index': idx,
        'Sample': label
    })
# Convert the list to a DataFrame
results_df = pd.DataFrame(results)
print(results_df)

df5 = df3[df4 == 1]
print(df5)


EA_model = keras.models.load_model("D:/Rashi/lncRNAs 24.07.24/Early vs Advanced/EA2.keras")
df6 = EA_model.predict(df5)
df6 =np.argmax(df6, axis = 1)
mapping = {0: "Advanced", 1: "Early"}
text_labels2 = np.vectorize(mapping.get)(df6)
print(text_labels2)

# Prepare a list to hold results
results2 = []
if len(text_labels2) != len(df5.index):
    raise ValueError("The length of text_labels must match the number of indices in df3.")
for idx, label in zip(df5.index, text_labels2):
    results2.append({
        'Index': idx,
        'Sample': label })
results_df2 = pd.DataFrame(results2)
print(results_df2)

df7 = df5[df6 == 0]
print(df7)

stages_model = keras.models.load_model("D:/Rashi/lncRNAs 24.07.24/Stages/stages7.keras")
df8 = stages_model.predict(df7)
df8 =np.argmax(df8, axis = 1)
mapping = {0: "Stage II", 1: "Stage III", 2: "Stage IV"}
text_labels3 = np.vectorize(mapping.get)(df8)
print(text_labels3)

# Prepare a list to hold results
results3 = []
if len(text_labels3) != len(df7.index):
    raise ValueError("The length of text_labels must match the number of indices in df3.")
for idx, label in zip(df7.index, text_labels3):
    results3.append({
        'Index': idx,
        'Sample': label })
results_df3 = pd.DataFrame(results3)
print(results_df3)


df_merged = pd.merge(results_df, results_df2, left_index=True, right_index=True, how='outer')
print(df_merged)

df_merged2 = pd.merge(df_merged, results_df3, left_index=True, right_index=True, how='outer')
print(df_merged2)

df_merged2.to_csv('Results.csv', index=True)
