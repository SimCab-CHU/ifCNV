import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import json
import re

from sklearn.ensemble import IsolationForest


#############################################
#            FUNCTIONS                      #
#############################################


def create_synthetic(norm_ref,med,N):
    synt = pd.DataFrame(index=norm_ref.index,columns=range(N))
    for j in range(N):
        for i in norm_ref.index:
            synt[j][i] = np.random.choice(norm_ref.loc[norm_ref.index==i].values.flatten())

        synt[j] = synt[j] * np.random.choice(med)
    synt.columns = ["sample_" + str(i) for i in range(synt.shape[1])]
    return synt
  
  
 def normalizeReads(reads,output_path=None,save=False):
    reads_norm=reads/reads.median(axis=0)
    reads = np.log(reads+1)
    if save==True:
        reads_norm.to_csv(output_path, sep="\t",index=True)
    return(reads_norm)
  
  
 def aberrantSamples(reads,conta='auto'):
    #reads = reads/np.sum(reads)
    
    tmp = np.percentile(reads, 99, axis = 0)/np.mean(reads, axis = 0)
    random_data = np.array(tmp).reshape(-1,1)
    clf = IsolationForest(contamination=conta).fit(random_data)
    preds = clf.predict(random_data)
    res_amp = np.array(reads.columns)[preds==-1]
    
    tmp = np.percentile(reads, 1, axis = 0)/np.mean(reads, axis = 0)
    random_data = np.array(tmp).reshape(-1,1)
    clf = IsolationForest(contamination=conta).fit(random_data)
    preds = clf.predict(random_data)
    res_del = np.array(reads.columns)[preds==-1]
    
    res = np.unique(np.concatenate((res_amp,res_del)))
    norm = reads.columns[~np.in1d(reads.columns,res)]
    
    return(res, norm)
  
  
  def aberrantAmpliconsPerSample(name,reads_norm,conta='auto',verbose=False):
    random_data = np.array(reads_norm[name]).reshape(-1,1)
    clf = IsolationForest(contamination=conta).fit(np.array(np.mean(reads_norm, axis = 1)).reshape(-1,1))
    preds = clf.predict(random_data)
    if verbose:
        print(name)
        print(np.array(reads_norm.index)[preds==-1])
    return(np.array(reads_norm.index)[preds==-1])
  
  def amplifEvalGene(reads,abSamples,gene,sample):
    reads_m = reads/reads.median(axis=0)
    sub = reads_m
    for i in abSamples:
        sub = sub.drop(labels=i,axis=1)
    reads_m = reads_m.filter(regex="^"+gene,axis=0)
    reads_m = reads_m[sample]   
    val = np.mean(reads_m)/np.mean(sub.mean())
    if val==np.inf:
        val = 100
    return val

def scoreAmplif(k,n,N):
    p = n/N
    x = np.log(1/((p**k)*(1-p)**(n-k)))*(k/n)
    # score = 1/(1+np.exp(-x))
    score = x/390 + 190/390
    
    return x
  
  
  
