import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

####################################################
#                  FUNCTIONS                       #
####################################################

def norm_ref(ref):
    med = ref.median(axis=0)
    norm_ref = ref/med
    return norm_ref, med

    
def create_synthetic(norm_ref,med,N):
    synt = pd.DataFrame(index=norm_ref.index,columns=range(N))
    for j in range(N):
        for i in norm_ref.index:
            synt[j][i] = np.random.choice(norm_ref.loc[norm_ref.index==i].values.flatten())

        synt[j] = synt[j] * np.random.choice(med)
    synt.columns = ["sample_" + str(i) for i in range(synt.shape[1])]
    return synt

def add_features(synt,sample,gene,factor,exon=None):

    for j in sample:
        if exon is None:
            pattern = gene
            tmp = [str(synt.index[i]).split("_")[0] for i in range(synt.shape[0])]
            synt.loc[[tmp[i]==pattern for i in range(len(tmp))],j] = synt.loc[[tmp[i]==pattern for i in range(len(tmp))]][j] * factor

        if exon is not None:
            pattern = gene + "_" + exon
            tmp = [str(synt.index[i]).split("_")[0] + "_" + str(synt.index[i]).split("_")[1] for i in range(synt.shape[0])]
            synt.loc[[tmp[i]==pattern for i in range(len(tmp))],j] = synt.loc[[tmp[i]==pattern for i in range(len(tmp))]][j] * factor

    return synt




####################################################
#                     MAIN                         #
####################################################

ref = pd.read_csv("path/to/experimental/data",sep="\t",index_col=0)

norm_ref, med = norm_ref(ref)

synt = create_synthetic(norm_ref,med, N)

synt = add_features(synt, ["sample_0","sample_1"],"EGFR",factor=10)
