import os
import json
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

from sklearn.ensemble import IsolationForest



###########################################
#               FUNCTIONS                 #
###########################################


def clean_reference(ref,outliers):
    for i in outliers:
        ref = ref.drop(labels=i,axis=1)

    return ref


def openJson(path,n):
    """
    Opens json files in path to create a reads matrix
    """
    tmp = os.listdir(path)
    tmp = np.array(tmp)[np.array([bool(re.findall("depths.json$",tmp[i])) for i in range(len(tmp))])]
    reads = np.zeros((n,len(tmp)))
    amplicons = ["" for x in range(n)]
    q=0
    for p in tmp:
        with open(path+p) as json_file:
            data = json.load(json_file)
            for i in range(n):
                amplicons[i] = data[i]['name']
                for j in data[i]['depths']:
                    reads[i,q] = int(data[i]['depths'][j]['min'])
        q=q+1
    reads = pd.DataFrame(data = reads,index=amplicons)
    reads.columns = [i.split('_')[0]+'_'+i.split('_')[1] for i in tmp]
    return reads

def sumLibraries(reads):
    samples = np.unique([i.split('_')[0] for i in reads.columns])
    reads_f = np.zeros((reads.shape[0],len(samples)))
    q=0
    for i in samples:
        sub = reads.filter(regex="^"+i)
        reads_f[:,q] = sub.sum(axis=1)
        q=q+1
    reads_f = pd.DataFrame(data = reads_f,index=reads.index)
    reads_f.columns = list(samples)
    return(reads_f)

def correctIndex(reads,correspondance):
    l = ["" for x in range(len(reads.index))]
    q=0
    for i in reads.index:
        l[q] = i[(len(i)-9):len(i)]
        q=q+1
    final = reads[[correspondance["amplicon"][0] in l[x] for x in range(len(l))]]
    for i in correspondance["amplicon"]:
        if i!=correspondance["amplicon"][0]:
            final = pd.concat([final,reads[[i in l[x] for x in range(len(l))]]])
    final.index = correspondance["gene_exon"] + "_" + correspondance["amplicon"]
    return final


def filterReads(reads,N,output_path=None):
    reads = reads.filter(regex="^(?!MSI)",axis=0)
    reads = reads.filter(regex="^(?!TN)")
    reads = reads.filter(regex="^(?!TP)")
    reads = reads.filter(regex="^(?!HD)")
    reads = reads.filter(regex="^(?!H2)")
    col = reads.columns
    reads = reads.loc[:,reads.sum(axis=0)>N]
    filtered_samples = col[~np.in1d(col,reads.columns)]
    if output_path is not None:
        reads.to_csv(output_path, sep="\t",index=True)
    return(reads, filtered_samples)


def normalizeReads(reads,output_path=None):
    reads_norm=reads/reads.sum(axis=0)
    if output_path is not None:
        reads_norm.to_csv(output_path+"/reads_norm.tsv", sep="\t",index=True)
    return(reads_norm)


def aberrantSamples(reads,conta='auto'):
    tmp = np.percentile(reads, 99, axis = 0)/np.mean(reads, axis = 0)
    random_data = np.array(tmp).reshape(-1,1)

    clf = IsolationForest(contamination=conta).fit(random_data)
    preds = clf.predict(random_data)
    res = np.array(reads.columns)[preds==-1]
    return(res)


def aberrantAmpliconsPerSample(name,reads_norm,verbose=False):
    random_data = np.array(reads_norm[name]).reshape(-1,1)
    clf = IsolationForest(contamination=0.01).fit(np.array(np.mean(reads_norm, axis = 1)).reshape(-1,1))
    preds = clf.predict(random_data)
    if verbose:
        print(name)
        print(np.array(reads_norm.index)[preds==-1])
    return(np.array(reads_norm.index)[preds==-1])



def amplifEvalGene(reads,abSamples,gene,sample):
    reads_m = reads/reads.median(axis=0)
    reads_m = reads_m.filter(regex="^"+gene,axis=0)
    sub = reads_m
    for i in abSamples:
        sub = sub.drop(labels=i,axis=1)
    reads_m = reads_m[sample]
    val = np.mean(reads_m)/np.mean(sub.mean())
    if val==np.inf:
        val = 100
    return val


def scoreAmplif(k,n,N,mu):
    p = n/N
    x = np.log(1/((p**k)*(1-p)**(n-k)))*(k/n)
    return x

def aberrantAmpliconsFinal(reads, reads_norm, abSamples,abSamples2,run,threshold):
    f = pd.DataFrame(columns=["run","name","gene","ratio","score"])

    q=0
    for name in abSamples2:
        abAmp = aberrantAmpliconsPerSample(name,reads_norm)
        if abAmp.shape!=(0,):
            genes = np.unique([i.split('_')[0] for i in abAmp])
            for gene in genes:
                r = re.compile(gene)
                abEx = list(filter(r.match, abAmp))
                exons1 = [i.split('_')[0]+"_"+i.split('_')[1] for i in abEx]
                tmp = reads.filter(regex="^"+gene,axis=0)
                exons2 = [i.split('_')[0]+"_"+i.split('_')[1] for i in tmp.index]

                score = scoreAmplif(len(abEx),tmp.shape[0],reads.shape[0],len(abEx)/tmp.shape[0])

                amplif = amplifEvalGene(reads_norm, abSamples, gene, name)

                if score>threshold and amplif>1.2:
                    f.loc[q] = [run,name,gene,amplif,score]
                    q=q+1


    return(f)



###########################################
#               MAIN                      #
###########################################
parser = argparse.ArgumentParser(description='Cette fonction sert a detecter les amplifications de MET dans les patients atteints de cancer du poumon')
parser.add_argument('-r', '--run', type=str, help='Numero du run')
args = parser.parse_args()

PATH = "/mnt/Bioinfo/BioTS/Results//ADIVaR/FDG_Juno_v3/"
output_path = "/mnt/chu-ngs/Labos/BioTS/SOMAT/DIAG/Juno/MET_amp/"
correspondance = pd.read_csv("/mnt/Bioinfo/BioTS/Projets/CNV/correspondance.txt",sep="\t")
run = args.run

if run is None:
   print("-r has no default value")
   ld = []
else:
   ld = os.listdir(PATH)
   ld = np.array(ld)[np.array([bool(re.findall(run,ld[i])) for i in range(len(ld))])]

N = 1536 #number of amplicons

res = pd.DataFrame(columns=["run","name","gene","ratio","score"])
k = 0

if len(ld)>0:
    path = PATH+ld[0]+'/data/'
    reads = openJson(path,N)
    reads = sumLibraries(reads)
    reads = correctIndex(reads,correspondance)
    final, filtered_samples = filterReads(reads,N*200,output_path=output_path+"reads_"+args.run+".tsv")
    q=0
    #filtered_samples = filtered_samples[[bool(re.search("^P", i)) for i in filtered_samples]]
    if len(filtered_samples)>0:
        for i in filtered_samples:
            tmp = [run,i,"-","Non Analysable","-"]
            res.loc[q] = tmp
            q=q+1
    final_norm = normalizeReads(final)
    abSamples = aberrantSamples(final)
    #allSamples = final.filter(regex="^P").columns

    ff = aberrantAmpliconsFinal(final,final_norm,abSamples,final.columns,run,1)
    if ff.shape[0]>0:
        res = res.append(ff)

    #res = res.loc[np.in1d(res["gene"],"MET"),:]
    res.index=range(res.shape[0])
    #negSamples = allSamples[~np.in1d(allSamples,res["name"])]
    #q = res.shape[0]
    #if len(negSamples)>0:
    #    for i in negSamples:
    #        q=q+1
    #        res.loc[q] = [run,i,"MET","Negatif","-"]

    res.to_csv(output_path+"CNV_Juno_all_"+run+".tsv", sep="\t",index=False)

else:
    print("Erreur. Verifier le nom du run.")
