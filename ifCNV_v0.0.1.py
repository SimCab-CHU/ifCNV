#
# Copyright (C) 2021
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Simon Cabello'
__copyright__ = 'Copyright (C) 2021'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 's-cabelloaguilar@chu-montpellier.fr'
__status__ = 'prod'


import os
import json
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import subprocess
import io
import argparse
from sklearn.ensemble import IsolationForest



###########################################
#               FUNCTIONS                 #
###########################################


def clean_reference(ref,outliers):
    for i in outliers:
        ref = ref.drop(labels=i,axis=1)

    return ref

def createReadsMatrix(pathToBam, bedFile, pathToBedtools, output=None, verbose=False):
    cmd = ["ls", pathToBam]
    res = subprocess.check_output(cmd)
    final=pd.DataFrame()

    for i in res.decode('utf-8').split("\n"):
        if i.endswith(".bam"):
            if verbose==True:
                print("Processing sample "+i[:-4]+"...")
            command = [
                pathToBedtools,
                "multicov",
                "-bams", pathToBam+"/"+i,
                "-bed", bedFile]

            res = subprocess.check_output(command)
            data = io.StringIO(res.decode("utf-8"))
            df = pd.read_csv(data, sep='\t',header=None)
            nam = i[:-4]
            final[nam] = df[len(df.columns)-1]
            if verbose==True:
                print(i[:-4]+" Done")
    final.index = list(df[3])

    if output is not None:
        if verbose==True:
            print("Reads matrix created !")
        final.to_csv(output,sep="\t")

    return(final)


def filterReads(reads,N,regtar=None,regsamp=None):
    col = reads.columns
    rows = reads.index
    if regtar is not None:
        reads = reads.filter(regex=regtar,axis=0)
    if regsamp is not None:
        reads = reads.filter(regex="^(?!"+regsamp+")")
    reads = reads.filter(regex="^(?!MSI)",axis=0)
    reads = reads.filter(regex="^(?!TN)")
    reads = reads.filter(regex="^(?!TP)")
    reads = reads.filter(regex="^(?!HD)")
    reads = reads.filter(regex="^(?!H2)")
    reads = reads.loc[reads.sum(axis=1)/len(reads.columns)>N,:]
    filtered_samples = col[~np.in1d(col,reads.columns)]
    filtered_targets = rows[~np.in1d(rows,reads.index)]
    return(reads, filtered_samples, filtered_targets)


def normalizeReads(reads):
    reads_norm=reads/reads.sum(axis=0)
    return(reads_norm)


def aberrantSamples(reads,conta='auto'):    
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
    norm = np.array(reads.columns[~np.in1d(reads.columns,res)])
    
    return(res, norm)



def aberrantAmpliconsPerSample(name,reads_norm,CNVneg,conta=0.01):
    random_data = np.array(reads_norm[name]).reshape(-1,1)
    norm = np.array(np.mean(reads_norm[CNVneg], axis = 1))
    clf = IsolationForest(contamination=conta).fit(norm.reshape(-1,1))
    preds = clf.predict(random_data)
    return(np.array(reads_norm.index)[preds==-1])

def scoreAmplif(k,n,N):
    p = n/N
    x = np.log(1/((p**k)*(1-p)**(n-k)))*(k/n)
    return x


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


def aberrantAmpliconsFinal(reads, reads_norm, CNVpos, CNVneg, scoreThreshold=10,ampThreshold=1.5,conta=0.01,mode="fast",run="ifCNV"):
    f = pd.DataFrame(columns=["Run","Sample name","Region","Reads ratio","Score"])
        
    if mode=="extensive":
        samples = [*CNVpos,*CNVneg]
    if mode=="fast":
        samples = CNVpos

    q=0
    for name in samples:       
        abAmp = aberrantAmpliconsPerSample(name,reads_norm,CNVneg,conta=conta)
        if abAmp.shape!=(0,):
            genes = np.unique([i.split('_')[0] for i in abAmp])
            for gene in genes:
                r = re.compile(gene)
                abEx = list(filter(r.match, abAmp))
                exons1 = [i.split('_')[0]+"_"+i.split('_')[1] for i in abEx]
                tmp = reads.filter(regex="^"+gene,axis=0)
                exons2 = [i.split('_')[0]+"_"+i.split('_')[1] for i in tmp.index]

                score = scoreAmplif(len(abEx),tmp.shape[0],reads.shape[0])
                amplif = amplifEvalGene(reads_norm, CNVneg, gene, name)

                if score>scoreThreshold and amplif>ampThreshold:
                    f.loc[q] = [run,name,gene,amplif,score]
                    q=q+1

    return(f)




###########################################
#               MAIN                      #
###########################################

parser = argparse.ArgumentParser(description='ifCNV')
parser.add_argument('-i', '--input', type=str, help='Path to the input bam folder')
parser.add_argument('-b', '--bed', type=str, help='Path to the bed file')
parser.add_argument('-t', '--bedtools', type=str, help='Path to bedtools')
parser.add_argument('-m', '--mode', type=str, default='fast' help='fast or extensive')
parser.add_argument('-min', '--minReads', type=str, default=100 help='Min mean reads per target')
parser.add_argument('-cs', '--contaSamples', default = "auto", help='Contamination parameter for the AberrantSamples function')
parser.add_argument('-ct', '--contaTargets', default = "auto", help='Contamination parameter for the AberrantTargets function')
parser.add_argument('-sT', '--scoreThreshold', type=int, default=5, help='Threshold on the localisation score')
parser.add_argument('-aT', '--ampThreshold', type=float, default=1.2, help='Threshold on the amplification ratio')
parser.add_argument('-rS', '--regSample', type=str, default="", help='A pattern for removing controls')
parser.add_argument('-rT', '--regTargets', type=str, help='A pattern for removing targets')
parser.add_argument('-v', '--verbose', type=str, help='A boolean, default ')
args = parser.parse_args()



reads = createReadsMatrix(pathToBam=args.input,bedFile=args.bed,pathToBedtools=args.bedtools,verbose=True)

filteredReads, filteredS, filteredT = filterReads(reads=reads, N=args.minReads, regsamp=args.regSample)

normReads = normalizeReads(filteredReads)

CNVpos, CNVneg = aberrantSamples(filteredReads,conta=args.contaSamples)

final = aberrantAmpliconsFinal(filteredReads,normReads,CNVpos,CNVneg,mode=args.mode)

final.to_csv(output, sep="\t")














