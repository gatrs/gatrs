from __future__ import division
import re
import os
import glob
import sys
import json
import datetime
import time
import threading
from multiprocessing import Pool, TimeoutError, Process
import dendropy
import math
import copy
import random
from tdeftcore import *
from scipy.stats.stats import pearsonr
from scipy.stats import kstest
from scipy.stats import spearmanr
lock=0
bstr=None
forceSeq=0
globalDists={}

def descend(arbor,node,depth):
    if len(node.child_nodes())<1:
        return depth
    node1=node.child_nodes()[0]
    depth=depth+1
    depth1=descend(arbor,node1,depth)
    depth2=0
    if len(node.child_nodes())>1:
        node2=node.child_nodes()[1]
        depth2=descend(arbor,node2,depth)
    depth=max(depth1,depth2)
    return depth

def getDepth(arbor):
        ns=arbor.taxon_namespace
        candidateIx=random.randrange(0,len(ns),1)
        candidateLabel=ns[candidateIx].label
        currentNode=arbor.find_node_with_taxon_label(candidateLabel)
        currentTaxon=currentNode.taxon
        currentLabel=currentTaxon.label
        parentNode=currentNode.parent_node
        while parentNode.parent_node is not None:
            parentNode=parentNode.parent_node
        depth=descend(arbor,parentNode,0)
        return depth

def filterZeroes(inp):
    outp={}
    for key in inp.keys():
        if inp[key]>0:
            outp[key]=inp[key]
    return outp

def lookup(observation,data,ranges,taxon):
    if taxon not in data[observation].keys():
        if taxon.replace("-"," ") in data[observation].keys():
            taxon=taxon.replace("-"," ")
        else:
            return False
    value=data[observation][taxon]
    minval=ranges[observation][0]
    maxval=ranges[observation][1]
    return float(minval),float(value),float(maxval)
        
def scoreByCluster(intree,ranges,data,observation="dlrat"):
    tree=copy.deepcopy(intree)
    ns=tree.taxon_namespace
    differences={}
    tree=copy.deepcopy(intree)
    for taxon in ns:
        if not lookup(observation,data,ranges,taxon.label):
            try:
                tree.prune_taxa_with_labels([taxon.label])
            except:
                pass 
    ns=tree.taxon_namespace
    treeHeight=math.log(len(ns),2)
    tau=int(round(treeHeight/2))
    for taxon in ns:
        neighbors=getLocalNeighborhoodTaxa(tree,taxon.label,tau)
        thisTaxon=lookup(observation,data,ranges,taxon.label)
        if not thisTaxon:
            continue
        difference=0
        counter=1
        taxInf=False
        for tx in neighbors:
            taxInf=lookup(observation,data,ranges,tx)
            if not taxInf:
                continue 
            difference=difference+math.sqrt((taxInf[1]-thisTaxon[1])**2)/(1+taxInf[2])
            counter=counter+1
        if difference==0:
            continue 
        differences[taxon.label]=(difference/counter)
    return(differences)

def compareGeneTree(genea,inPath,minTreeSize=0,morphRanges={},morphData={}):
    tree1=dendropy.Tree.get(path=inPath+genea+".nwk", schema="newick")
    tree2=dendropy.Tree.get(path=inPath+"species.nwk", schema="newick")
    treeHeight=min(getDepth(tree1),getDepth(tree2))
    ns=tree1.taxon_namespace
    ns2=tree1.taxon_namespace
    if len(ns)<minTreeSize or len(ns2)<minTreeSize:
        return {},0
    treeHeight=math.log(len(ns),2)
    tau=int(round(treeHeight/2))
    print(genea)
    score,dist=computeDistanceBetweenTwoTrees(tree1,tree2,tau,1,1,-1,0)
    clusterScores={}
    for cat in morphData.keys():
        cScore=scoreByCluster(tree1,morphRanges,morphData,cat)
        if len(cScore)==0:
            pass
            continue 
        scs=0
        for sct in cScore.keys():
            scs=scs+cScore[sct]
        if len(cScore)/len(ns) < 0.1:
            pass
            continue
        clusterScores[cat]=scs/len(cScore)
    score=filterZeroes(score)
    return [score,dist,clusterScores]

def adjustForTreeEntropy(comps,dists):
    for gene in comps.keys():
        if gene not in dists.keys():
            print("No adjustment found")
            continue
        adjustBy=dists[gene]
        for taxon in comps[gene].keys():
            if (comps[gene][taxon]-adjustBy)**2>0.1:
                comps[gene][taxon]=100*comps[gene][taxon]/adjustBy
            else:
                del(comps[gene][taxon])
    return comps

def checkClustering(comps,clusters,corr,pv,dists):
    threshold = 5 
    clusterdata={}
    for gene in clusters.keys():
        if gene not in clusterdata.keys():
            clusterdata[gene]={}
        for obs in clusters[gene].keys():
            if obs not in clusterdata[gene].keys():
                clusterdata[gene][obs]={}
            if clusters[gene][obs] < threshold:
                clusterdata[gene][obs]["cluster"]=clusters[gene][obs]
                clusterdata[gene][obs]["treeDist"]=dists[gene]
    return clusterdata

def compareAllGenes(inPath=".\\",minTreeSize=0,morphRanges={},morphData={}):
    fileList= glob.glob(inPath+"*.nwk")
    fileCount=len(fileList)
    allComps={}
    allDists={}
    clusterScores={}
    workers=[]
    fnlist=[]
    for i in range(0, fileCount):
        fileName=fileList[i].replace(inPath,"").replace(".nwk","")
        if fileName != "species":
            allComps[fileName],allDists[fileName],clusterScores[fileName]=compareGeneTree(fileName,inPath,minTreeSize,morphRanges,morphData)
    return allComps,clusterScores,allDists

def findLargest(inp):
    mysort={}
    for gene in inp.keys():
        for taxon in inp[gene].keys():
            if inp[gene][taxon] < 1:
                continue
            if inp[gene][taxon] not in mysort.keys():
                mysort[inp[gene][taxon]]=[]
            mysort[inp[gene][taxon]].append([gene,taxon])
    allBuckets=sorted(mysort.keys())
    for size in allBuckets:
        print("--------- "+str(size))
        print(mysort[size])

def correl(a,b):
    correls={}
    pvals={}
    for gene in a.keys():
        lista=[]
        listb=[]
        for morph in b.keys():
          for mtax in b[morph].keys():
            gtax=mtax.replace(" ","-")
            if gtax in a[gene].keys():
                lista.append(a[gene][gtax]**2)
                listb.append(b[morph][mtax]**2)
          if len(lista)<2 or len(listb)<2:
                continue 
          minCoverage=max(5,len(a[gene])*0.1)
          if len(lista)<minCoverage:
              continue
          corr,pv=pearsonr(lista,listb)
          spearman,spv=spearmanr(lista,listb)
          if spv==pv:
              pass
              pv=spv
          if gene not in correls.keys():
            correls[gene]={}
          if gene not in pvals.keys():
            pvals[gene]={}
          correls[gene][morph]=corr       
          pvals[gene][morph]=pv
    return correls,pvals    

def readCSVFileToDictList(fileName):
    import csv
    f=open(fileName)
    source=f.readlines()
    f.close()
    fileData=csv.DictReader(source)
    asList=[]
    for item in fileData:
        asList.append(item)
    return asList

def filterValues(record):
    thresholds={
        "Gill rakers on lower limb":[-1,-1],
        "Gill rakers on upper limb":[-1,-1],
        "Gill rakers total":[4,-1],
        "Anal fin count":[-1,-1],
        "Anal fin spines":[1,-1],
        "Barbels":[2,-1],
        "Dorsal fin count":[-1,-1],
        "Dorsal finlets":[-1,-1],
        "Lateral lines":[-1,-1],
        "Pectoral fin spines":[2,-1],
        "Pelvic fin spines":[2,-1],
        "preanal":[2,-1],
        "Size (cm)":[-1,-1],
        "Ventral finlets":[-1,-1]}
    if record["Observation"] not in thresholds.keys():
        return True 
    if thresholds[record["Observation"]][0]<0 and thresholds[record["Observation"]][1]<0:
        return False 
    if thresholds[record["Observation"]][0]>0:
        if record["Value"]<thresholds[record["Observation"]][0]:
            return False 
    if thresholds[record["Observation"]][1]>0:
        if record["Value"]>thresholds[record["Observation"]][1]:
            return False 
    return True 

def parseMorphObservations(fileName):
    data=readCSVFileToDictList(fileName)
    out={}
    ranges={}
    for record in data:
        if not filterValues(record):
            continue 
        if record["Observation"] not in out.keys():
            out[record["Observation"]]={}
        if record["Observation"] not in ranges.keys():
            ranges[record["Observation"]]=[float(record["Value"]),float(record["Value"])]
        if record["Value"] > ranges[record["Observation"]][1]:
            ranges[record["Observation"]][1]=float(record["Value"])
        if record["Value"] < ranges[record["Observation"]][0]:
            ranges[record["Observation"]][0]=float(record["Value"])
        diff=math.sqrt((float(record["Value"])-float(record["Mean"]))**2)
        out[record["Observation"]][record["Taxon"]]=diff
    return out,ranges

def parseMorphDisplacements(fileName):
    data=readCSVFileToDictList(fileName)
    out={}
    for record in data:
        if record["Observation"] not in out.keys():
            out[record["Observation"]]={}
        out[record["Observation"]][record["Taxon"]]=float(record["Displacement"])
    return out    
        
def get(data,gene,obs,key):
    if data == None:
        return ""
    if gene not in data.keys():
        return ""
    if obs not in data[gene].keys():
        return ""
    if key not in data[gene][obs].keys():
        return ""
    return str(data[gene][obs][key])

morphRanges={}
allMorphDists=parseMorphDisplacements(".\\morphology\\finalDisplacements2.csv")
allMorphComps,morphRanges=parseMorphObservations(".\\morphology\\finalObservations.csv")
minTreeSize=10
allGeneComps,clusters,allDists=compareAllGenes(".\\nwk\\",minTreeSize,morphRanges,allMorphComps)
corr,pv=correl(allGeneComps,allMorphDists)
allpvs=[]
for gene in corr.keys():
    for morph in corr[gene].keys():
        allpvs.append(pv[gene][morph])
sortedpvs=sorted(allpvs)
ct=0
criticals=[]
fdr=0.01
threshold=0.05
tests=len(sortedpvs)
for apv in sortedpvs:
    ct=ct+1
    criticals.append((tests/ct)*fdr)
largestp=0
for i in range(0,len(criticals)):
    if sortedpvs[i]<criticals[i]:
               largestp=sortedpvs[i]
print threshold
threshold=largestp
clusterData=checkClustering(allGeneComps,clusters,corr,pv,allDists)
print("Pvalues < 0.000005 are reported as 0.0")
print("Gene,Observation,Pvalue,Correlation")
outfile=open("g_out.csv","w")
outfile.write("Locus,Observation,pval,corr,cluster,dist,sval,spval\n")
for gene in corr.keys():
    for morph in corr[gene].keys():
        absol=round(math.sqrt(corr[gene][morph]**2),4)
        minCorr=0.1
        if pv[gene][morph] < 0.25:
            if get(clusterData,gene,morph,"cluster") == "":
                sfval=pv[gene][morph]
                spval=str(math.sqrt(sfval))
                sval=str(sfval)
            else:
                sfval=(1+float(get(clusterData,gene,morph,"cluster")))*pv[gene][morph]
                sval=str(sfval)
                spval=str(float(get(clusterData,gene,morph,"cluster"))*math.sqrt(pv[gene][morph]))
            outfile.write(gene+","+morph + "," + str(pv[gene][morph])+"," + str(absol)+","+get(clusterData,gene,morph,"cluster")+","+get(clusterData,gene,morph,"treeDist")+","+sval+","+spval+"\n")
outfile.close()
