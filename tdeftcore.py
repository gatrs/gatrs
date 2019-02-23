import re
import os
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
forceSeq=0
globalDists={}

def descendBranch(rootNode,levels,allNodes):
    if levels < 1:
        return allNodes
    if rootNode.child_nodes() == []:
        allNodes.append(rootNode)
        return allNodes
    for child in rootNode.child_nodes():
        if child.child_nodes() == []:
            allNodes.append(child)
        else:
            allNodes=descendBranch(child,levels-1,allNodes)
            
    return allNodes

def findNeighbors(outData,myNode,levels,ref=0):
  subNodes=[]
  currentNode=myNode
  emptyNode=dendropy.Node
  for i in range(0,levels):
    if currentNode is None:
        continue
    parentNode=currentNode.parent_node
    if parentNode is None:
        continue
    for child in parentNode.child_nodes():
      if child != currentNode:
        if child.child_nodes() != []:
            newNodes=descendBranch(child,levels-i-1,[])
            for eachNode in newNodes:
                subNodes.append(eachNode)
        else:
            subNodes.append(child)
    currentNode=parentNode
  if ref != 0:
      outData[ref]=subNodes
  return subNodes

def getLocalNeighborhood(myTree,taxonString,levels):
    myNode=myTree.find_node_with_taxon_label(taxonString)
    if myNode==None:
        return None
    outData=[]
    return findNeighbors(outData,myNode,levels)

def getNamespace(tree):
    strs=[]
    for item in tree.taxon_namespace:
        strs.append(str(item).replace("'",""))
    return strs

def getLocalNeighborhoodTaxa(myTree,taxonString,levels):
    taxa=[]
    neighborhood=getLocalNeighborhood(myTree,taxonString,levels)
    if neighborhood is not None:
        for node in neighborhood:
            taxa.append(node.taxon.label)
    return taxa

def scoreDifference(taxonInfo,ns1,ns2):
    penalty=0
    if len(taxonInfo)<2:
        return 0
    for neighbor in taxonInfo[0]:
        if neighbor not in taxonInfo[1]:
            if neighbor in ns2:
                penalty=penalty+10
    for neighbor1 in taxonInfo[1]:
        if neighbor1 not in taxonInfo[0]:
            if neighbor1 in ns1:
                penalty=penalty+10
    sampleSize=max(1,len(taxonInfo[0]),len(taxonInfo[1]))
    return int(round(penalty/(sampleSize**2),0))
    
def computeOneSidedDistance(ns1,rawneighborhoods,tree2,levels):
    if type(tree2)==type("abc"):
        tree2=dendropy.Tree.get_from_string(tree2, schema="newick")
    neighborhoods=copy.deepcopy(rawneighborhoods)
    ns2=getNamespace(tree2)
    for taxon2 in ns1:
        if taxon2 not in neighborhoods.keys():
            neighborhoods[taxon2]=[]
        neighborhoods[taxon2].append(getLocalNeighborhoodTaxa(tree2,taxon2,levels))
    distances={}
    totalDistance=0
    for taxon in neighborhoods.keys():
        thisDistance=scoreDifference(neighborhoods[taxon],ns1,ns2)
        distances[taxon]=thisDistance
        totalDistance=totalDistance+thisDistance
    maxSize=len(ns1)
    return totalDistance

def computeDistanceBetweenTwoTrees(tree1,tree2,levels):
    ns1=getNamespace(tree1)
    neighborhoods={}
    for taxon1 in ns1:
        if taxon1 not in neighborhoods.keys():
            neighborhoods[taxon1]=[]
        neighborhoods[taxon1].append(getLocalNeighborhoodTaxa(tree1,taxon1,levels))            
    return computeOneSidedDistance(ns1,neighborhoods,tree2,levels)
 

