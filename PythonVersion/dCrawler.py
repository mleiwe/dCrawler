# -*- coding: utf-8 -*-
"""
dCrawler is an unsupervised clustering algorithm which includes all points within a threshold distance (d) of a centroid. For a detailed description see https://github.com/mleiwe/dCrawler/blob/main/Schema_For_dCrawler.pdf

Inputs
InputMatrix - 
Thd - 

Outputs
ClusterIDs

"""


## Step 0 - Initialisation

# Required libraries
import numpy as np
from scipy.spatial.distance import cdist

# Create functions
def mnl_dCrawler_Crawl(InputMatrix, IntermediateMatrix, Thd, ClusterNum, ClusterIDs):
    #Check any points unallocated?
    nPoints = len(InputMatrix) #The number of points
    NumNaN = np.sum(np.isnan(ClusterIDs)) #Number of unallocated points
    while NumNaN > 0:
        #Initialisation
        i = 0
        while i < nPoints:
            if np.isnan(ClusterIDs[i]):#Is i already part of a cluster? - If it is it will not be a nan
                tPoint = InputMatrix[i,:]
                IntermediateMatrix[i] = np.nan
                
                #Assign to new cluster
                PrelimCluster = Cluster(Point, i, ClustNum)
                
                #Crawl for a cluster starts here
                chk = False
                while
                #Are there any preliminary cluster points outside the distance?
                ClustPoints = PrelimCluster.Points
                ClustCentroid = PrelimCluster.Centroid
                Distances = mnl_MeasureDistances_MatrixToPoint(ClustPoints, ClustCentroid)
                idx = Distances > Thd
                
                if sum(sum(idx)) > 0:
                    temp = ClustPoints * idx
                    idx = temp != np.zeros((1,ndim))
                    PointsToRemove = ClustPoints[idx,:]
                    for p in PointsToRemove:
                        PrelimCluster.RemovePointAndPos(p)
                        
                    PrelimCluster.Update_Cluster()
                
                #Are there any unallocated points within Thd?
                Distances = mnl_MeasureDistances_MatrixToPoint(IntermediateMatrixMatrix, PrelimCluster.Centroid)
                idx = Distances <= Thd
                
                if sum(sum(idx))>0: #If there are any more points to add into the cluster
                    #Identify the closest point and add to cluster
                    minI = np.nanargmin(Distances)
                    NewPoint = IntermediateMatrix[minI]
                    PrelimCluster.AppendPointAndPosition(NewPoint, minI)
                    PrelimCluster.Update_Centroid()
                else:
                    #Assign the Cluster Number to the ClusterIDs
                    idx = PrelimCluster.Pos
                    ClusterIDs[idx] = ClusterNum
                    #Now assign the next lowest numnber to the new ClusterNum  
                    ClusterNum = mnl_FindMinimumMissingInterger(ClusterIDs)                                                     
            else:
                i += 1   
            ## 
            #Insert Progress Bar here at a later date 
            ## 
            NumNaN = np.sum(np.isnan(ClusterIDs)) #Update the number of NaNs remaining   
    return IntermediateMatrix,ClusterIDs

def mnl_MeasureDistances_MatrixToPoint (Matrix,ClusterPos):
    return cdist(Matrix, [ClusterPos], metric='euclidean')

def mnl_FindMinimumMissingInterger(ClusterIDs):
    MxVal = np.max(ClusterIDs)
    MnVal = MxVal
    while MnVal == MxVal:
        for Val in range(MxVal):
            idx = ClusterIDs == Val
            if sum(idx) < 0:
                MnVal = Val
        MnVal = MxVal + 1
    return MnVal
# Create classes
class Cluster():
    def __init__(self, Point, ListPos, ClustNum):
        self.Points = [Point]
        self.Centroid = Point
        self.Pos = ListPos
        self.Number = ClustNum
    
    def Update_Centroid(self):
        self.Centroid = np.mean(self.Points, axis=0)
    
    def AppendPointAndPosition(self, NewPoint, NewPos):
        self.Points = np.append(self.Points, [NewPoint], axis=0)
        self.Pos = np.append(self.Pos, NewPos)
    
    def RemovePointAndPos(self, OldPoint):
        #Mask of points to keep
        a = self.Points != OldPoint.all(axis=1)
        #Now only keep the ones that belong to the cluster
        self.Points = self.Points[a]
        self.Pos = self.Pos[a]
    
    def ConvertToDict(self):
        
        

# Create variables
nPoints = len(InputMatrix) #The number of points
nDim = len(InputMatrix[0]) #The number of dimensions
ClusterNum = 0 #The initial number to start the clustering from
IntermediateMatrix = InputMatrix #IntermediateMatrix to be filled with NaNs as points are allocated
ClusterIDs = np.empty(nPoints) * np.nan

def dCrawler(InputMatrix,Thd)
## Step 1 - Initial Crawl
IntermediateMatrix, ClusterIDs = mnl_dCrawler_Crawl(InputMatrix, IntermediateMatrix, Thd, ClusterNum, ClusterIDs)

## Step 2 - Adjust Phase

## Step 3 - Merge Phase
