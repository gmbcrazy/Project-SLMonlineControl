#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import suite2p
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
import glob
##python -m pip install -U scikit-image
from skimage import measure
from natsort import natsorted

def BinList_PSTHHeatMap(MultiBinFileList,BaselineInd,ResponseInd,ops0):

   # print(BaselineInd)
   # print(ResponseInd)
    nplanes=ops0['nplanes']
    #print(nplanes)
    baseMap=np.zeros((len(MultiBinFileList),nplanes,ops0['Ly'],ops0['Lx']))
    ResponseMap=np.zeros((len(MultiBinFileList),nplanes,ops0['Ly'],ops0['Lx']))
    for TrialI,Trial in enumerate(MultiBinFileList):
        #print(Trial)
        rawBin = suite2p.io.BinaryFile(Ly=ops0['Ly'],Lx=ops0['Lx'], filename=Trial)
        #print(rawBin.shape)
        FramePerPlane=np.floor(rawBin.shape[0]/ops0['nplanes'])
        TotalFrameNeed=np.int32(np.floor(FramePerPlane)*nplanes)
        #print(TotalFrameNeed)
    #print(rawBin.shape)
        nplanes=ops0['nplanes']
    
        for plane_id in range(0,nplanes):
            #print(plane_id)
                #print(range(0+plane_id,TotalFrameNeed,nplanes))
            #print(np.array(range(0+plane_id,TotalFrameNeed,nplanes)))
            plane_data = rawBin[range(0+plane_id,TotalFrameNeed,nplanes),:,::-1]
        #print(plane_data.shape)
            plane_data=np.double(plane_data)

            baseMap[TrialI,plane_id,:,:]=np.mean(plane_data[BaselineInd,:,:],axis=0)
            ResponseMap[TrialI,plane_id,:,:]=np.mean(plane_data[ResponseInd,:,:],axis=0)

    #DimEnd=len(baseMap.shape)
    baseMap=np.mean(baseMap,axis=0)
    ResponseMap=np.mean(ResponseMap,axis=0)
    ResponseMap=ResponseMap-baseMap
    ResponseMap=np.transpose(ResponseMap,(0,2,1))
    return ResponseMap


def plotCellCenter(ax, Ly, cellCenter, Radius, colorCell, LineWidth):

    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (cellCenter.shape[0], 1))

    if np.shape(Radius) == 1:
        Radius = np.tile(Radius, cellCenter.shape[0])
    
   # fig, ax = plt.subplots()
    if len(cellCenter.shape)==1:
       cellN=1
       CoorDim=len(cellCenter)
    else:
       CoorDim=cellCenter.shape[1]

    #print(cellN)
    # Loop through each center point
    if cellN==1:

       x = cellCenter[1]
       y = Ly-cellCenter[0]

       print(cellCenter)
        # Plot a circle around the center with the specified radius, color, and line width
       circle = Circle((x, y), Radius, color=colorCell[0], linewidth=LineWidth, fill=False)
       ax.add_patch(circle)
       plt.show()
    else:
        for i in range(cellCenter.shape[0]):
            # Extract x and y coordinates of the current center
            x = cellCenter[i, 1]
            y = cellCenter[i, 0]

        # Plot a circle around the center with the specified radius, color, and line width
            circle = Circle((x, y), Radius[i], color=colorCell[i], linewidth=LineWidth, fill=False)
            ax.add_patch(circle)
            plt.show()





# Example usage:
# cellCenter = np.array([[x1, y1], [x2, y2], ...])
# Radius = np.array([r1, r2, ...])
# colorCell = np.array([[r1, g1, b1], [r2, g2, b2], ...]) or a single color
# plotCellCenter(cellCenter, Radius, colorCell)
def plotCellCenter3D(ax,cellCenter, Radius, colorCell, LineWidth):
    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (cellCenter.shape[0], 1))

    if np.shape(Radius) == 1:
        Radius = np.tile(Radius, cellCenter.shape[0])
    
    #print(Radius.shape)
    plt.show()
    if len(cellCenter.shape)==1:
       cellN=1
       CoorDim=len(cellCenter)
    else:
       CoorDim=cellCenter.shape[1]
    
    if CoorDim == 2 and cellN > 1:
        cellCenter = np.hstack((cellCenter, np.zeros((len(cellCenter), 1))))
    elif CoorDim == 2 and cellN == 1:
        cellCenter=np.append(cellCenter,0)
    else:
        cellCenter=cellCenter

    if cellN==1:
       theta = np.linspace(0, 2*np.pi, 50)
       x = cellCenter[1] + Radius * np.cos(theta)
       y = cellCenter[0] + Radius * np.sin(theta)
       z = np.zeros_like(x) + cellCenter[2]
       ax.plot(x, y, z, color=colorCell[0], linewidth=LineWidth, linestyle='-')
       plt.show()
    else:
    
        for i in range(cellN):
            theta = np.linspace(0, 2*np.pi, 50)
            x = cellCenter[i, 1] + Radius[i] * np.cos(theta)
            y = cellCenter[i, 0] + Radius[i] * np.sin(theta)
            z = np.zeros_like(x) + cellCenter[i, 2]
            ax.plot(x, y, z, color=colorCell[i], linewidth=LineWidth, linestyle='-')
            plt.show()
    
   # ax.set_xlabel('X-axis')
   # ax.set_ylabel('Y-axis')
    ##ax.set_zlabel('Z-axis')
   # ax.set_title('Cell Centers with Circles (3D)')


# Iterate over X and Y values and search for matching files
#PointFile=[None]*len(PointI_values)
def PSTHFromBin(matching_files,statCell,plane_idx,SinnL,TrialNum):
    PSTH=np.zeros((SingL,TrialNum))
    for Ind,Trial in enumerate(matching_files):
        rawBin = suite2p.io.BinaryFile(Ly=ops0['Ly'],Lx=ops0['Lx'], filename=Trial)
        print(rawBin.shape)
        print(Trial)
        plane_data = rawBin[range(0+plane_idx,SingL*nplanes,nplanes),:,::-1]
    
        stat_after_extraction, F, Fneu, F_chan2, Fneu_chan2 = suite2p.extraction_wrapper(statCell, f_reg=plane_data,
                                                          f_reg_chan2 = plane_data,ops=ops0)
        dF = F.copy() - ops1['neucoeff']*Fneu

        dF = suite2p.extraction.preprocess(
        F=dF,
        baseline=ops1['baseline'],
        win_baseline=ops1['win_baseline'],
        sig_baseline=ops1['sig_baseline'],
        fs=ops1['fs'],
        prctile_baseline=ops1['prctile_baseline']
        )
        PSTH[:,Ind]=np.double(F)

    return PSTH
    #preStim=range(0,10)
   # print(PSTH.shape)
    
#baseLine=np.mean(np.mean(PSTH[preStim,:],axis=1),axis=0)
#Response=np.mean(PSTH,axis=1)
#print(Response.shape)
#baseLine=np.tile(baseLine.T,(SingL,1)).T
#print(baseLine.shape)
#plt.imshow(np.mean(PSTH,axis=1)-baseLine,cmap='seismic')
#np.shape(Response-baseLine)


#fig, ax = plt.subplots(1,3)
#ax[0].imshow(baseMap[:,:,0],cmap='seismic')


# Iterate over X and Y values and search for matching files
#PointFile=[None]*len(PointI_values)
def PSTHFromBin(matching_files,statCell,plane_idx,SinnL,TrialNum):
    PSTH=np.zeros((SingL,TrialNum))
    for Ind,Trial in enumerate(matching_files):
        rawBin = suite2p.io.BinaryFile(Ly=ops0['Ly'],Lx=ops0['Lx'], filename=Trial)
        print(rawBin.shape)
        print(Trial)
        plane_data = rawBin[range(0+plane_idx,SingL*nplanes,nplanes),:,::-1]
    
        stat_after_extraction, F, Fneu, F_chan2, Fneu_chan2 = suite2p.extraction_wrapper(statCell, f_reg=plane_data,
                                                          f_reg_chan2 = plane_data,ops=ops0)
        dF = F.copy() - ops1['neucoeff']*Fneu

        dF = suite2p.extraction.preprocess(
        F=dF,
        baseline=ops1['baseline'],
        win_baseline=ops1['win_baseline'],
        sig_baseline=ops1['sig_baseline'],
        fs=ops1['fs'],
        prctile_baseline=ops1['prctile_baseline']
        )
        PSTH[:,Ind]=np.double(F)

    return PSTH
    #preStim=range(0,10)
   # print(PSTH.shape)
    
#baseLine=np.mean(np.mean(PSTH[preStim,:],axis=1),axis=0)
#Response=np.mean(PSTH,axis=1)
#print(Response.shape)
#baseLine=np.tile(baseLine.T,(SingL,1)).T
#print(baseLine.shape)
#plt.imshow(np.mean(PSTH,axis=1)-baseLine,cmap='seismic')
#np.shape(Response-baseLine)


#fig, ax = plt.subplots(1,3)
#ax[0].imshow(baseMap[:,:,0],cmap='seismic')



def plotCellCenter3DGroup(cellCenterGroup, Radius, *args):
    if len(args) < 1:
        colorCell = plt.cm.jet(range(len(cellCenterGroup)))
        LineWidth = 1
    elif len(args) == 1:
        colorCell = args[0]
        LineWidth = 1
    elif len(args) == 2:
        colorCell = args[0]
        LineWidth = args[1]
    else:
        # Handle other cases if needed
        pass

    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (len(cellCenterGroup), 1))

    #print(Radius.shape)
    if Radius.shape[0] == 1:
        Radius = np.tile(Radius, len(cellCenterGroup))
        
    #print(Radius.shape)
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    for iGroup in range(len(cellCenterGroup)):
        if cellCenterGroup[iGroup] is not None:
            plotCellCenter3D(cellCenterGroup[iGroup], Radius[iGroup], colorCell[iGroup], LineWidth)

    plt.show()
    return ax

def plotCellCenter3D(ax,cellCenter, Radius, colorCell, LineWidth):
    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (cellCenter.shape[0], 1))

    if np.shape(Radius) == 1:
        Radius = np.tile(Radius, cellCenter.shape[0])
    
    #print(Radius.shape)
    plt.show()
    if len(cellCenter.shape)==1:
       cellN=1
       CoorDim=len(cellCenter)
    else:
       CoorDim=cellCenter.shape[1]
    
    if CoorDim == 2:
        cellCenter = np.hstack((cellCenter, np.zeros((len(cellCenter), 1))))

    for i in range(cellN):
        theta = np.linspace(0, 2*np.pi, 50)
        x = cellCenter[i, 1] + Radius[i] * np.cos(theta)
        y = cellCenter[i, 0] + Radius[i] * np.sin(theta)
        z = np.zeros_like(x) + cellCenter[i, 2]
        ax.plot(x, y, z, color=colorCell[i], linewidth=LineWidth, linestyle='-')
        plt.show()
    
   # ax.set_xlabel('X-axis')
   # ax.set_ylabel('Y-axis')
    ##ax.set_zlabel('Z-axis')
   # ax.set_title('Cell Centers with Circles (3D)')


def Suite2pCellIDMap(ops,stat,iscell):
    #"""
    #Get cellIDMap from Suite2p processed data.
    #"""


    cellValid = iscell[:, 0]
    stat=stat[iscell[:, 0]==1]
    #cellIDMap = np.zeros_like(ops['meanImg'])
    cellIDMap = np.zeros((ops['Ly'],ops['Lx']))
    MedCenter = []

    m, n = cellIDMap.shape

    validCellList = np.where(cellValid)[0]

    cellID=0
    for statCell in stat:
        cellID = cellID + 1
        roiPix = np.ravel_multi_index((statCell['ypix'], statCell['xpix']), cellIDMap.shape)
        cellIDMap.flat[roiPix] = cellID

        MedCenter.append([max(min(int(np.median(statCell['ypix'])), m), 1),
                          max(min(int(np.median(statCell['xpix'])), n), 1)])

    MedCenter = np.array(MedCenter)

    temp = np.setdiff1d(cellIDMap, 0)
    CellPixCount = np.histogram(temp, bins=np.arange(temp.min(), temp.max() + 2))[0]

    return cellIDMap, CellPixCount, MedCenter


def plotCellBoundary(cellIDMap, cellIDmark, *args):
    # Check the number of input arguments
    if len(args) < 1:
        colorCell = plt.cm.jet(range(len(cellIDmark)))
        LineWidth = 3
    elif len(args) == 1:
        colorCell = args[0]
        LineWidth = 3
    elif len(args) == 2:
        colorCell = args[0]
        LineWidth = args[1]
    else:
        # Handle other cases if needed
        pass

    # If colorCell is a single row, replicate it for each cellIDmark
    if colorCell.shape[0] == 1:
        colorCell = [colorCell] * len(cellIDmark)

    # Get cell boundaries using the CellIDMap2Boundary function
    cellBoundary = CellIDMap2Boundary(cellIDMap, cellIDmark)

    # Loop through each cell ID
    for i in range(len(cellIDmark)):
        # Loop through each boundary of the current cell ID
        for iB in range(len(cellBoundary[i])):
            # Plot the cell boundary using specified color and LineWidth
            plt.plot(cellBoundary[i][iB][:, 1], cellBoundary[i][iB][:, 0],
                     color=colorCell[i], linewidth=LineWidth)

    plt.show()  # Show the plot


def CellIDMap2Boundary(cellIDMap, cellCandidate=None):
    if cellCandidate is None:
        # Identify unique non-zero cell IDs
        cellCandidate = np.unique(cellIDMap)
        cellCandidate = cellCandidate[cellCandidate != 0]
        cellCandidate = np.sort(cellCandidate)

    cellBoundary = []

    # Loop through each cell candidate
    for icell in cellCandidate:
        # Create a temporary binary map for the current cell
        tempMap = np.zeros_like(cellIDMap)
        tempMap[cellIDMap == icell] = 1

        # Find the boundary of the cell using bwboundaries
        boundaries = measure.find_contours(tempMap, 0.5)

        cellBoundary.append(boundaries)

    return cellBoundary


def MultiMatrix3DPlotZ(ax,Data, ZPlot,vLim,Alpha):
    """
    Plot 3D matrices along the Z dimension.
    """
    X, Y = np.meshgrid(np.arange(Data.shape[1]), np.arange(Data.shape[0]))
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    vmax=vLim[1]
    vmin=vLim[0]
    for k in range(Data.shape[2]):
        Z = np.ones_like(Data[:, :, k]) * ZPlot[k]  # Specify Z-values for each surface
        temp=Data[:, :, k]
        ax.plot_surface(X, Y, Z, facecolors=plt.cm.seismic((Data[:, :, k]-vmin)/(vmax-vmin)), alpha=Alpha,shade=False)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()
    return ax

def plotCellCenter3DGroup(cellCenterGroup, Radius, *args):
    if len(args) < 1:
        colorCell = plt.cm.jet(range(len(cellCenterGroup)))
        LineWidth = 1
    elif len(args) == 1:
        colorCell = args[0]
        LineWidth = 1
    elif len(args) == 2:
        colorCell = args[0]
        LineWidth = args[1]
    else:
        # Handle other cases if needed
        pass

    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (len(cellCenterGroup), 1))

    #print(Radius.shape)
    if Radius.shape[0] == 1:
        Radius = np.tile(Radius, len(cellCenterGroup))
        
    #print(Radius.shape)
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    for iGroup in range(len(cellCenterGroup)):
        if cellCenterGroup[iGroup] is not None:
            plotCellCenter3D(cellCenterGroup[iGroup], Radius[iGroup], colorCell[iGroup], LineWidth)

    plt.show()
    return ax

def plotCellCenter3D(ax,cellCenter, Radius, colorCell, LineWidth):
    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (cellCenter.shape[0], 1))

    if np.shape(Radius) == 1:
        Radius = np.tile(Radius, cellCenter.shape[0])
    
    #print(Radius.shape)
    plt.show()
    if len(cellCenter.shape)==1:
       cellN=1
       CoorDim=len(cellCenter)
    else:
       CoorDim=cellCenter.shape[1]
    
    if CoorDim == 2:
        cellCenter = np.hstack((cellCenter, np.zeros((len(cellCenter), 1))))

    for i in range(cellN):
        theta = np.linspace(0, 2*np.pi, 50)
        x = cellCenter[i, 1] + Radius[i] * np.cos(theta)
        y = cellCenter[i, 0] + Radius[i] * np.sin(theta)
        z = np.zeros_like(x) + cellCenter[i, 2]
        ax.plot(x, y, z, color=colorCell[i], linewidth=LineWidth, linestyle='-')
        plt.show()
    
   # ax.set_xlabel('X-axis')
   # ax.set_ylabel('Y-axis')
    ##ax.set_zlabel('Z-axis')
   # ax.set_title('Cell Centers with Circles (3D)')


