#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import suite2p
from suite2p.registration import register
from suite2p.registration import rigid
from suite2p.io import tiff
from skimage.registration import phase_cross_correlation
import numpy as np
import scipy.io
from scipy.signal import convolve, convolve2d
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
import glob
##python -m pip install -U scikit-image
from skimage import measure
from natsort import natsorted
from collections import defaultdict
import re
import time
import FastBin_Suite2p as FBS
# Iterate over X and Y values and search for matching files
#PointFile=[None]*len(PointI_values)
def NeuroFromBin(matching_files,statCell,plane_idx,SingL,ops0):
    
    TrialNum=len(matching_files)
    Fsig=np.zeros((len(statCell),SingL,TrialNum))
    DeltaF=np.zeros((len(statCell),SingL,TrialNum))  
    spks=np.zeros((len(statCell),SingL,TrialNum))
    nplanes=ops0['nplanes']
    Invalid=[]
    maxIneed=np.max(SingL)

    for Ind,Trial in enumerate(matching_files):
        rawBin = suite2p.io.BinaryFile(Ly=ops0['Ly'],Lx=ops0['Lx'], filename=Trial)
        if maxIneed*nplanes>rawBin.shape[0]:
           Invalid.append(Ind)
           print(Trial+' is too short')
           continue
        plane_data = rawBin[range(0+plane_idx,SingL*nplanes,nplanes),:,:]
        stat_after_extraction, F, Fneu, F_chan2, Fneu_chan2 = suite2p.extraction_wrapper(statCell, f_reg=plane_data,f_reg_chan2 = plane_data,ops=ops0)
        dF = F.copy() - ops0['neucoeff']*Fneu
        dF = suite2p.extraction.preprocess(F=dF, baseline=ops0['baseline'],win_baseline=ops0['win_baseline'],sig_baseline=ops0['sig_baseline'],fs=ops0['fs'],prctile_baseline=ops0['prctile_baseline'])
        spksTemp = suite2p.extraction.oasis(F=dF, batch_size=ops0['batch_size'], tau=ops0['tau'], fs=ops0['fs'])
        Fsig[:,:,Ind]=np.double(F)
        DeltaF[:,:,Ind]=np.double(dF)
        spks[:,:,Ind]=np.double(spksTemp)
    #print(Ind)
    Fsig=np.delete(Fsig,Invalid,axis = 2)
    DeltaF=np.delete(DeltaF,Invalid, axis = 2)
    spks=np.delete(spks,Invalid, axis = 2)
    return Fsig, DeltaF, spks
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

def BinList_PSTHHeatMap(MultiBinFileList,BaselineInd,ResponseInd,ops0):


    nplanes=ops0['nplanes']
    baseMap=np.zeros((len(MultiBinFileList),nplanes,ops0['Ly'],ops0['Lx']))
    ResponseMap=np.zeros((len(MultiBinFileList),nplanes,ops0['Ly'],ops0['Lx']))
    Invalid=[]
    for TrialI,Trial in enumerate(MultiBinFileList):
        rawBin = suite2p.io.BinaryFile(Ly=ops0['Ly'],Lx=ops0['Lx'], filename=Trial)
        maxIneed=np.max([np.max(BaselineInd),np.max(ResponseInd)])
        if maxIneed*nplanes>rawBin.shape[0]:
           Invalid.append(TrialI)
           print(Trial+' is too short')
           continue
      
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

    baseMap=np.delete(baseMap,Invalid,axis = 0)
    ResponseMap=np.delete(ResponseMap,Invalid, axis = 0)
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
       #plt.show()
    else:
        for i in range(cellCenter.shape[0]):
            # Extract x and y coordinates of the current center
            x = cellCenter[i, 1]
            y = cellCenter[i, 0]

        # Plot a circle around the center with the specified radius, color, and line width
            circle = Circle((x, y), Radius[i], color=colorCell[i], linewidth=LineWidth, fill=False)
            ax.add_patch(circle)
            #plt.show()





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






def PointLaser_files(file_names):
    # Regular expression to extract the point and laser level from the file name
    pattern = re.compile(r"Laser(\d+\.\d+)[a-zA-Z]Point(\d+)\.bin")

    # Nested dictionary to store the results
    # Defaultdict with a lambda to create another defaultdict for nested structure
    organized_files = defaultdict(lambda: defaultdict(list))

    

    for file_name in file_names:
        match = pattern.search(file_name)
        if match:
            # Extract laser level and point number
            laser_level = match.group(1)
            point_number = match.group(2)

            # Organize file into the dictionary
            organized_files[point_number][laser_level].append(file_name)

    # Convert inner defaultdicts to regular dictionaries for easier use
    OrganizedFile = {point: dict(lasers) for point, lasers in organized_files.items()}
    return OrganizedFile
    
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


def get_new_files(previous_files, current_files):
    # Return files that are in current_files but not in previous_files
    if previous_files is None:
        return [file for file in currentfiles]
    else:
        return [file for file in current_files if file not in previous_files]

def ImgNormalize(Image1):
    mean = np.mean(Image1)
    std_dev = np.std(Image1)
    Image2 = (Image1 - mean) / std_dev
    return Image2

def monitor_folderBinFiles(folder_path, reference_Data, plane_idx, ops0, TimeTh=15, ExistingBinFile=None, folder_keywords=None):
    # Initialize previous files list with ExistingBinFile if provided, else empty list
    if ExistingBinFile is None:
        previous_files = []
    elif ExistingBinFile=='Current':
        previous_files = glob.glob(folder_path + '/*.bin')
    else:
        previous_files = ExistingBinFile
    
    # Initialize a time marker for the timing threshold
    time_marker = time.time()
    Pixel_ShiftAll=[]
    corrValue=[]
    fileList=[]
    while True:
          pattern = folder_path + '/*'
          if folder_keywords:
             pattern += folder_keywords
          pattern += '*.bin'
          current_files = glob.glob(pattern)
          new_files = get_new_files(previous_files, current_files)
        #print(new_files)
          if new_files:
            # Reset the time marker since new files have been found
             print(new_files)             
            # Process new files
             for file in new_files:       
                 time.sleep(10)
                # Construct the full path to the file
                #file_path = os.path.join(folder_path, file)
                #print(file_path)
                 TestBin, FramePerPlane, TotalFrameNeed = FBS.LoadBin(file, ops0)
                 Test_data = TestBin[range(0 + plane_idx, TestBin.shape[0], ops0['nplanes']), :, :]
                 TestPlane = np.mean(Test_data, axis=0)
                
                # Calculate pixel shift
                 Img1=ImgNormalize(reference_Data)
                 Img2=ImgNormalize(TestPlane)
                 pixel_shift, junk1, junk2 = phase_cross_correlation(Img1,Img2)
                
                 image_product = np.fft.fft2(Img1) * np.fft.fft2(Img2).conj()
                 cc_image = np.fft.fftshift(np.fft.ifft2(image_product))
 
                 Pixel_ShiftAll.append(pixel_shift)
                 corrValue.append(np.max(np.max(cc_image.real))/ops0["Lx"]/ops0["Ly"])
                 fileList.append(file)
                # Check if pixel shift exceeds threshold
                #if np.max(pixel_shift) > threshold:
                 print(f"{file} Pixel shift detected! Offset: {pixel_shift}")
                 time_marker = time.time() 
                 
                #return  # Or take other appropriate action
                
          else:
            # If no new files are found, check the time since the last new file
              print(f"No new files detected {np.int16(TimeTh-time.time()+time_marker)} seconds left")
              if time.time() - time_marker > TimeTh:
                 print("Stopping monitoring.")
                 break
        
        # Update previous files list
          previous_files = current_files
        
        # Check if the 'S' key is pressed to manually stop monitoring
        #if keyboard.is_pressed('s'):
            #print("Monitoring stopped by user.")
           # break
        
        # Wait for a while before checking again
          time.sleep(5)
    
    nFile=len(fileList)
    Pixel_ShiftAll=np.array(Pixel_ShiftAll)
    corrValue=np.array(corrValue)
        #if not nFile > 0:

    return Pixel_ShiftAll,corrValue,fileList




def smooth(data, smooth):
    if not isinstance(data, np.ndarray):
        raise ValueError("Data must be a numpy array.")
    if not isinstance(smooth, (list, tuple, np.ndarray)) or not all(isinstance(x, (int, float)) for x in smooth):
        raise ValueError("Smooth must be a list or tuple of integers or floats.")

    if len(smooth) == 0 or all(s == 0 for s in smooth):
        return data  # No smoothing required

    if data.ndim not in [1, 2]:
        raise ValueError("Smoothing applies only to vectors or matrices.")
    if data.ndim == 1:
        # Vector smoothing
        if data.shape[0] == 1:
            data = data.T
        size = len(data)
        size = min(size, 1001)
        r = np.linspace(-1, 1, 2 * size + 1)
        smooth_factor = smooth[0] / size
        smoother = np.exp(-0.5 * (r / smooth_factor) ** 2)
        smoother /= np.sum(smoother)
        smoothed = np.convolve(data, smoother, mode='same')
        return smoothed
    else:
        # Matrix smoothing
        if len(smooth) == 1:
            smooth = [smooth[0], smooth[0]]  # Same smoothing for both dimensions

        vertical_size = data.shape[0]
        horizontal_size = data.shape[1]
        vertical_size = min(vertical_size, 1001)
        horizontal_size = min(horizontal_size, 1001)

        if smooth[0] > 0:
            r = np.linspace(-1, 1, 2 * vertical_size + 1)
            v_smooth_factor = smooth[0] / vertical_size
            v_smoother = np.exp(-0.5 * (r / v_smooth_factor) ** 2)
            v_smoother /= np.sum(v_smoother)
        if smooth[1] > 0:
            r = np.linspace(-1, 1, 2 * horizontal_size + 1)
            h_smooth_factor = smooth[1] / horizontal_size
            h_smoother = np.exp(-0.5 * (r / h_smooth_factor) ** 2)
            h_smoother /= np.sum(h_smoother)

        if smooth[0] == 0:
            smoothed = convolve(data, h_smoother[None, :], mode='same')
        elif smooth[1] == 0:
            smoothed = convolve(data, v_smoother[:, None], mode='same')
        else:
            smoothed = convolve2d(data, np.outer(v_smoother, h_smoother), mode='same')

        return smoothed
