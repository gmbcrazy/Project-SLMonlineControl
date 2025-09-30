#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import suite2p
import cv2
from suite2p.registration import register
from suite2p.registration import rigid
from suite2p.io import tiff
from skimage.registration import phase_cross_correlation
from scipy.ndimage import gaussian_filter
import numpy as np
import scipy.io
from scipy.signal import convolve, convolve2d
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import glob
##python -m pip install -U scikit-image
from skimage import measure
from collections import defaultdict
import re
import time
import FastBin_Suite2p as FBS
import os
import fnmatch
# Iterate over X and Y values and search for matching files
#PointFile=[None]*len(PointI_values)


## This is used for getting neural signals from 1 bin file, including multiple Zseries/SLM trials 
def NeuroFromBin(matching_file,statCell,Cellplane,SingL,ops1,DelFrameRep,refMeanImg=None):
    """
    Extract neural signals from a bin file, including multiple Z-series or SLM trials.

    Parameters:
    - matching_file: The file path for the bin file containing imaging data.
    - statCell: A list of dictionaries containing ROI information for each cell.
    - Cellplane: An array indicating which plane each cell belongs to.
    - SingL: The total number of imaging frames per trial.
    - ops0: A dictionary containing various options and parameters for Suite2p extraction.
    - DelFrameRep: A list of frame indices to be deleted from the imaging sequence.
    - refMeanImg: (optional) A reference image for alignment. If provided, pixel shifts are calculated.

    Returns:
    - Fsig: The fluorescence signals for each cell.
    - DeltaF: The delta F signals for each cell.
    - spks: The inferred spike signals for each cell.
    """
        # Calculate the effective number of frames after removing repetitions
        
    DataL=SingL-len(DelFrameRep)
    Fsig=np.zeros((len(statCell),DataL))
    DeltaF=np.zeros((len(statCell),DataL))  
    spks=np.zeros((len(statCell),DataL))
    nplanes=ops1['nplanes']
    

    rawBin, FramePerPlane, TotalFrameNeed=FBS.LoadBin(matching_file,ops1)
    PlaneAll=np.unique(Cellplane)

    for planeID in PlaneAll:
         # Extract data for the current plane
        plane_data = rawBin[range(planeID,SingL*nplanes,nplanes),:,:]
        plane_data = np.delete(plane_data, DelFrameRep, axis=0)
    
            # If a reference image is provided, calculate the pixel shift
        yxShift=np.zeros((3,2))
        if refMeanImg is not None:
           current_mean_img = np.mean(plane_data, axis=0)
           yxShift[planeID,:] = Pixel_Shift(refMeanImg[planeID, :, :], current_mean_img, 1)
           if not np.all(yxShift == 0):
             # Apply reverse pixel shift to the plane data
              plane_data = apply_affine_transform_to_all_frames(plane_data, -yxShift[planeID,:])


        Ind=np.where(Cellplane == planeID)[0]
        statCelltemp=statCell[Ind]
        stat_after_extraction, F, Fneu, F_chan2, Fneu_chan2 = suite2p.extraction_wrapper(statCelltemp, f_reg=plane_data,f_reg_chan2 = plane_data,ops=ops1)
        dF = F.copy() - ops1['neucoeff'] * Fneu
        dF = suite2p.extraction.preprocess(
            F=dF,
            baseline=ops1['baseline'],
            win_baseline=ops1['win_baseline'],
            sig_baseline=ops1['sig_baseline'],
            fs=ops1['fs'],
            prctile_baseline=ops1['prctile_baseline']
        )
        spksTemp = suite2p.extraction.oasis(F=dF, batch_size=ops1['batch_size'], tau=ops1['tau'], fs=ops1['fs'])
        Fsig[Ind,:]=np.double(F)
        DeltaF[Ind,:]=np.double(dF)
        spks[Ind,:]=np.double(spksTemp)

    return Fsig, DeltaF, spks,yxShift
   #preStim=range(0,10)
   # print(PSTH.shape)
    
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

def plotCellCenter3D(ax, cellCenter, Radius, colorCell, LineWidth):
    if len(colorCell) == 1:
        colorCell = np.tile(colorCell, (cellCenter.shape[0], 1))

    if np.shape(Radius) == 1:
        Radius = np.tile(Radius, cellCenter.shape[0])

    # Ensure cellN is defined before use
    if len(cellCenter.shape) == 1:
        cellN = 1
        CoorDim = len(cellCenter)
    else:
        cellN = cellCenter.shape[0]  # Set cellN properly
        CoorDim = cellCenter.shape[1]

    if CoorDim == 2 and cellN > 1:
        cellCenter = np.hstack((cellCenter, np.zeros((len(cellCenter), 1))))
    elif CoorDim == 2 and cellN == 1:
        cellCenter = np.append(cellCenter, 0)
    else:
        cellCenter = cellCenter

    if cellN == 1:
        theta = np.linspace(0, 2 * np.pi, 50)
        x = cellCenter[1] + Radius * np.cos(theta)
        y = cellCenter[0] + Radius * np.sin(theta)
        z = np.zeros_like(x) + cellCenter[2]
        ax.plot(x, y, z, color=colorCell[0], linewidth=LineWidth, linestyle='-')
    else:
        for i in range(cellN):
            theta = np.linspace(0, 2 * np.pi, 50)
            x = cellCenter[i, 1] + Radius[i] * np.cos(theta)
            y = cellCenter[i, 0] + Radius[i] * np.sin(theta)
            z = np.zeros_like(x) + cellCenter[i, 2]
            ax.plot(x, y, z, color=colorCell[i], linewidth=LineWidth, linestyle='-')

import numpy as np
import matplotlib.pyplot as plt

def plotCellCenter(ax, cellCenter, Radius, colorCell, LineWidth):
    """
    Plots cell centers as circles in 2D.
    
    Parameters:
        ax : matplotlib axis
            The axis on which to plot.
        cellCenter : numpy.ndarray
            Array of shape (N,2) with (x,y) positions of cell centers.
        Radius : float or numpy.ndarray
            Radius of circles. Can be scalar or array of shape (N,).
        colorCell : numpy.ndarray
            Color information for each circle.
        LineWidth : float
            Line width of the circles.
    """
    if len(colorCell.shape) == 1:
        colorCell = np.tile(colorCell, (cellCenter.shape[0], 1))

    if np.isscalar(Radius):
        Radius = np.full(cellCenter.shape[0], Radius)

    for i in range(cellCenter.shape[0]):
        circle = plt.Circle((cellCenter[i, 1], cellCenter[i, 0]), Radius[i], 
                            edgecolor=colorCell[i], facecolor='none', linewidth=LineWidth)
        ax.add_patch(circle)

def multi_planes_2d_show(img, cell_boundary, pos3d, pos3d_label, zdepth, color_cell, img_clim, plot_param=None):
    """
    Displays 2D images of multi-plane data and overlays cell boundaries and positions.
    
    Parameters:
        img (numpy.ndarray): Image data, can be 2D or 3D (multi-plane).
        cell_boundary (list): Cell boundary data to overlay on the image.
        pos3d (numpy.ndarray): 3D positions of cells (X, Y, Z).
        pos3d_label (list): Labels for cell positions.
        zdepth (list): Z-depth values corresponding to each plane.
        color_cell (numpy.ndarray or list): Colors to use for cell boundaries.
        img_clim (tuple): Limits for image display (contrast adjustment).
        plot_param (dict, optional): Plotting parameters. Defaults will be used if not provided.
    
    Returns:
        fig, axes: Matplotlib figure and axes handles for further customization.
    """
    dim = img.shape
    dim_pos = pos3d.shape
    
    if plot_param is None:
        plot_param = {
            'RowPlot': True,
            'RowColNum': 1,
            'RowColID': 1,
            'EdgeParam': [0.06, 0.1, 0.06, 0.06, 0.06, 0.06],
            'CellCenterWidth': 1,
            'CellBoundaryWidth': 0.5,
            'PlotCenter': True
        }
    
    # Ensure color_cell is a numpy array
    color_cell = np.array(color_cell)
    if color_cell.ndim == 1:
        color_cell = np.tile(color_cell, (dim_pos[0], 1))
    
    pos_z = pos3d[:, 2] if dim_pos[1] == 3 else []
    
    if len(dim) == 2:
        fig, ax = plt.subplots()
        ax.imshow(img.T, cmap='gray', clim=img_clim)
        
        if cell_boundary:
            plotCellBoundary(cell_boundary, color_cell, plot_param['CellBoundaryWidth'])
        
        if pos3d_label:
            label_cell_center(pos3d[:, [1, 0]], pos3d_label, color_cell)
        
        ax.set_xticks([])
        ax.set_yticks([])
        plt.show()
        return fig, ax
    
    elif len(dim) == 3:
        fig, axes = plt.subplots(1, dim[2], figsize=(dim[2] * 5, 5))  # Adjust figure size
        
        for iplane in range(dim[2]):
            #print(f"Processing plane {iplane}, image shape: {img.shape}")
            ax = axes[iplane] if dim[2] > 1 else plt.gca()
            ax.imshow(img[:, :, iplane], cmap='gray', vmin=img_clim[0], vmax=img_clim[1])
            

            if pos3d.size > 0:
                indices = np.where(np.abs(pos3d[:, 2] - zdepth[iplane]) < 0.1)[0]
                
                if len(indices) > 0 and cell_boundary:
                    plotCellBoundary([cell_boundary[i] for i in indices], color_cell[indices, :], plot_param['CellBoundaryWidth'])
                
                if len(indices) > 0 and pos3d_label:
                    label_cell_center(pos3d[indices][:, [1, 0]], [pos3d_label[i] for i in indices], color_cell[indices, :])
                
                if plot_param['PlotCenter']:
                    cell_center = pos3d[indices][:, [1, 0]]
                    celln = color_cell.shape[1] if color_cell.ndim > 1 else 1
                    if cell_center.shape[1] == 2 and celln > 1:
                        if cell_center.ndim == 1:
                            cell_center = np.append(cell_center, 0)
                        else:
                            cell_center = np.hstack((cell_center, np.zeros((len(cell_center), 1))))
                        
                        cellN = len(cell_center)  # Ensure cellN is properly initialized
                        plotCellCenter(ax, cell_center, np.full(cellN, 9) if isinstance(9, (int, float)) else 9, color_cell[indices, :], plot_param['CellCenterWidth'])
            
            ax.set_xticks([])
            ax.set_yticks([])
    
    for ax in axes.flat:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_facecolor('black')  # Set background to black to blend in
    fig.patch.set_alpha(0)  # Fully transparent figure background
    return fig, axes

def label_cell_center(positions, labels, color_cell):
    for (x, y), label, color in zip(positions, labels, color_cell):
        plt.text(x, y, str(label), color=color, fontsize=8, ha='center', va='center')

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


def amp_normalize_dim(data, norm_dim, prc_th):
    """
    Normalize each slice independently along a specified dimension.

    Parameters:
    - data: 3D NumPy array (e.g., n x height x width).
    - norm_dim: The dimension along which to normalize (0, 1, or 2).
    - prc_th: Tuple (low percentile, high percentile) for normalization.

    Returns:
    - Normalized 3D NumPy array.
    """
    sub_data = get_data_dim(data, norm_dim)  # Extract slices along norm_dim

    # Normalize each slice independently
    for i in range(len(sub_data)):
        sub_data[i] = amp_normalize(sub_data[i], prc_th)

    return reconstruct_from_data_dim(sub_data, norm_dim, data.shape)  # Rebuild full data

def get_data_dim(data, norm_dim):
    """
    Extract slices along a specified dimension.

    Returns:
    - List of 2D slices along the specified dimension.
    """
    return [np.take(data, i, axis=norm_dim) for i in range(data.shape[norm_dim])]

def reconstruct_from_data_dim(sub_data, norm_dim, original_shape):
    """
    Reconstructs a 3D array from slices.

    Returns:
    - Reconstructed NumPy array.
    """
    return np.stack(sub_data, axis=norm_dim).reshape(original_shape)

def amp_normalize(slice_data, prc_th):
    """
    Apply percentile-based normalization to a 2D slice.

    Returns:
    - Normalized slice.
    """
    low_t, high_t = np.percentile(slice_data, prc_th)

    if high_t - low_t < 1e-8:
        high_t = low_t + 1e-8  # Prevent division by zero

    slice_data = np.clip(slice_data, low_t, high_t)  # Clip extreme values
    return (slice_data - low_t) / (high_t - low_t)  # Normalize



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





def gaussian_smooth_3D(image, sigma, smooth_dim=0):
    """
    Applies 2D Gaussian smoothing to each slice along the specified dimension of a 3D image.

    Parameters:
    - image: 3D NumPy array (n x height x width).
    - sigma: Gaussian smoothing parameter (standard deviation).
    - smooth_dim: The dimension along which to apply 2D smoothing (default: 0).

    Returns:
    - smoothed_image: 3D NumPy array after applying Gaussian smoothing.
    """
    smoothed_image = np.zeros_like(image)

    if smooth_dim == 0:
        for i in range(image.shape[0]):  # Loop over first dimension (n)
            smoothed_image[i, :, :] = gaussian_filter(image[i, :, :], sigma=sigma)
    elif smooth_dim == 1:
        for i in range(image.shape[1]):  # Loop over second dimension (height)
            smoothed_image[:, i, :] = gaussian_filter(image[:, i, :], sigma=sigma)
    elif smooth_dim == 2:
        for i in range(image.shape[2]):  # Loop over third dimension (width)
            smoothed_image[:, :, i] = gaussian_filter(image[:, :, i], sigma=sigma)
    else:
        raise ValueError("Invalid smooth_dim. Must be 0, 1, or 2.")

    return smoothed_image


##Not used
def smooth3D_dim1(data, smooth_par):
    """
    Applies a smoothing function along the third dimension of a 3D NumPy array.

    Parameters:
    - data: 3D NumPy array (nPlane or nFrame, Height, Width)
    - smooth: Function to apply for smoothing 

    Returns:
    - smoothed: 3D NumPy array after applying smoothing function.
    """
    smoothed = np.zeros_like(data)

    for i in range(data.shape[2]):  # Loop over the third dimension
        smoothed[i,:, :] = smooth(data[i,:, :],smooth_par)

    return smoothed


def smoothDec(data, smooth):
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
##Not used

def apply_affine_transform_to_all_frames(frames, yxShift, planeID=None):
    """
    Apply the same affine transform to all frames in a 3D stack.

    Parameters:
    - frames: A 3D numpy array of shape (n_frames, height, width).
    - yxShift: A 2D array with shape (n_planes, 2), where each row represents the [y_shift, x_shift] for each plane. If planeID is None, yxShift is expected to be a 1x2 vector.
    - planeID: (optional) The plane ID that determines which shift to apply. If None, yxShift is used as a single shift value.

    Returns:
    - Transformed frames as a 3D numpy array of the same shape.
    """
    # Determine the shift values
    if planeID is not None:
        shift_y, shift_x = yxShift[planeID]
    else:
        shift_y, shift_x = yxShift  # Assume yxShift is a 1x2 vector if planeID is None

    # Define the transformation matrix for OpenCV's warpAffine
    correctShift_mat = np.array([
        [1, 0, shift_x],  # Note the order: OpenCV uses (x, y) not (y, x)
        [0, 1, shift_y]
    ], dtype=np.float32)

    # Invert the sign to make sure the direction matches `stat_pixel_shift`
    correctShift_mat[0, 2] = -correctShift_mat[0, 2]
    correctShift_mat[1, 2] = -correctShift_mat[1, 2]

    # Get the dimensions of a frame
    n_frames, height, width = frames.shape

    # Apply warpAffine to all frames using a list comprehension and stack them back into a single array
    transformed_frames = np.stack(
        [cv2.warpAffine(frame, correctShift_mat, (width, height), borderMode=cv2.BORDER_CONSTANT, borderValue=0) for frame in frames]
    )

    return transformed_frames

def apply_affine_transform_to_ImgPlane(ImgPlane, yxShift):
    """
    Apply specific affine transform to each frame in a 3D stack based on its index.

    Parameters:
    - ImgPlane: A 3D numpy array of shape (n_ImgPlane, height, width) or a 2D array of shape (height, width).
    - yxShift: A 2D array with shape (n_ImgPlane, 2) or a 1D array of shape (2) if ImgPlane are 2D.

    Returns:
    - Transformed ImgPlane as a 3D numpy array of the same shape if input is 3D,
      or a 2D array if the input is 2D.
    """
    # Handle the case where ImgPlane is a 2D array
    is_2d = False
    if ImgPlane.ndim == 2:
        ImgPlane = np.expand_dims(ImgPlane, axis=0)
        yxShift = np.expand_dims(yxShift, axis=0)
        is_2d = True

    n_ImgPlane, height, width = ImgPlane.shape

    transformed_ImgPlane = np.stack([
        cv2.warpAffine(
            frame, 
            np.array([
                [1, 0, -yxShift[frameID, 1]],  # Apply x shift (negative sign for correct direction)
                [0, 1, -yxShift[frameID, 0]]   # Apply y shift (negative sign for correct direction)
            ], dtype=np.float32),
            (width, height),
            borderMode=cv2.BORDER_CONSTANT,
            borderValue=0
        ) for frameID, frame in enumerate(ImgPlane)
    ])

    # Return as 2D array if the input was 2D
    if is_2d:
        return transformed_ImgPlane[0]

    return transformed_ImgPlane

def Pixel_Shift(rmg1,rmg2,nplanes):

    pixel_shifts = []
    if nplanes >1:
        for plane_idx in range(nplanes):
            shift, junk1, junk2 = phase_cross_correlation(rmg1[plane_idx, :, :], rmg2[plane_idx, :, :])
            pixel_shifts.append(shift)
    else:
        pixel_shifts, junk1, junk2 = phase_cross_correlation(rmg1, rmg2)

    pixel_shifts = np.array(pixel_shifts)
    return pixel_shifts




def find_all_folders_keywords(work_folder, folder_keyword, subdepth=None):
    """
    Find all folders in 'work_folder' whose names contain 'folder_keyword'.

    Parameters
    ----------
    work_folder : str
        Base directory path to search (UNC/network paths OK).
    folder_keyword : str
        Substring to match in folder names (case-sensitive).
    subdepth : int | None, optional
        Maximum depth (in folder levels below 'work_folder') to traverse.
        - 1 => only immediate subfolders of work_folder (e.g., 'SL0864', 'SL0893').
        - 2 => immediate subfolders AND their children (e.g., the date folders under 'SL0864').
        - None (default) => search all nested subfolders with no depth limit.

    Returns
    -------
    list[str]
        Full paths to matching folders.

    Notes
    -----
    Depth is counted in directory levels below 'work_folder'. No files are returned.
    """
    if subdepth is not None:
        if not isinstance(subdepth, int) or subdepth < 1:
            raise ValueError("subdepth must be a positive integer or None")

    matches = []
    # topdown=True allows us to prune traversal by editing 'dirs' in-place
    for root, dirs, _ in os.walk(work_folder, topdown=True):
        # current depth of 'root' relative to base (0 at base)
        rel = os.path.relpath(root, work_folder)
        depth = 0 if rel == '.' else rel.count(os.sep) + 1

        # Record matches among immediate children of 'root'
        for d in dirs:
            if folder_keyword in d:
                matches.append(os.path.join(root, d))

        # If we've reached the maximum allowed depth, stop descending further
        if subdepth is not None and depth >= subdepth:
            dirs[:] = []  # prune deeper traversal

    return matches



def get_exp_data_folder(work_folder, folder_keyword, included_list):
    """
    Function to find an experimental data folder containing specific files.
    
    Inputs:
        work_folder (str): Root directory to search in.
        folder_keyword (str): Keyword to filter folders.
        included_list (list): List of required files/folders within the target folder.
    
    Output:
        str: Path to the first folder that meets the criteria, empty if none found.
    
    Example usage:
        workspace_folder = get_exp_data_folder('C:/Data', 'Experiment', ['Data', 'Results.mat'])
    """
    workspace_folder = ""
    result_paths = find_all_folders_keywords(work_folder, folder_keyword)
    
    for folder in result_paths:
        # Ensure folder is not empty
        if len(os.listdir(folder)) > 2:
            items = {item.name for item in os.scandir(folder)}  # Convert to set for fast lookup
            if sum(1 for pattern in included_list if any(pattern in item for item in items)) >= len(included_list):
                return folder  # Return the first matching folder
    
    return workspace_folder  # Return empty string if no folder matches

