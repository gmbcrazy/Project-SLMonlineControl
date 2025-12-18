#!/usr/bin/env python
# coding: utf-8

# 03252024, Lu Zhang for fast processing .bin data recorded from PrairieLink using Suite2p package without motion correction

# Import necessary packages for data manipulation, plotting, file management, and more.
#from __future__ import annotations
import suite2p  # Suite2p for calcium imaging data analysis
import numpy as np  # NumPy for numerical operations
import os  # OS module for interacting with the file system
import glob  # Glob for filename pattern matching
from pathlib import Path
from natsort import natsorted  # Natsort for natural sorting
 # Explicit garbage collection to free memory

from typing import Mapping, Any, Tuple
import pandas as pd
from scipy.interpolate import interp1d


#import Suite2p_QualityControl as QC 



def configLoad(file_path,ymlName,opsName):
    ymlPath=os.path.join(file_path,ymlName)
    opsPath=os.path.join(file_path,opsName)
    ymlSet=read_yaml(ymlPath)
    ops=np.load(opsPath,allow_pickle=True).item()
    ops=opsAddyaml(ymlSet,ops)
    return ops,ymlSet

def read_yaml(file_path):
    # Open the file and read lines
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Remove empty lines and strip whitespace
    lines = [line.strip() for line in lines if line.strip()]
    
    # Initialize the results dictionary
    results = {}
    
    # Parse each line
    for line in lines:
        # Ignore if this line is a comment
        if line.startswith('#'):
            continue
        
        # Find the separator between key and value
        sep_index = line.find(':')
        if sep_index == -1:  # Skip if no separator found
            continue
        
        # Extract key and value, trimming whitespace
        key = line[:sep_index].strip()
        value = line[sep_index+1:].split('#')[0].strip()  # Ignore comments
        
        # Attempt to convert value to a numeric type
        try:
            value = float(value) if '.' in value else int(value)
        except ValueError:
            pass  # Keep as string if conversion fails
        
        # Store the key and value
        results[key] = value
    
    return results

def opsAddyaml(yaml,ops):
    ops.update(yaml)
    ops['Lx']=yaml['SLM_Pixels_X']  ##Width of a image
    ops['Ly']=yaml['SLM_Pixels_Y']  ##Height of a image
    ops['xrange']=[0,ops['Lx']]
    ops['yrange']=[0,ops['Ly']]

    substrings=ops['ETL'].split()
    ops['ETL'] = [float(substring) for substring in substrings]

    ops['nplanes']=len(ops['ETL'])
    return ops


def suite2pInitiate(ops0,ConfigFolder):
    SaveFolder=os.path.join(ops0['save_path0'], 'suite2p')
    print('Processed data would be saved in'+ SaveFolder)
    SLMSettingPath = Path(ConfigFolder) / 'SLMsetting.yml'
    # Check if SLMsetting.yml exists
    if SLMSettingPath.exists():
        print(f"File found: {SLMSettingPath}")
    
    # Define destination folder and path
        CopySLMSettingPath = Path(ops0['save_path0']) / 'CurrentSLMsetting.yml'
        print(f"Copying to: {CopySLMSettingPath}")
    
    # Copy file
        CopySLMSettingPath.write_bytes(SLMSettingPath.read_bytes())
        print("Current SLMsetting copied successfully!")
    else:
        print("Error: SLMsetting.yml does not exist in the ConfigFolder.")
    
    if not os.path.exists(SaveFolder):
       os.makedirs(SaveFolder)
    else:
       print("Folder already exists.")
    return SaveFolder 


def LoadBin(binFile,ops1):
    nplanes=ops1['nplanes']
    if isinstance(binFile, str):
       rawBin = suite2p.io.BinaryFile(Ly=ops1['Ly'],Lx=ops1['Lx'], filename=binFile)
    elif not isinstance(binFile, str) and len(binFile)>1:
       print('More than one bin files found!!')
       return None
    else:
       rawBin = suite2p.io.BinaryFile(Ly=ops1['Ly'],Lx=ops1['Lx'], filename=binFile[0])
    
    nFrame = rawBin.shape[0] // nplanes

#Make sure each plane has same number of frames.
    FramePerPlane=np.floor(rawBin.shape[0]/ops1['nplanes'])
    TotalFrameNeed=np.int32(np.floor(FramePerPlane)*ops1['nplanes'])
    rawBin=rawBin[range(0,TotalFrameNeed),:,:]
    #print(f"Total Frames {TotalFrameNeed} is found")
    return rawBin, FramePerPlane, TotalFrameNeed

# Define the process_plane function to process individual imaging planes
def process_plane(plane_idx, rawBin, SaveFolder, ops1):
    nplanes=ops1['nplanes']
    print('Processing plane ' + str(plane_idx))
    # Construct the file path for saving data of the current plane
    fpath = os.path.join(SaveFolder, f'plane{plane_idx}/')
    # Create the save directory if it does not exist
    if not os.path.exists(fpath):
        os.makedirs(fpath)
    else:
        print("Saving folder already exists.")
    # Extract the plane data from the raw binary file
    FramePerPlane=np.floor(rawBin.shape[0]/ops1['nplanes'])
    TotalFrameNeed=np.int32(np.floor(FramePerPlane)*ops1['nplanes'])
    plane_data = rawBin[range(0 + plane_idx, TotalFrameNeed, nplanes), :, :]
    # Perform cell detection using suite2p
    ops, stat = suite2p.detection_wrapper(f_reg=plane_data, ops=ops1)
    # Update ops dictionary with additional information
    #nFrame = rawBin.shape[0] // nplanes
    ops['nframes'] = FramePerPlane
    ops['meanImg'] = np.mean(plane_data, axis=0)
    ops = suite2p.registration.register.enhanced_mean_image(ops)
    # Save the ops dictionary to a numpy file
    np.save(os.path.join(fpath, "ops.npy"), ops)
    # Extract signals from identified cells
    stat_after_extraction, F, Fneu, F_chan2, Fneu_chan2 = suite2p.extraction_wrapper(stat, f_reg=plane_data, f_reg_chan2=plane_data, ops=ops1)
    if len(stat) > 0:
        # Classify cells
        classfile = suite2p.classification.builtin_classfile
        iscell = suite2p.classify(stat=stat_after_extraction, classfile=classfile)
        # Preprocess signals for deconvolution
        dF = F.copy() - ops['neucoeff'] * Fneu
        dF = suite2p.extraction.preprocess(
            F=dF,
            baseline=ops1['baseline'],
            win_baseline=ops1['win_baseline'],
            sig_baseline=ops1['sig_baseline'],
            fs=ops1['fs'],
            prctile_baseline=ops1['prctile_baseline']
        )
        # Identify spikes in the fluorescence signal
        spks = suite2p.extraction.oasis(F=dF, batch_size=ops1['batch_size'], tau=ops['tau'], fs=ops['fs'])
        # Save processed data to files
        np.save(os.path.join(fpath, "stat.npy"), stat)
        np.save(os.path.join(fpath, "F.npy"), F)
        np.save(os.path.join(fpath, "Fneu.npy"), Fneu)
        np.save(os.path.join(fpath, "iscell.npy"), iscell)
        np.save(os.path.join(fpath, "ops.npy"), ops)  # Note: ops.npy is saved again, possibly to update it with new information
        np.save(os.path.join(fpath, "spks.npy"), spks)
        # Prepare data for MATLAB compatibility
        ops_matlab = ops
        ops_matlab['save_path'] = fpath
        suite2p.io.save_mat(ops, stat, F, Fneu, spks, iscell, redcell=-1)
    else:
        print("No ROIs found, only ops.npy file saved")

# Define the CombinePlanes function to combine data from multiple planes after processing
def CombinePlanes(SaveFolder, ops0):
    # Combine the data only if there are multiple planes
    if ops0['nplanes'] > 1:
        result = suite2p.io.save.combined(SaveFolder, save=True)
        # Extract combined data
        stat, ops, F, Fneu, spks, iscell_0, iscell_1, redcell_0, redcell_1, hasred = result
        iscell = np.stack((iscell_0, iscell_1), axis=-1)
        fpath = os.path.join(SaveFolder, 'combined')
        print(ops0['nplanes'],' planes combined')
        ops['save_path0']=SaveFolder
    if not os.path.exists(fpath):
        os.makedirs(fpath)
    else:
        print("Folder already exists.")
        ops['save_path0']=SaveFolder
        np.save(os.path.join(fpath, "F.npy"), F)
        np.save(os.path.join(fpath, "Fneu.npy"), Fneu)
        np.save(os.path.join(fpath, "spks.npy"), spks)
        np.save(os.path.join(fpath, "ops.npy"), ops)
        np.save(os.path.join(fpath, "stat.npy"), stat)
        np.save(os.path.join(fpath, "iscell.npy"), iscell)
        ops_matlab=ops
        ops_matlab['save_path']=fpath
        suite2p.io.save_mat(ops, stat, F, Fneu, spks, iscell,redcell=-1)
    


def PostMannual(SaveFolder, ops0):

    CombinePath=os.path.join(SaveFolder, 'combined');
    iscellCombined = np.load(os.path.join(CombinePath, "iscell.npy"))
    statCombined = np.load(os.path.join(CombinePath, "stat.npy"),allow_pickle=True)
    opsCombined = np.load(os.path.join(CombinePath, "ops.npy"),allow_pickle=True)
    print(len(statCombined),'of units including', np.int16(np.sum(iscellCombined[:,0])),'identified cells in combined planes')
    nplanes=ops0['nplanes']
##extract the plane ID into a numpy variable
    UnitPlane=[]
    for index, Unit in enumerate(statCombined):
        UnitPlane.append(Unit['iplane'])
    UnitPlane=np.array(UnitPlane)


    statPlane=[None]*nplanes
    opsPlane=[None]*nplanes
    for plane_idx in range(nplanes):
    #print(plane_idx)
        fpath=os.path.join(SaveFolder, f'plane{plane_idx}/');
        iscell=iscellCombined[UnitPlane==plane_idx,:]
        np.save(os.path.join(fpath, "iscell.npy"), iscell) ##update the iscell for each plane from the combined data
        opsPlane[plane_idx]= np.load(os.path.join(fpath, "ops.npy"),allow_pickle=True)

    ##asign information of stat for each plane to the combined planes
        statPlane[plane_idx] = np.load(os.path.join(fpath, "stat.npy"),allow_pickle=True)
        statCombined[UnitPlane==plane_idx]=statPlane[plane_idx]

    stat=statCombined
    np.savez(os.path.join(CombinePath, "statUpdate.npz"), stat=statCombined,UnitPlane=UnitPlane)




def copy_and_merge_files(output_file, input_files):
    """
    Merges multiple binary files into a single binary file.

    Parameters
    ----------
    output_file: str
        Path to the output binary file.
    input_files: list of str
        Paths to the input binary files to merge.
    """
    with open(output_file, 'wb') as fout:
        for input_file in input_files:
            with open(input_file, 'rb') as fin:
                while True:
                    data = fin.read()  # Read in chunks of 64KB
                    if not data:
                        break
                    fout.write(data)

    print(f"Merged binary file saved to {output_file}")

import os
import numpy as np
from suite2p.io.tiff import imread

def load_ref_ch(data_folder, ChID=1):
    """
    Load a specific reference channel (1 or 2) from multi-channel TIFF files using Suite2p.
    
    Parameters:
    - data_folder (str): The folder containing TIFF files.
    - ChID (int): The channel ID to load (1 or 2). Ch2 data is identified by filenames containing 'Ch2'.
    
    Returns:
    - channel_data (numpy.ndarray): Concatenated data from the specified channel across all TIFF files.
      The output shape is (nPlane, Ly, Lx).
    """
    if ChID not in [1, 2]:
        raise ValueError("ChID must be 1 or 2.")

    # List all TIFF files that match the specified channel
    tiff_files = [os.path.join(data_folder, f) for f in os.listdir(data_folder) 
                  if f.lower().endswith('.tif') and f'Ch{ChID}' in f]

    if not tiff_files:
        print(f"No TIFF files found for Channel {ChID} in the specified folder.")
        return None

    # Load the specified channel from all TIFF files
    channel_data = []
    for tiff_file in tiff_files:
        img_data = imread(tiff_file)  # Load the TIFF data
        if img_data is None:
            print(f"Failed to load image: {tiff_file}")
            continue

        # Handle multi-page or multi-channel TIFF files
        if img_data.ndim == 4:  # Expected shape: (num_frames, num_channels, height, width)
            plane_data = img_data[ChID-1, 0, :, :]  # Select the specified channel and the first plane only
            channel_data.append(plane_data)
        elif img_data.ndim == 3 and img_data.shape[0] in [2, 3]:
            img_data = img_data[ChID-1, :, :]  # Select the correct channel if multi-channel
            channel_data.append(img_data)
        elif img_data.ndim == 2:  # Already a 2D image
            channel_data.append(img_data)
        else:
            print(f"Unexpected image dimensions for {tiff_file}: {img_data.shape}")
            continue

    # Concatenate data along the plane dimension to get (nPlane, Ly, Lx)
    if channel_data:
        channel_data = np.stack(channel_data, axis=0)
        print(f"Loaded data shape: {channel_data.shape}, Expected nPlane: {len(tiff_files)}")
        return channel_data
    else:
        print("No valid channel data loaded.")
        return None


