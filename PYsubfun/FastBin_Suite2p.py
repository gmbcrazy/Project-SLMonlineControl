#!/usr/bin/env python
# coding: utf-8

# 03252024, Lu Zhang for fast processing .bin data recorded from PrairieLink using Suite2p package without motion correction

# Import necessary packages for data manipulation, plotting, file management, and more.
import suite2p  # Suite2p for calcium imaging data analysis
import numpy as np  # NumPy for numerical operations
import os  # OS module for interacting with the file system
import glob  # Glob for filename pattern matching
from natsort import natsorted  # Natsort for natural sorting
 # Explicit garbage collection to free memory


def LoadBin(binFile,SaveFolder,ops1):
    nplanes=ops1['nplanes']
    if not os.path.exists(SaveFolder):
       os.makedirs(SaveFolder)
    else:
       print("Folder already exists.")
       os.path.exists(binFile[0])
    if len(binFile)>1:
       print('More than one bin files found!!')
    else:
       rawBin = suite2p.io.BinaryFile(Ly=ops1['Ly'],Lx=ops1['Lx'], filename=binFile[0])
       nFrame = rawBin.shape[0] // nplanes

#Make sure each plane has same number of frames.
    FramePerPlane=np.floor(rawBin.shape[0]/ops1['nplanes'])
    TotalFrameNeed=np.int32(np.floor(FramePerPlane)*ops1['nplanes'])
    print(f"Total Frames {TotalFrameNeed} is found")
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
    plane_data = rawBin[range(0 + plane_idx, rawBin.shape[0], nplanes), :, :]
    # Perform cell detection using suite2p
    ops, stat = suite2p.detection_wrapper(f_reg=plane_data, ops=ops1)
    # Update ops dictionary with additional information
    nFrame = rawBin.shape[0] // nplanes
    ops['nframes'] = nFrame
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

    if not os.path.exists(fpath):
        os.makedirs(fpath)
    else:
        print("Folder already exists.")

        np.save(os.path.join(fpath, "F.npy"), F)
        np.save(os.path.join(fpath, "Fneu.npy"), Fneu)
        np.save(os.path.join(fpath, "spks.npy"), spks)
        np.save(os.path.join(fpath, "ops.npy"), ops)
        np.save(os.path.join(fpath, "stat.npy"), stat)
        np.save(os.path.join(fpath, "iscell.npy"), iscell)
        ops_matlab=ops;
        ops_matlab['save_path']=fpath;
        suite2p.io.save_mat(ops, stat, F, Fneu, spks, iscell,redcell=-1)
    


def PostMannual(SaveFolder, ops0):

    CombinePath=os.path.join(SaveFolder, 'combined');
    iscellCombined = np.load(os.path.join(CombinePath, "iscell.npy"))
    statCombined = np.load(os.path.join(CombinePath, "stat.npy"),allow_pickle=True)
    opsCombined = np.load(os.path.join(CombinePath, "ops.npy"),allow_pickle=True)
    print(len(statCombined),'of units including', np.int16(np.sum(iscellCombined[:,0])),'identified cells in combined planes')
    nplanes=ops0['nplanes']
##extract the plane ID into a numpy variable
    UnitPlane=[];
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




