import os
import re
from pathlib import Path
import xml.etree.ElementTree as ET
import shutil
from shutil import ignore_patterns
from datetime import datetime


class Animal:
    """
    Class for each animal profile including necessary information
    Here is a list of the variable and the associated column
    Birth is DOB
    Name is Animal 
    Files is Tseries not to be confused with files variable above
    Mag is Mag
    Age is Date2
    Recording Date is Date
    FileTime is the time associated with the label specific to each batch of recordings spit out
    """
    def __init__(self,birth,name,files,mag,age,recordingDate):
        self.birth = birth
        self.name = name
        self.files = files
        self.mag = mag
        self.age = age
        self.recordingDate = recordingDate


#Gets paths for XML files 
#Takes input lists of Animal type and folders (typically drives)
def getXMLPath(listOfAnimals: list , folderpath: list) -> list:
    allAnimalsList = []
    for folder in folderpath:
        for animal in listOfAnimals:
            date = str(int(animal.recordingDate))
            name = animal.name
            age = animal.age
            folderLocBeg = folder + date[:8] + "_" + name + "_" + age.lower()
            listOfFolders = scanDirTseries(folderLocBeg)
            allAnimalsList.append(listOfFolders)
    return allAnimalsList 


"""
Returns a list of all folders that are Tseries inside a specified folder
Takes in the animal folder name i.e. mouse9-2_p22 for example
Ignores any file or folder that does not begin with Tseries
"""
def getTseriesList(parentFolder: str):
    dir = parentFolder
    files = os.listdir(dir)
    filesClean = []
    for file in files:
        fullDir = Path(str(dir) + "/" + file)
        if (os.path.isdir(fullDir)) and (file[0:7].lower() == 'tseries'):
            filesClean.append(str(fullDir))
    return filesClean


"""
Slightly tweaked version of function above that gets Tseries locations
Function then gets the main xml file in those locations and appends them to a list
Similar to other function takes in a string of a folder
"""
def scanDirTseries(folder: str):
    dir = folder
    filesClean = []
    if os.path.isdir(dir):
        files = os.listdir(dir)
        for file in files:
            fullDir = Path(str(dir) + "/" + file)
            strDir = str(fullDir)
            fullDirFile = str(fullDir) + "/" + strDir[-25:] + ".xml"
            if (os.path.isfile(fullDirFile)):
                filesClean.append(str(fullDirFile))
    return filesClean


#Assign parameters here
def createListFilesAndParams(folderpaths: list, tiffFolders: list, savedisk=None) -> list:
    #Iterate over each individual folder location and then add parameters and folder location to list 
    folderpaths = folderpaths
    tiffFolders = tiffFolders
    db = []
    for drive in folderpaths:
        for folder in tiffFolders:
            filepath = Path(str(drive) + str(folder))
            savepath = Path("I:/Analyzed/" + str(folder) + "-results")
            if('savedisk' in locals()) and (savedisk != None):
                savepath = Path(str(savedisk) + "Analyzed/" + str(folder) + "-results")
            #Check to make sure directory exists
            if not os.path.isdir(filepath):
                print('Error: No Directory Exists at ' + str(filepath))
                continue
            else:
                prevDetectionsLoc = Path(str(savepath) + "/" + "suite2p")
                prevDetectionsMov = Path(str(savepath) + "_bak" + datetime.now().strftime("_%m-%d-%y--%H-%M-%S"))
                if os.path.isdir(prevDetectionsLoc) and not os.path.isdir(prevDetectionsMov):
                    shutil.copytree(prevDetectionsLoc,prevDetectionsMov, ignore=ignore_patterns('*.bin'))
                #Uncomment and add an entry to any parameter you wish to change
                parameters = {
                    'data_path' : [str(filepath)],
                    'look_one_level_down' : True,
                    #'fast_disk' : "/Users/holdennn/Documents/BinTemps",
                    #'delete_bin'
                    #'mesoscan'
                    #'bruker'
                    #'bruker_bidirectional'
                    #'h5py
                    #'h5py_key'
                    # nwb_file'
                    #'nwb_driver'
                    #'nwb_series'
                    'save_path0' : str(savepath),
                    #'save_folder'
                    #'subfolders'
                    #'move_bin' : True,
                    #'nplanes'
                    #'nchannels'
                    #'functional_channel'
                    'tau' : 1.5,
                    'fs' : 30,
                    #'force_sktiff'
                    #'frames_include'
                    #'multiplane_parallel'
                    #'ignore_flyback' : [-1],
                    #'preclassify' : 0.7,
                    'save_mat' : True,
                    #'save_NWB'
                    #'combined'
                    #'aspect'
                    #'do_biphase'
                    #'bidiphase'
                    #'bidi_corrected'
                    #'do_registration'
                    #'two_step_registration'
                    #'keep_movie_raw'
                    #'nimg_init' : 5000,
                    #'batch_size' : 4500,
                    #'maxregshift'
                    #'align_by_chan'
                    'reg_tif' : 1,
                    #'reg_tif_chan2'
                    #'subpixel'
                    #'smooth_sigma_time'
                    #'smooth_sigma'
                    #'th_badframes' : 0.9,
                    #'norm_frames'
                    #'force_refImg'
                    #'pad_fft'
                    #'nonrigid'
                    #'block_size'
                    #'snr_thresh'
                    #'maxregshiftNR'
                    #'1Preg'
                    #'spatial_hp_reg'
                    #'pre_smooth'
                    #'spatial_taper'
                    #'roidetect'
                    #'spikedetect'
                    #'sparse_mode'
                    #'spatial_scale'
                    #'connected'
                    #'nbinned' : 20000,
                    #'max_iterations' : 50,
                    #'threshold_scaling'
                    #'max_overlap'
                    #'high_pass'
                    #'spatial_hp_detect'
                    #'denoise'
                    #'anatomical_only'
                    #'diameter'
                    #'cellprob_threshold' : 0.2,
                    #'flow_threshold'
                    #'spatial_hp_cp'
                    #'pretrained_model'
                    #'soma_crop'
                    #'neuropil_extract'
                    #'inner_neuropil_radius'
                    #'min_neuropil_pixels'
                    #'lam_percentile'
                    #'allow_overlap'
                    'use_builtin_classifier' : True,
                    #'classifier_path', 
                    #'chan2_thres'
                    #'baseline'
                    #'win_baseline'
                    'sig_baseline' : 20.0,
                    #'prctile_baseline'
                    'neucoeff' : 0.7
                }
                db.append(parameters)
    return db


#Parses an XML files absolute path and returns the float value of that files magnification
def getMag(XMLFilePath) -> float:
    tree = ET.parse(XMLFilePath)
    value = None
    for shard in tree.iterfind('PVStateShard'):
        value = shard.find("PVStateValue[@key='opticalZoom']").attrib['value']
    return float(value)