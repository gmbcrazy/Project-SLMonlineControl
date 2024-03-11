%% Process calcium data Part 2
% --- Ana last revision 08/24/2022 ---

clear
close all

% --- Set paths ---
processedDir        = 'E:\Hari\11282022\statest\';%'C:\Data\Hari\statest3\';
%processedDirWhisk   = 'C:\Data\Hari\whisks\';%'C:\Data\Hari\whisks\';
excelfilename       = 'C:\test\2022_VIPSSTSummaryCOPY.xlsx';
processedDirSuite2P = 'E:\Hari\11282022\suite2pfinal\';% 'C:\Data\Hari\suite2pfinal\';
expType             = 'chronic';
expBatch            = 'chronic';
rowsOfInterest      = [18]; % % Exceptions lines 9-21-24-49; CHANGE HERE; BATCH PROCESSING


% ---
% processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\Botox_Analyzed\';
% processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\Botox_AnalyzedWhisk\';
% excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2021_Opto&BotoxSummary.xlsx';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\2_Botox\';
% expType             = 'botox';
% expBatch            = 'botox';
% rowsOfInterest      = [39 41 42 43 44 47 49 50 55 57]; % [23 25 36 37] included in 'opto' at this point (I then did some folder transfer); no exceptions; CHANGE HERE; BATCH PROCESSING
% rowsOfInterest      = [25 37 39 41 42 43 47 49 55 57]; % at some later point when I decided to run group


% ---
% processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\Neuromod_Analyzed\';
% processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\Neuromod_AnalyzedWhisk\';
% excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2021_NeuromodSummary.xlsx';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\4_Neuromod\';
% expType             = 'neuromod';
% expBatch            = 'neuromod';
% rowsOfInterest      = [3 5 7 18 19 21 22 35 36 41 43 45 46 51 53 55 59 61 63 64]; % Exceptions lines 18, then 21 (not handled yet) and 57 (not handled yet); CHANGE HERE; BATCH PROCESSING; 


% ---
% processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\SCE_Analyzed\';
% processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\SCE_AnalyzedWhisk\';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\5_Starters\';
% excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2021_SCESummary.xlsx';
% expType             = 'SCE';
% expBatch            = 'SCE';
% rowsOfInterest      = [16:18 20 22 24:25 27:34 36:38 40:41 43 45 47 49]; % Exception line 17; CHANGE HERE; BATCH PROCESSING


% ---
% processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\Opto_Analyzed\';
% processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\Opto_AnalyzedWhisk\';
% excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2021_Opto&BotoxSummary.xlsx';
% expType             = 'opto';

% expBatch            = 'archt_tdTom1';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\3_Opto1\';
% rowsOfInterest      = [16 18 20 22 23 24 25 27 28 30 33 35 36 37]; % [44 50] included in 'botox' at this point (I then did some folder transfer) CHANGE HERE; BATCH PROCESSING
% 
% expBatch            = 'archt_tdTom2';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\1_Opto2\';
% rowsOfInterest      = [66:67 69 71:72 75:77 80 82]; % Exception lines 66 (not handled yet); CHANGE HERE; BATCH PROCESSING
% 
% expBatch            = 'jaws';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\9_Opto3\'; 
% rowsOfInterest      = [92 94]; % CHANGE HERE; BATCH PROCESSING
% 
% expBatch            = 'archt_tdTom3';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\10_Opto4\';
% rowsOfInterest      = [89 99 100 103 105 106 109 111 112 113]; % Exception lines 99; CHANGE HERE; BATCH PROCESSING


% ---
% processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\Dreadds_Analyzed\';
% processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\Dreadds_AnalyzedWhisk\';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\8_Dreadds\';
% excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2022_DreaddsSummary.xlsx';
% expType             = 'neuromod';
% expBatch            = 'dreadds';
% rowsOfInterest      = [9]; % CHANGE HERE; BATCH PROCESSING


% ---
% processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\SST_Analyzed\';
% processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\SST_AnalyzedWhisk\';
% processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\11_SST_Ana\';
% excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2022_InterneuronsSummary.xlsx';
% expType             = 'chronic';
% expBatch            = 'SST';
% rowsOfInterest      = [2:7]; % CHANGE HERE; BATCH PROCESSING

%% --- Merge Ana and Suite2P data, and calculate deltaFoF -
% Get F, Fneu, spk, redCell from Suite2P ouptup, and add to Ana vars file

Merge(processedDir, processedDirSuite2P, excelfilename, rowsOfInterest)

%% --- Align Ca/Locomotion and Whisk files ---
% This function assumes that all previous processing is correct
% Will take into account as trial number the lenght of the variable vRec

%AlignBatch(processedDir,processedDirWhisk, excelfilename,
%rowsOfInterest)%%original
AlignBatch(processedDir, excelfilename, rowsOfInterest)%%By Hari
%% --- Verify cam delay ---
% In the future, add warning for short cam delays; save manually

%[camDelays] = checkCamDelay(processedDir, excelfilename, rowsOfInterest,
%expBatch);%%By Hari 104

%% --- Get experimental groups ---
% Inspect all data; remove any unwanted trials from groups manually for the
% moment (excpetions indicated above and in ProcessedData_Inspection.docx)
%groupTrials(processedDir, excelfilename, rowsOfInterest);%%By Hari

%% --- Some tests ---
% % clear
% % close all
% % 
% % processedDir        = 'Z:\UFNC15\suite2p&AnaFINAL\Opto_AnalyzedNeuropil09\';
% % processedDirWhisk   = 'Z:\UFNC15\suite2p&AnaFINAL\Opto_AnalyzedWhisk\';
% % suite2p09           = 'Z:\UFNC15\suite2pFINALOPTO09\';
% % excelfilename       = 'Z:\UFNC15\suite2p&AnaFINAL\2021_Opto&BotoxSummary.xlsx';
% % 
% % expBatch            = 'archt_tdTom1';
% % processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\3_Opto1\';
% % rowsOfInterest      = [16 18 20 22 23 24 27 28 30 33 35 36];
% % Merge09(processedDir, processedDirSuite2P, suite2p09, excelfilename, rowsOfInterest)
% % AlignBatch(processedDir, processedDirWhisk, excelfilename, rowsOfInterest)
% % 
% % expBatch            = 'botox';
% % rowsOfInterest      = [44 50];
% % processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\2_Botox\';
% % Merge09(processedDir, processedDirSuite2P, suite2p09, excelfilename, rowsOfInterest)
% % AlignBatch(processedDir, processedDirWhisk, excelfilename, rowsOfInterest)
% % 
% % expBatch            = 'archt_tdTom2';
% % processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\1_Opto2\';
% % rowsOfInterest      = [66:67 69 71:72 75:77 80];
% % Merge09(processedDir, processedDirSuite2P, suite2p09, excelfilename, rowsOfInterest)
% % AlignBatch(processedDir, processedDirWhisk, excelfilename, rowsOfInterest)
% % 
% % expBatch            = 'archt_tdTom3';
% % processedDirSuite2P = 'Z:\UFNC15\suite2pFINAL\10_Opto4\';
% % rowsOfInterest      = 103; % [89 99 100 103 105 106 109 111 112 113];
% % Merge09(processedDir, processedDirSuite2P, suite2p09, excelfilename, rowsOfInterest)
% % AlignBatch(processedDir, processedDirWhisk, excelfilename, rowsOfInterest)

