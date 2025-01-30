%% Process calcium data Part 1
% --- Ana last revision 01/16/2022 ---

% --- Prepare data for analysis ---
% 0) Experiment == Animal + Date + ExpDivision
% 1) Summarize experiments in a standardized excelsheet (text)
%    Go through the manual annotations, actual files (Fiji-based visual
%    inspection) and histology and determine experiments and trials to be
%    included in the analysis
%    
%    I have 2 types of excel files: with TseriesOpto_not_used (opto) and
%    with TseriesOpto (Neuromod and Dreadds)

% 2) Download experiments (Ca and Whisk folders from server) - be aware of
%    local drives capacity (run batches)
% 3) Check the Ca aquisition software version in the metadata
% 4) Convert Ca files using the software version-matched Bruker utility
% 5) Carefully set working directories and parameters

clear
close all

% --- Set paths ---
tempCaDir =  'C:\Data\Hari\'; % 'D:\06132022_2'% CHANGE HERE; BATCH PROCESSING
processedDir = 'C:\Data\Hari\statest'; % CHANGE HERE; BATCH PROCESSING

excelfilename = 'C:\test\2021Opto&BotoxSummary'; % CHANGE HERE; BATCH PROCESSING
rowsOfInterest = [13]; % CHANGE HERE; BATCH PROCESSING

%% --- Metadata extraction and simple motion correction ---

MotionCorrection2020MorePadding_updateHari(tempCaDir, processedDir, excelfilename, rowsOfInterest); % excelSheetNum,

% - Once ready, inspect each experiment's output; focus on identifying
% potential experimental issues (such as liquid loss - immersion lens -
% leading to decreased F, or Z-drifts; Z-drifts are generally the worst
% problem that can occur); decide if certain parts of the recording can be
% truncated or if the experiment has to be discarded
% - Proceed with Suite2P processing; so far, Ana has not written batch
% processing scripts for Suite2P; simply prepare the tempCaDir by removing
% all the data that should not be taken into consideration and do "manual"
% batch processing, 10 experiments at a time)
