function [fSpeed,fStim,timeStampCa_Plane]=PV_VolExtract_MultiFolder(confSet,fileID)


% PV_SpeedExtract extracts Voltage signals including speed, Stim, and timestamps from photovoltaic speed recording systems.
% This function aligns the extracted speed data with timestamps from calcium imaging data,
% facilitating the analysis of movement alongside neuronal activity across different planes.
%
% Inputs:
% confSet - A struct containing configuration settings and paths for data extraction.
%
% Outputs:
% fSpeed - A matrix containing the processed speed data for each imaging plane.
% fStim - A matrix containing the processed stim data (triggering whisking stimuli) for each imaging plane.
% timeStampCa_Plane - A matrix containing the timestamps for calcium imaging data for each plane.

% The following commented-out code block was intended for selecting a trial based on directory structure.
% It has been disabled for this version.

% Read XML configuration and trial data.
fSpeed=[];
fStim=[];
timeStampCa_Plane=[];

for iFile=1:length(fileID)
    [fSpeedTemp,fStimTemp,timeStampCa_PlaneTemp]=PV_VolExtract(confSet,fileID(iFile));
    fSpeed=[fSpeed;fSpeedTemp];
    fStim=[fStim;fStimTemp];
    timeStampCa_Plane=[timeStampCa_Plane;timeStampCa_PlaneTemp];
end