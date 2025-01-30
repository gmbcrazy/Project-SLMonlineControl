%% Locomotion event onsets and offsets (special rules)
% --- Ana 08/24/2022 ---

clear
close all



%processedDir = 'Z:\UFNC15\suite2p&AnaFINAL\SCE_Analyzed\';
%excelfilename = 'Z:\UFNC15\suite2p&AnaFINAL\2021_SCESummary.xlsx';
%rowsOfInterest = [16:18 20 22 24:25 27:34 36:38 40:41 43 45 47 49]; % CHANGE HERE; ONE BY ONE OR BATCH

processedDir        = 'C:\Data\Hari\statest2\';
excelfilename       = 'C:\test\2022_VIPSSTSummaryCOPY.xlsx';
rowsOfInterest      = [17];

% --- Detection of locomotion events ---
[~, ~] = speedDetect2021(processedDir, excelfilename, rowsOfInterest);
%[~, ~, ~] = speedDetect2021(processedDir,
%excelfilename,rowsOfInterest);original
% [fTrial, fTrialBehav1, fTrialSpeed] = speedDetect2021(processedDir, excelfilename, rowsOfInterest);
% f = fields(fTrialSpeed);

%% --- Check ---
%%commented by Hari
% i = 1;
% 
% clf
% fieldname = char(f(i));
% plot(fTrial.(fieldname).fSpeed - min(fTrial.(fieldname).fSpeed))
% hold on
%plot(fTrial.(fieldname).fLEDCL)
%plot(fTrial.(fieldname).fWhisk1 / max(fTrial.(fieldname).fWhisk1)+2)
%Lines(fTrialBehav1.(fieldname).onsets2(fTrialBehav1.(fieldname).runInd == 1), [], 'g');
%Lines(fTrialBehav1.(fieldname).offsets2(fTrialBehav1.(fieldname).runInd == 1), [], 'g');
% Lines(fTrialSpeed.(fieldname).onsets1, [], 'y');
% Lines(fTrialSpeed.(fieldname).offsets1, [], 'y');
% Lines(fTrialSpeed.(fieldname).onsets3, [], 'k');
% Lines(fTrialSpeed.(fieldname).offsets3, [], 'r');
% Lines(fTrialSpeed.(fieldname).forPlotting, [0 1], 'm');

%ylim([0 5])
