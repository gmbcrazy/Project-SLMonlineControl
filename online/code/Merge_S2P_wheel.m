%%


function Merge_S2P_wheel()
% --- Get Suite2P data ---
clc;    % Clear the command window.
%close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

%--------------------------------------------------------------------------------
% Get folder name where the images are that we want to average.
% Option #1:
folder = fileparts(which('cameraman.tif')); % Determine where MATLAB demo images folder is.
% Option #2:
folder = 'D:\My Pictures\Illusions'; % Specify some particular folder.
% Option #3:
folder = uigetdir(pwd, 'Select folder');
% folder will be 0 (a double) if they click cancel.
% folder will be the path (a string) if they clicked OK.
if folder == 0
	% Clicked cancel.  Exit program.
	return;
end
% Comment out whichever folder selection options above that you don't want to use.

% Make sure the folder actually exists.
if ~isdir(folder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', folder);
	uiwait(warndlg(errorMessage));
	return;	
else
	fprintf('Ca++ information in following folder:\n          %s', folder);
end

Cal_varPath = [dir(fullfile(folder,'*.mat'))];

Cal_var=load([folder '\' Cal_varPath(end).name]);

%%--- Extract Suite2P parameters ---
F=Cal_var.F;
Fneu=Cal_var.Fneu;
spks=Cal_var.spks;
Fsubtracted = F - 0.7*Fneu;

% --- Calculate Baseline ---
%from suite2p
 [NT , ~] = size(Fsubtracted);
% determine and subtract the running baseline
 win         = 60; % in sec
 Fbase       = zeros(size(Fsubtracted), 'single');
 % if getOr(ops, 'runningBaseline', 0)
 Fbase   = Fsubtracted;
   %    if getOr(ops, 'zBaseline', 0)
   %        [~, ix] = sort(ops.zdrift);
    %        Fbase = Fbase(ix, :);
    %    end
    fR=30; %%fr in Hz
ntBase  = 2*ceil(win * fR/2)+1;
 Fbase   = cat(1, Fbase((ntBase-1)/2:-1:1, :), Fbase, Fbase(end:-1:end-(ntBase-1)/2, :));
        
 Fbase   = my_conv2(Fbase, 20, 1);
 Fbase   = movmin(Fbase, ntBase,1);
 Fbase   = movmax(Fbase, ntBase,1);
 Fbase   = Fbase((ntBase-1)/2 + [1:NT], :);
        
        %   if getOr(ops, 'zBaseline', 0)
        %       Fbase(ix,:) = Fbase;
        %   end
        % end
deltaFoF    = Fsubtracted - Fbase; % THIS IS NOT deltaFoF; I JUST DID NOT WANT TO CHANGE ALL MY CODES
        
        % Fsort     = my_conv2(F1, ceil(ops.fs), 1);
        % Fsort     = sort(Fsort, 1, 'ascend');
        % baselines = Fsort(ceil(NT/20), :);
        % Fbase2    = bsxfun(@times, ones(NT,1), baselines);
        % F1        = F1 - Fbase2;
        % Fbase     = Fbase + Fbase2;
        
        % F1        = F1./max(std(F1,1,1), Fbase);
        
        % --- Normalization ---
        % Option 1 - DeltaF/F (can lead to a lot of distortions)
        % deltaFoF = zeros(size(Fsubtracted));
        % for c = 1:numberOfCells
        %   deltaFoF(:,c) = (Fsubtracted(:,c)-(F0(:,c)))./abs(F0(:,c));
        % end
        %
        % Norm - option 2 - From suite2p
        % normalize signal
  sd         = 1/sqrt(2) * std(deltaFoF(2:end, :) - deltaFoF(1:end-1, :), [], 1);
  deltaFoF   = bsxfun(@rdivide, deltaFoF , 1e-12 + sd); % THIS IS NOT deltaFoF; I JUST DID NOT WANT TO CHANGE ALL MY CODES

        % --- Append vars ---
  %save([folder ], 'iscell', 'redcell', 'F', 'Fneu', 'Fbase', 'deltaFoF', 'spks', 'fR', '-append');
  %%-------speed detect--------
   excelfilename='TSeries-10192022-1257-003_Cycle00001_VoltageRecording_001';
  L = readmatrix([folder '\' excelfilename]);
  L=L(1:end-50001,:);%% need to check later about it
 fr=30;
 srateV=10000;
 tcal=(1:size(deltaFoF,2))/fr;
 twmov=(1:size(L,1))/srateV;
  %%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%
avg_win=333;
loop1=1;
LmovV=L(:,2);

for i=1:333:size(L,1)-333
Lmovtemp=LmovV(i:i+333-1,1);
LmovVT(loop1)=mean(Lmovtemp);
loop1=loop1+1;
end
t
%%%%
 

