function imgRead_rioMakrPoistion()

processedDir =[ 'D:\h06142022_2\LS408\TSeries-06142022-0918-003'];
% Dir_path=['D:\h06132022_2\06132022_SL172\TSeries-06132022-0928-002'];
% cd(Dir_path)
% files = dir('*.tif');
% 
% for i=1:402:length(files)
%      imageName=files.name(i);    
% end

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
	fprintf('Averaging images in the following folder:\n          %s', folder);
end

%--------------------------------------------------------------------------------
% Get a list of all TIF, PNG, BMP, and JPG files.
imageFiles = [dir(fullfile(folder,'*.TIF')); dir(fullfile(folder,'*.PNG')); dir(fullfile(folder,'*.BMP')); dir(fullfile(folder,'*.jpg'))];
%%%for to know the total plan and how many avg
fileparts_end=imageFiles(end).name;
fileparts_end= regexp(fileparts_end,'\d+','match');
fileparts_end=str2double(fileparts_end);
%--------------------------------------------------------------------------------
% Now do the averaging in a loop
numberOfImages = length(imageFiles);
theyreColorImages = false;
kcount=1;
%%%
for ii=1:numberOfImages
    fullFileName = fullfile(folder, imageFiles(ii).name);
	fprintf('About to read %s\n', fullFileName);
	thisImage=imread(fullFileName);
    thisImage_total(:,:,ii)=thisImage;
   % kcount=kcount+1;
end
%%%%%%%%%%%%%%%%%%%%%
%thisImagetest=imread('C:\Data\Hari\08112022\HS8999\max_projection_withoutMC.tif');
ScanAmp_V_FOV_X=[-2.085852 2.085852];
ScanAmp_V_FOV_Y=[-2.332365 2.302270];
FOVsize_PX=512;
LUTx = linspace(ScanAmp_V_FOV_X(1), ScanAmp_V_FOV_X(2), FOVsize_PX);
LUTy = linspace(ScanAmp_V_FOV_Y(1), ScanAmp_V_FOV_Y(2), FOVsize_PX);

%Pointsxy=load('C:\Data\Hari\08112022\POINT.mat');
Pointsxy=load('E:\09022022\TSeries-09022022-1211-001\New folder\1\POINT.mat');
PointsX=Pointsxy.points.X;
PointsY=Pointsxy.points.Y;
%figure;
%imshow(thisImagetest,[]);
figure;
imshow(thisImage_total(:,:,1),[]);
hold on
for j=1:size(PointsX,2)
    Xtemp=PointsX(j);
    Ytemp=PointsY(j);
    %imshow(thisImage,[])
    %hold on
    plot(Xtemp,Ytemp,'Or','MarkerSize',10)
%hold on
end

%%%%%%%%%%%%%%%%%%%%%

%%%%



x