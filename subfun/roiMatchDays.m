
function [Mapping,ROIInfo]=roiMatchDays(Data,Param)

% Description: This function performs registration and matching of regions of interest (ROIs) between consecutive days of imaging data.
% Input:
%   - Data: loaded from suite2p processed data, struct array containing imaging data for multiple days, including mean images, cell maps, and cell information.


%   - Param.OverlapTh, >= Param.OverlapThpix of overlap ratio to define initial legal matching cells 
%   - Param.PalphaTh, Significance Threshold P values to define matched ROIs across days.
%   - Param.PalphaTh, Significance Threshold P values to define matched ROIs across days.

%   - Param.PixFromMedCenter, Number of pixels from the median center of each ROI to get the ROI neighborhood.
%   - Param.NormPerc, % Set normalization percentiles for amplitude normalizationd.
                      % NormPerc=[20 99], lower and higher limit to
                      % calculate the amplitude, see AmpNormalize.m for more details
%   - Param.ShuffleNum, Define the number of permutations to perform

% Output:
%   - Mapping: Mapping matrix representing the correspondence between ROIs of consecutive days, based on significance threshold of Permutation test.
%   - ROIList: List of ROIs after matching, containing information about each ROI.

if nargin==1
  Param.OverlapTh=0.5;
  Param.PalphaTh=0.01;
  Param.PixFromMedCenter=20;
  Param.NormPerc=[20 99];
  % Param.ShuffleNum=10000;
end


% Set normalization percentiles for amplitude normalization
NormPerc=Param.NormPerc;
NanVec=[0 0];
%Used as Data = AmpNormalize(Data, prcTh, NanVec);

OverlapTh=Param.OverlapTh;
PixFromMedCenter=Param.PixFromMedCenter;


% DisX=[-1:0.05:1]; % Bins for distribution of ROI similarity measures;
% Palpha=[0.05 0.01 0.005 0.001]; %% P values to extract permutation test significance threshold of ROI neighourhood matching
% Th=prctile(corrPerm,(1-Palpha)*100);
% ShuffleNum=Param.ShuffleNum; % Define the number of permutations to perform
% PalphaTh=Param.PalphaTh;    % Significance Threshold P values to define matched ROIs across days.



% Loop through consecutive days for ROI matching
 % Set the first day as reference or Fixed day, register the next day (Moving day) to the Fixed day 
for FixI=1:length(Data)-1
    MappingTemp=[];

    MovingI=FixI+1;
     % Extract information for the fixed and moving days
    Fixed=Data(FixI).ops.meanImg;
    FixedCell=Suite2pCellIDMap(Data(FixI));
    FixedCellInfo=Suite2pCellInfo(Data(FixI));
    FixedCellID=1:sum(Data(FixI).iscell(:,1));


    Moving=Data(MovingI).ops.meanImg;
    MovingCell=Suite2pCellIDMap(Data(MovingI));
    MovingCellInfo=Suite2pCellInfo(Data(MovingI));
    MovingCellID=1:sum(Data(MovingI).iscell(:,1));  

    % Normalize images and display side by side for visual inspection
     Moving=AmpNormalize(Moving,NormPerc);
     Fixed=AmpNormalize(Fixed,NormPerc);

    % figure;
    % imshowpair(Fixed,Moving,"montage")


    % Perform image registration using histogram matching and demons algorithm
    Moving2 = imhistmatch(Moving,Fixed);
    [D,MovingReg] = imregdemons(Moving2,Fixed,[500 400 200],...
    'AccumulatedFieldSmoothing',1.3);  %%D is the non-linear transform of pixel from Moving to Fixed

    %% Using the D transform to transform Moving CellID map to Fixed imaging
    MovingCellReg = imwarp(MovingCell,D,'Interp','nearest','SmoothEdges',false);


    clear MovingFromFix OverlapPixN MappingPixN
    % Initialize the MappingPixN matrix to store the number of pixels of paired cells
    % 1st column: Number of pixels in the corresponding Fixed cell
    % 2nd column: Number of pixels in the corresponding Moving cell
    MappingPixN=zeros(length(MovingCellID),2);
 
    % Initialize MovingFromFix matrix to store the mapping relationship between Fixed and Moving cells
    % 1st column: Corresponding Fixed CellID
    % 2nd column: Moving CellID
    MovingFromFix(:,2)=MovingCellID(:);   

    % Loop through Moving cells for mapping information
    for ic=1:length(MovingCellID)
         % Identify pixels in Moving cell that correspond to the current MovingCellID
        CellInd=MovingCellReg==MovingCellID(ic);
         MappingPixN(ic,2)=sum(sum(CellInd));

         %
         FixedOverlap=FixedCell(CellInd); % Extract the corresponding Fixed cells based on the identifie pixels in the Moving cell;
         % multiple fixed cells could overlapwith Moving cells
         
         [OverlapCount,FixUnitID]=histcounts(categorical(FixedOverlap(:)')); % Count the occurrences of each Fixed cell in the overlap
         if ~isempty(OverlapCount)
         FixUnitID=cellfun(@str2double, FixUnitID);

         % Identify the Fixed cell with the maximum overlap count
         [OverlapPixN(ic,1),I]=max(OverlapCount); 
         FixUnitID=FixUnitID(I);
         MovingFromFix(ic,1)=FixUnitID;

         end
    end

%%Delete non-matching Cells
Invalid=find(MovingFromFix(:,1)==0);
MovingFromFix(Invalid,:)=[];
OverlapPixN(Invalid)=[];
MappingPixN(Invalid,:)=[];
%

%% sort Mapping pair with Fixed CellID
[~,I]=sort(MovingFromFix(:,1));
MappingTemp=MovingFromFix(I,:);
OverlapPixN=OverlapPixN(I);
% MappingPixN=MappingPixN(I,:);


MappingPixN(:,1)=double(FixedCellInfo.npix(MappingTemp(:,1)));
% MappingPixN(:,2)=double(MovingCellInfo.npix(MappingTemp(:,2)));

OverlapRmin=min(OverlapPixN./MappingPixN,[],2);

Invalid=OverlapRmin<OverlapTh;
MappingTemp(Invalid,:)=[];
OverlapRmin(Invalid)=[];

%% in case of a Moving Cell is asigned to >=2 cells in Fixed Cells
% For each Fixed Cell with multiple assignments, choose the one with the highest overlap ratio
clear tempC tempCell CheckCell
[tempC,tempCell]=histcounts(categorical(MappingTemp(:,1)));
CheckList=find(tempC>=2);
tempCell=cellfun(@str2double, tempCell);
CheckCell=tempCell(CheckList);
%Cells with highest overlap rato is chosen as the correct matching pair.
RepeatCellI=[];
if ~isempty(CheckCell)
   for ic=1:length(CheckCell)
       RepeatI=find(MappingTemp(:,1)==CheckCell);
       [~,I]=max(OverlapRmin(RepeatI))
       temp=setdiff(RepeatI,RepeatI(I));
       RepeatCellI=[RepeatCellI;temp(:)];
       clear temp
   end

   MappingTemp(RepeatCellI,:)=[];
   OverlapRmin(RepeatCellI)=[];
end

OverFixedCellInfo=FixedCellInfo(MappingTemp(:,1),:);
OverMovingCellInfo=MovingCellInfo(MappingTemp(:,2),:);

% Extract the neighborhood of the mapped Fixed and Moving cells.
% FixRoiN: Neighborhood of Fixed cells.
% FixRoiB: Boundary coordinates of Fixed cells.
[FixRoiN,FixRoiB]=ROI2neighbour(Fixed,FixedCell,MappingTemp(:,1),PixFromMedCenter);
[MovingRoiN,MovingRoiB]=ROI2neighbour(Moving,MovingCell,MappingTemp(:,2),PixFromMedCenter); % Extract the neighborhood of the registered Moving cells.

% Extract the neighborhood of the registered Moving cells.
[MovingRoiRegN,MovingRoiRegB]=ROI2neighbour(MovingReg,MovingCellReg,MappingTemp(:,2),PixFromMedCenter);


% Loop through each Moving cell and calculate the correlation with the corresponding Fixed cell.
clear corrMatch;
for i=1:size(MovingRoiN,3)
   [~,corrMatch(i)]=imregcorr(AmpNormalize(MovingRoiN(:,:,i)),AmpNormalize(FixRoiN(:,:,i)));
   temp1=MovingRoiN(:,:,i);
   temp2=FixRoiN(:,:,i);
   corrMatch(i)=corr(temp1(:),temp2(:));
end


% Extract the neighborhood of all Fixed cells and their boundary coordinates.
[FixRoiNAll,FixRoiBAll]=ROI2neighbour(Fixed,FixedCell,FixedCellID,PixFromMedCenter);
% Extract the neighborhood of all Moving cells and their boundary coordinates.
[MovingRoiNAll,MovingRoiBAll]=ROI2neighbour(Moving,MovingCell,MovingCellID,PixFromMedCenter);


Mapping{FixI}=MappingTemp;

ROIInfo(FixI).FixRoiN=FixRoiN;
ROIInfo(FixI).FixRoiB=FixRoiB;
ROIInfo(FixI).MovingRoiN=MovingRoiN;
ROIInfo(FixI).MovingRoiB=MovingRoiB;
ROIInfo(FixI).MovingRoiRegN=MovingRoiRegN;
ROIInfo(FixI).MovingRoiRegB=MovingRoiRegB;
ROIInfo(FixI).corrMatch=corrMatch;
ROIInfo(FixI).FixRoiNAll=FixRoiNAll;
ROIInfo(FixI).FixRoiBAll=FixRoiBAll;
ROIInfo(FixI).MovingRoiNAll=MovingRoiNAll;
ROIInfo(FixI).MovingRoiBAll=MovingRoiBAll;
ROIInfo(FixI).FixedCellID=FixedCellID;
ROIInfo(FixI).MovingCellID=MovingCellID;
ROIInfo(FixI).MovingImg=Moving;
ROIInfo(FixI).FixImg=Fixed;
ROIInfo(FixI).MovingImgReg=MovingReg;
ROIInfo(FixI).NonLinearTransform=D;



end




end