function varargout=PSTHmapSampleCal(BinFile,PSTHparam,confSet,varargin)

  if nargin>3
     ROIrowcol=varargin{1};
     CellPlane=varargin{2};
  end

         PreData = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PreInd);
         PostData= Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PostInd);
         PSTHtemp=squeeze(mean(PostData,3)-mean(PreData,3));
         PSTHtemp=SmoothDecDim3(PSTHtemp,PSTHparam.SmoothSD);
  varargout{1}=PSTHtemp;
  if nargin>3       
         PreCell=squeeze(PreData(:,:,:,CellPlane));
         PostCell=squeeze(PostData(:,:,:,CellPlane));
         linearIdx = sub2ind([confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X], ROIrowcol(:,1), ROIrowcol(:,2));

         for t = 1:length(PSTHparam.PreInd)
              PreTSeries(:, t) = PreCell(linearIdx + (t-1)*confSet.SLM_Pixels_Y*confSet.SLM_Pixels_X);
         end

         for t = 1:length(PSTHparam.PostInd)
              PostTSeries(:, t) = PostCell(linearIdx + (t-1)*confSet.SLM_Pixels_Y*confSet.SLM_Pixels_X);
         end

    ROITseries=[PreTSeries PostTSeries];
    varargout{2}=ROITseries;
  end
end