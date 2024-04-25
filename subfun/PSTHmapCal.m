function PSTHtemp=PSTHmapSampleCal(BinFile,PSTHparam,confSet,stat)


         PreData = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PreInd);
         PostData= Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PostInd);
         PSTHtemp=squeeze(mean(PostData,3)-mean(PreData,3));
         PSTHtemp=SmoothDecDim3(PSTHtemp,PSTHparam.SmoothSD);


         

end