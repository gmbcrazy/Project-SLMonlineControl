function ImageSeq=FrameIndMultiTiffs(TiffFolder,nPlane,Ind)

% TiffFolder='E:\LuRecording\03082024_MouseMei03_VGlutAi203SLM\TSeries-03082024-1105-003\'
% TiffFile=dir([TiffFolder '*TSeries*Ch2*.ome.tif'])

% % nPlane=3;
% for iPlane=1:nPlane
%     for iframe=1:length(TiffFile)
%         X(:,:,iframe,iPlane)=imread([TiffFile(iframe).folder '\' TiffFile(iframe).name],iPlane);
%     end
%     MeanImage(:,:,iPlane) = mean(X(:,:,:,iPlane),3);
% end
% 

% nPlane=3;
for iframe=1:length(Ind)
    % iframe
    TiffFile=dir([TiffFolder '*TSeries*Cycle*0' num2str(Ind(iframe)) '_Ch2_*.ome.tif']);
    if length(TiffFile)==0
       continue;
    end

    for iPlane=1:nPlane
        ImageSeq(:,:,iframe,iPlane)=imread([TiffFile(1).folder '\' TiffFile(1).name],iPlane);
    end
end

