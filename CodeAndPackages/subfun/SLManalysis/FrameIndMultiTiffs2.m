function ImageSeq=FrameIndMultiTiffs2(TiffFolder,nPlane,Ind,Channel)

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
if nargin==3
   Channel=2;
end

for iframe=1:length(Ind)
    % iframe
    for iPlane=1:nPlane

        paddedStr = padNumber(Ind(iframe), 5);
        TiffFile=dir([TiffFolder '*TSeries*Cycle' paddedStr '_Ch' num2str(Channel) '_00000' num2str(iPlane) '.ome.tif']);
        if length(TiffFile)==0
            continue;
        end

        ImageSeq(:,:,iframe,iPlane)=imread([TiffFile(1).folder '\' TiffFile(1).name]);
    end
end

