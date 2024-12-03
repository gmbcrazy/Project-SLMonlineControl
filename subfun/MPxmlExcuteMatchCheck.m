function [XMLTable,UnMatchI]=MPxmlExcuteMatchCheck(FileGenerateInfo,XMLTable,AllTestPoints3D,confSet)

MPxmlFiles=dir([FileGenerateInfo.tifFolder '*MarkPoints*.xml']);
framesMarkStimuli=[];
clear EXcuteMP=[];
    for i = 1:length(MPxmlFiles)
        % Extract the cycle number from the MarkPoints.xml file name
        xmlFileName = MPxmlFiles(i).name;
        xmlTokens = regexp(xmlFileName, 'Cycle(\d{5})', 'tokens');
        if ~isempty(xmlTokens)
            markCycle = str2double(xmlTokens{1}{1});
            framesMarkStimuli=[framesMarkStimuli;markCycle];

            % Find the next cycle number after the MarkPoints cycle
%             nextCycle = min(find(cycleNumbers > markCycle));
            % if ~isempty(nextCycle)
            %     % Find the corresponding file for the next cycle and plane 001
            %     nextFile = sprintf('TSeries-%%*%05d_Cy*_%03d.ome.tif', nextCycle, uniquePlanes(1));
            %     nextFileMatch = dir([folderPath, nextFile));
            %     if ~isempty(nextFileMatch)
            %         framesAfterStimuli{end+1} = nextFileMatch.name;
            %     end
            % end
            [tbl,StimID{i},StimPowerTemp]=MPxml2Table([MPxmlFiles(i).folder '\' MPxmlFiles(i).name]);
            XYPosPixel=ExcuteMPxmlXYtoPixel([tbl.X(1) tbl.Y(1)],confSet);
            EXcuteMP(i,:)=[markCycle XYPosPixel];
        end
        % StimuliPower=[StimuliPower;MPxml2yaml([MPxmlFiles(i).folder '\' MPxmlFiles(i).name])];
        

    end


UnMatchI=find(sum(abs(AllTestPoints3D(XMLTable(:,2),1:2)-EXcuteMP(:,2:end)),2)>4);
if isempty(UnMatchI)
   disp('MarkPoints.Xml excuted successfully');
   UnMatchingFiles={};
else
   disp('UnMatching MarkPoints.Xml detected');
   UnMatchingFiles={MPxmlFiles.name};
   UnMatchingFiles=UnMatchingFiles(UnMatchI);
   for iXML=1:size(EXcuteMP,1)
       temp1=sum(abs(EXcuteMP(iXML,2:3)-AllTestPoints3D(:,1:2)),2);
       [temp2,I1]=min(temp1)
       if temp2<=4
          XMLTable(iXML,2)=I1;
       end
   end

end