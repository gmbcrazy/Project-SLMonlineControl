XMLTable=[XMLTable;tempXMLTable];
[indexVector, stimulusIDVector, PrePostStimuliVector] = getPSTHFrames(PVparam.InterMPRepetition, PreSLMCal, PostSLMCal);
tempPSTHmap = CalMultiPSTHBin(ExpFileInfo(CountExp).binFile, confSet, indexVector, stimulusIDVector, prePostStimuliVector);
PSTHmap=cat(PSTHmap,tempPSTHmap,3);
CountExp=CountExp+1;
OutPng=[SumDataFolder 'CheckFunExp' num2str(CountExp) '.png'];



[averagedPSTHmap,nSample] = AveragePSTHByGroupLaser(XMLTable, PSTHmap,TotalGroupIDs);
Zdepth=confSet.ETL;



figure;
for iFun = 1:length(nSample)
    if iFun<length(nSample)
       FinalFunScore(Group(iGroup).Indices, 1) = iGroup;
       Pos3D=FinalPos3D(Group(iFun).Indices,:);
       RealMPInd=find((isnan(FinalFunScore(Group(iGroup).Indices,2)))==0);
       AddMPInd=find((isnan(FinalFunScore(Group(iGroup).Indices,2)))==1);

     for iplane = 1:nPlane
        if nSample(iFun)>0
           subplotLU(length(nSample), nPlane, iFun, iplane);
            
            % Display the current plane's image
            imagesc(averagedPSTHmap(:,:,iplane,iFun)');
          
            I = find(abs(Pos3D(:,3) - Zdepth(iplane)) < 0.1);
            I1=intersect(I,RealMPInd);
            I2=intersect(I,AddMPInd);
            if ~isempty(I1)
                plotCellCenter(Pos3D(I1,[2 1]), 7, [0.1 0.9 0.1],1);    
            end

            if ~isempty(I2)
                plotCellCenter(Pos3D(I2,[2 1]), 7, [0.9 0.1 0.1],1);    
            end
        end
     end
    else
     for iplane = 1:nPlane
        if nSample(iFun)>0
           subplotLU(length(nSample), size(Img,3), iFun, iplane);
            
            % Display the current plane's image
            imagesc(averagedPSTHmap(:,:,iplane,iFun)');

        end
     end




    end


end

papersizePX=[0 0 nPlane*4 length(nSample)*4]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,OutPng,'png');

