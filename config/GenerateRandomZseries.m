clear all
SampleFre=6.9;  %Sampling rate in Bruker  <<-------------------------------------------------------------------------------------Edited to adjust current imaging setting
% SaveFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\';
SaveFolder='C:\Users\User\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\';

ZIntBySecond=[6;10]; %Interval for a Zseries in Tsereis by seconds<<-------------------------------------------------------------Edited to adjust current experiment setting
ZIntByFrame=round(ZIntBySecond*SampleFre);%%Interval for a Zseries in Tsereis by frames
FrameStep=5; %Multiple intervals would be generated between ZIntByFrame(1) and ZIntByFrame(1), with step of 5 frames<<-----------Edited to adjust current experiment setting


ZIntCandidate=ZIntByFrame(1):FrameStep:ZIntByFrame(2); %Candidate ZInt 

SLMnum=10;  %%number of SLM in a Tseries<<---------------------------------------------------------------------------------------Edited to adjust current experiment setting
ZNum=SLMnum+1;%%1st Z is without MP stimuli; The rest were synchronized with MarkPoints.
 
SequenceN=5;%%Number of Tseries in total, real experiment would rotate in the Tseries generated.<<-------------------------------Edited to adjust current experiment setting
FunGroupN=3;%%Number of SLM target groups.<<-------------------------------------------------------------------------------------Edited to adjust current experiment setting
FunGroupID=1:FunGroupN;

clear SequenceTemp TSeries MarkTWhiskStim MarkTZeroPower

PreSLMframeMin=40; %Minal number of frames in the 1st Zseries with each Tseries;<<-----------------------------------------------Edited to adjust current experiment setting
SLMZeroStimPro=0.109;%Probability that SLM stimuli is with Zero Power ;<<--------------------------------------------------------Edited to adjust current experiment setting
SLMWithWhiskStimPro=0.5;%Probability that SLM stimuli is together with a whisk stimuli ;<<---------------------------------------Edited to adjust current experiment setting
                        %Currently, whisk stimuli probility is the same for zeroPower SLM group and non-zero Power SLM group. 

RepStimNumTh=2;         %maximal number of the same target groups continously stimulate.<<---------------------------------------Edited to adjust current experiment setting
FrameN=550;             %frame number of a single plane recorded in a file, not including the 1st frame synchronized with SLM.<<-Edited to adjust current experiment setting
MaxDiffNum=2;           %Maximal Population differnece between SLM across different target groups.  
WhiskRepMaxTh=3;        %Maximal number of repeated SLM together with whisker stimuli.


TotalStimN=SLMnum*SequenceN;
ZeroN=round(SLMZeroStimPro*TotalStimN);
WhiskN=round(SLMWithWhiskStimPro*TotalStimN);


[TotalSequence,CountID] = generate_stimulus_sequence(SLMnum*SequenceN, FunGroupID, RepStimNumTh,MaxDiffNum);
SeqZero=GenerateZeroSeq(TotalSequence,FunGroupID,ZeroN);
SeqWhisk=GenerateWiskFromSeqAndZero(TotalSequence,SeqZero,SLMWithWhiskStimPro,WhiskRepMaxTh);



TSeriesFunMat=reshape(TotalSequence,SLMnum,SequenceN);
SeqZeroMat=reshape(SeqZero,SLMnum,SequenceN);
SeqWhiskMat=reshape(SeqWhisk,SLMnum,SequenceN);




TotalFrame=FrameN+SLMnum;  %%Ensure the total frame has FrameN+SLMnum frames for each plane, noted that SLMnum frames will be deleted
TSeriesCount=1;

while TSeriesCount <= SequenceN
    FrameSequenceTemp=ZIntCandidate(randi(length(ZIntCandidate),SLMnum,1));
    SequenceFunIDTemp=TSeriesFunMat(:,TSeriesCount);
    SequenceZeroTemp=SeqZeroMat(:,TSeriesCount);
    SequenceWhiskTemp=SeqWhiskMat(:,TSeriesCount);
    Pre1stZ=TotalFrame-sum(FrameSequenceTemp);
    if Pre1stZ>=PreSLMframeMin
       TSeries(:,TSeriesCount)=[Pre1stZ;FrameSequenceTemp(:)];
       TSeriesXMLFileID(:,TSeriesCount)=[0;SequenceFunIDTemp(:)];
       TSeriesZeroXMLFile(:,TSeriesCount)=[0;SequenceZeroTemp(:)];    %1st Z is not combined with whiskerStim for sure
       TSeriesVolOutSyn(:,TSeriesCount)=[0;SequenceWhiskTemp(:)]; 
       TSeriesBrukerTBL{TSeriesCount}=[TSeries(:,TSeriesCount) TSeriesXMLFileID(:,TSeriesCount) TSeriesVolOutSyn(:,TSeriesCount)];
       
       TSeriesBruker(TSeriesCount).Reps = TSeries(:,TSeriesCount);
       TSeriesBruker(TSeriesCount).SynMP = TSeriesXMLFileID(:,TSeriesCount)>0;       
       TSeriesBruker(TSeriesCount).SynMPFunGroup = TSeriesXMLFileID(:,TSeriesCount);       
       TSeriesBruker(TSeriesCount).VolOut = TSeriesVolOutSyn(:,TSeriesCount);       
       TSeriesBruker(TSeriesCount).ZeroPower = TSeriesZeroXMLFile(:,TSeriesCount);       

        % Construct TSeriesBrukerTBL table
        % Construct TSeriesBrukerTBL table
        TSeriesBrukerTBL{TSeriesCount} = table(TSeries(:, TSeriesCount), ...
            TSeriesXMLFileID(:, TSeriesCount) > 0, ...
            TSeriesXMLFileID(:, TSeriesCount), ...
            TSeriesVolOutSyn(:, TSeriesCount), ...
            TSeriesZeroXMLFile(:, TSeriesCount), ...
            'VariableNames', {'Reps', 'SynMP', 'SynMPFunGroup', 'VolOut', 'PowerZero'});
       TSeriesCount=TSeriesCount+1;
    end
end
save([SaveFolder 'SpontBeh5T_Z' num2str(ZNum) 'Frame' num2str(FrameN) '-' date '.mat']);
% imagesc([TotalSequence SeqZero SeqWhisk])


