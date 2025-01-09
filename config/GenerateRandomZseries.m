clear all
SampleFre=6.9;  %Sampling rate in Bruker
SaveFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\';
SaveFolder='C:\Users\User\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\';

ZIntBySecond=[6;10]; %Interval for a Zseries in Tsereis by seconds
ZIntByFrame=round(ZIntBySecond*SampleFre);%%Interval for a Zseries in Tsereis by frames
FrameStep=5; %Multiple intervals would be generated between ZIntByFrame(1) and ZIntByFrame(1), with step of 5 frames


ZIntCandidate=ZIntByFrame(1):FrameStep:ZIntByFrame(2); %Candidate ZInt 

SLMnum=10;  %%number of SLM in a Tseries
ZNum=SLMnum+1;%%1st Z is without MP stimuli; The rest were synchronized with MarkPoints.
 
SequenceN=5;
FunGroupN=3;
FunGroupID=1:FunGroupN;

clear SequenceTemp TSeries MarkTWhiskStim MarkTZeroPower

PreSLMframeMin=40; %At least 40 frames before SLM;
SLMWithWhiskStimPro=0.5;
SLMZeroStimPro=0.109;

RepStimNumTh=3;

% sequence = generate_stimulus_sequence(SLMnum, FunGroupID, RepStimNumTh)
TotalStimN=SLMnum*SequenceN;
ZeroN=round(SLMZeroStimPro*TotalStimN);
WhiskN=round(SLMWithWhiskStimPro*TotalStimN);


[TotalSequence,CountID] = generate_stimulus_sequence(SLMnum*SequenceN, FunGroupID, RepStimNumTh,2);
SeqZero=GenerateZeroSeq(TotalSequence,FunGroupID,ZeroN);
SeqWhisk=GenerateWiskFromSeqAndZero(TotalSequence,SeqZero,SLMWithWhiskStimPro);



TSeriesFunMat=reshape(TotalSequence,SLMnum,SequenceN);
SeqZeroMat=reshape(SeqZero,SLMnum,SequenceN);
SeqWhiskMat=reshape(SeqWhisk,SLMnum,SequenceN);

imagesc([TotalSequence SeqZero SeqWhisk])


FrameN=550;
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

       TSeriesCount=TSeriesCount+1;

    end
end

save([SaveFolder 'Z' num2str(ZNum) 'Frame' num2str(FrameN) '-' date '.mat']);
