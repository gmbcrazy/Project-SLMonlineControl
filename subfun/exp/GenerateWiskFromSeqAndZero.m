function SeqWhisk=GenerateWiskFromSeqAndZero(TotalSequence,SeqZero,SLMWithWhiskStimPro)
SeqWhisk=zeros(size(TotalSequence));
ZeroI=find(SeqZero==1);
StimI=find(SeqZero==0);
WiskZeroN=round(length(ZeroI)*SLMWithWhiskStimPro);
WiskStimN=round(length(StimI)*SLMWithWhiskStimPro);

AddWZero=randperm(length(ZeroI),WiskZeroN);
AddWStim=randperm(length(StimI),WiskStimN);

SeqWhisk(ZeroI(AddWZero))=1;
SeqWhisk(StimI(AddWStim))=1;

