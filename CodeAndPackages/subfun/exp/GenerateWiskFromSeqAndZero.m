function SeqWhisk=GenerateWiskFromSeqAndZero(TotalSequence,SeqZero,SLMWithWhiskStimPro,WhiskRepMaxTh)
SeqWhisk=zeros(size(TotalSequence));
ZeroI=find(SeqZero==1);
StimI=find(SeqZero==0);
WiskZeroN=round(length(ZeroI)*SLMWithWhiskStimPro);
WiskStimN=round(length(StimI)*SLMWithWhiskStimPro);



while 1
    SeqWhisk=zeros(size(TotalSequence));

    AddWZero=randperm(length(ZeroI),WiskZeroN);
    AddWStim=randperm(length(StimI),WiskStimN);

    SeqWhisk(ZeroI(AddWZero))=1;
    SeqWhisk(StimI(AddWStim))=1;

    SameSeqPeriod = MarkToPeriod(SeqWhisk);
    if isempty(SameSeqPeriod)
       break
    end
    SameSeqL = SameSeqPeriod(2, :) - SameSeqPeriod(1, :) + 1;
    if max(SameSeqL)<=WhiskRepMaxTh
       break;
    end
end
