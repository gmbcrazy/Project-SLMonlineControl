function PVparam=BrukerTBLtoPVparm(TseriesTBL,nPlane)



PVparam.InterMPRepetition=TseriesTBL.Reps(:)';
frameRepetition=sum(PVparam.InterMPRepetition); %%Total repepitions of Z series in T series setting;
PVparam.maxFrame=nPlane*frameRepetition;
PVparam.BreakPointFrame=PVparam.InterMPRepetition(1:end-1)*nPlane;
% PVparam.InterMPFrame=[40 60 30 20]*nPlane;
PVparam.TrialMPSwitch=length(PVparam.InterMPRepetition)-1;
PVparam.nPlane=nPlane;

PVparam.MPFunGroup=TseriesTBL.SynMPFunGroup(2:end);   %1st Z is not synchronized with MP
PVparam.MPPowerZero=TseriesTBL.PowerZero(2:end);      %1st Z is not synchronized with MP
PVparam.TseriesTBL=TseriesTBL;



