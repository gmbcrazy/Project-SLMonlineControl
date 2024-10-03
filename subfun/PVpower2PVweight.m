function PVweight=PVpower2PVweight(PVpower, refPVPower)

HigherThanMax=find(PVpower > refPVPower);
if ~isempty(HigherThanMax)
   disp('Warning, the maximal PVpower needed is higher than the fix PV power');

end
PVweight=PVpower./refPVPower;

PVweight=ceil(PVweight*100)/100;
PVweight(PVweight>1)=1;