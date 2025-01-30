function quick_roi_Ftraces()
Ftraces=load('C:\Work\Analysis\07192022\TSeries-07192022-1023-006\suite2p\plane0\Fall.mat');
FtracesF=Ftraces.F;

for i=1:size(FtracesF,1)
    FtracesFtemp=FtracesF(i,:);
    FtracesFtempMax=max(FtracesFtemp);
    FtracesFtempN=FtracesFtemp/FtracesFtempMax;
    FtracesFNT(i,:)=FtracesFtempN;

end
sampleSize=1:5000;
fs=30;
ys=16.6667;
xS=sampleSize/30;
FtracesF_G1=FtracesFNT([83,10,125,67,115,207,9],:);
figure;
subplot(7,1,1,'replace')
plot(xS,FtracesF_G1(1,:),'r')
hold on
plot(16.6667,0:0.1:1,'|k')
hold on
subplot(7,1,2,'replace')
plot(xS,FtracesF_G1(2,:),'r')
hold on
plot(16.6667,0:0.1:1,'|k')

hold on
subplot(7,1,3,'replace')
plot(xS,FtracesF_G1(3,:),'r')
hold on
plot(16.6667,0:0.1:1,'|k')

hold on
subplot(7,1,4,'replace')
plot(xS,FtracesF_G1(4,:),'r')
hold on
plot(16.6667,0:0.1:1,'|k')

hold on
subplot(7,1,5,'replace')
plot(xS,FtracesF_G1(5,:),'r')
hold on
plot(16.6667,0:0.1:1,'|k')

hold on
subplot(7,1,6,'replace')
plot(xS,FtracesF_G1(6,:),'r')
hold on
plot(16.6667,0:0.1:1,'|k')
hold on
%%%
FtracesF_G2=FtracesFNT([163,48,166,32,50,121,164,291,26,5,61,11,249,30],:);
figure;
subplot(13,1,1)
plot(xS,FtracesF_G2(1,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,2)
plot(xS,FtracesF_G2(2,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,3)
plot(xS,FtracesF_G2(3,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,4)
plot(xS,FtracesF_G2(4,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,5)
plot(xS,FtracesF_G2(5,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,6)
plot(xS,FtracesF_G2(6,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,7)
plot(xS,FtracesF_G2(7,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,8)
plot(xS,FtracesF_G2(8,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,9)
plot(xS,FtracesF_G2(9,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,10)
plot(xS,FtracesF_G2(10,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,11)
plot(xS,FtracesF_G2(11,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,12)
plot(xS,FtracesF_G2(12,:),' c')
hold on
plot(16.7167,0:0.1:1,'|k')
subplot(13,1,13)
plot(xS,FtracesF_G2(13,:),'c')
hold on
plot(16.7167,0:0.1:1,'|k')
%%%
FtracesF_G3=FtracesFNT([28,75],:);
figure;
subplot(2,1,1)
plot(xS,FtracesF_G3(1,:),'g')
hold on
plot(16.7667,0:0.1:1,'|k')
subplot(2,1,2)
plot(xS,FtracesF_G3(2,:),'g')
hold on
plot(16.7667,0:0.1:1,'|k')
%%
FtracesF_G4=FtracesFNT([77,16,104,235,31],:);
figure;
subplot(5,1,1)
plot(xS,FtracesF_G4(1,:),'b')
hold on
plot(16.8167,0:0.1:1,'|k')
subplot(5,1,2)
plot(xS,FtracesF_G4(2,:),'b')
hold on
plot(16.8167,0:0.1:1,'|k')
subplot(5,1,3)
plot(xS,FtracesF_G4(3,:),'b')
hold on
plot(16.8167,0:0.1:1,'|k')
subplot(5,1,4)
plot(xS,FtracesF_G4(4,:),'b')
hold on
plot(16.8167,0:0.1:1,'|k')
subplot(5,1,5)
plot(xS,FtracesF_G4(5,:),'b')
hold on
plot(16.8167,0:0.1:1,'|k')
%%%
FtracesF_G5=FtracesFNT([2,44,141,249,96,207,52,39,149,43,28,59,33,99,102],:);
figure;
subplot(15,1,1)
plot(xS,FtracesF_G5(1,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,2)
plot(xS,FtracesF_G5(2,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,3)
plot(xS,FtracesF_G5(3,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,4)
plot(xS,FtracesF_G5(4,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,5)
plot(xS,FtracesF_G5(5,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,6)
plot(xS,FtracesF_G5(6,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,7)
plot(xS,FtracesF_G5(7,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,8)
plot(xS,FtracesF_G5(8,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,9)
plot(xS,FtracesF_G5(9,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,10)
plot(xS,FtracesF_G5(10,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,11)
plot(xS,FtracesF_G5(11,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,12)
plot(xS,FtracesF_G5(12,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,13)
plot(xS,FtracesF_G5(13,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,14)
plot(xS,FtracesF_G5(14,:),'m')
hold on
plot(16.8667,0:0.1:1,'|k')
subplot(15,1,15)
plot(xS,FtracesF_G5(15,:),'m' )
hold on
plot(16.8667,0:0.1:1,'|k')
%%
x;