%Test LoadRegRefFile;


FileType=0;
RefFile=[];
[ops, refImg] = LoadRegRefFile(RefFile, FileType);
MultiMatrix3DHeatmap(refImg)


FileType=0;
RefFile='E:\TestRegRefFile\TSeries-11142024-0927-003\';
[ops, refImg] = LoadRegRefFile(RefFile, FileType);
MultiMatrix3DHeatmap(refImg)




FileType=1;
RefFile='E:\LuSLMOnlineTest\SL0777-Ai203\11142024\Data\suite2p\';
[ops, refImg] = LoadRegRefFile(RefFile, FileType);
MultiMatrix3DHeatmap(refImg)


FileType=2;
RefFile='F:\LuSLMOnlineTest\SL0242-Ai203\09042024\TSeries-09042024-0929-001.bin';

[ops, refImg] = LoadRegRefFile(RefFile, FileType, [512,512,3,60]);
MultiMatrix3DHeatmap(refImg)


fameNum=7500*3;
nplane=3;



fileIDraw = fopen(RefFile,'r');
fileIDreg = fopen(['F:\LuSLMOnlineTest\SL0242-Ai203\09042024\Reg-002.bin'], 'wb');
Ly= 512;
Lx= 512;
for i=1:fameNum
     frame = fread(fileIDraw, [Ly, Lx], 'uint16');
     plane=mod(i-1,nplane)+1;
     tic
     [regFrame,dv(i,:),cv] = return_offsets_phasecorr(single((frame)),ops{plane});
     fwrite(fileIDreg,regFrame, 'uint16');
     toc
end

fclose(fileIDraw)
fclose(fileIDreg)


figure;
plot(dv)
legend({'x shift','y shfit'})
ylabel('Pixels')
xlabel('Time in Frames')
