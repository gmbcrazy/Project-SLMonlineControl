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
RefFile='E:\TestRegRefFile\TSeries-11142024-0927-003.bin';

[ops, refImg] = LoadRegRefFile(RefFile, FileType, [512,512,3,50]);
MultiMatrix3DHeatmap(refImg)
