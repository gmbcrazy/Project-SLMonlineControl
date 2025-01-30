
fileList=dir('F:\LuSLMOnlineTest\04022024\DataDelete1stFrame\TSeries*')

for i=1:length(fileList)
    if fileList(i).isdir
       a=dir([fileList(i).folder '\' fileList(i).name '\*Cycle00001_Ch2*.ome.tif']);
       if ~isempty(a)
          delete([a.folder '\' a.name]);
       end
    end
end