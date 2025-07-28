
function FileList=SearchFile(path,FileKeyWord, varargin)


FileList=dir([path '*' FileKeyWord]);

if nargin==3
   LookingtoSubFolder=varargin{1};
else
   LookingtoSubFolder=1;

end

FileList1=dir(path);

for i=3:length(FileList1)
    if FileList1(i).isdir&&LookingtoSubFolder==1
       FileListTemp=SearchFile([FileList1(i).folder '\' FileList1(i).name '\'],FileKeyWord);
     
       FileList=[FileList;FileListTemp];
    end
end
