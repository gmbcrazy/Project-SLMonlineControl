drive = 'F:\';%'C'
user = 'Lu';


dateOrigin=datestr(now,'mmddyyyy');

% dateOrigin=datestr(now,'01252024');

DataFolder='F:\LuSLMOnlineTest\';
DataFolder=[DataFolder dateOrigin '\'];
% mkdir(subfolder);

i = 0;
db = [];
i = i+1;
db(i).mouse_name        = '';
db(i).date              = dateOrigin;
db(i).expts             = 1;%%change by Hari [1:3,6] to [1]
db(i).nchannels         = 1;
db(i).gchannel          = 1; 
db(i).nplanes           = 3;
db(i).planesToProcess   = db(i).nplanes;

ops0.RootStorage=DataFolder;
% ops0.RootStorage = [drive ':\Data\' user '\'];
% ops0.RootStorage = '/Users/henrydalgleish/Desktop/';
