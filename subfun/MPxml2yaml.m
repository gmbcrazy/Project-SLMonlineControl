function laserPower=MPxml2yaml(MPxmlFile)
% This function reads information from an XML file and converts it into a YAML structure.

% Open the XML file for reading
fid=fopen(MPxmlFile);

% % Find and extract the maximum voltage parameters for X and Y axes
% output=FindKeywords(fid,'maxVoltage');
% str=fgetl(fid);
% pattern = '"XAxis" value="(-?[\d.]+)"';
% tokens = regexp(str, pattern, 'tokens');
% Xmax=str2num(tokens{1}{1});
% pattern = '"YAxis" value="(-?[\d.]+)"';
% str=fgetl(fid);
% tokens = regexp(str, pattern, 'tokens');
% Ymax=str2num(tokens{1}{1});
% 
% % Find and extract the minimum voltage parameters for X and Y axes
% fseek(fid,0,'bof');
% output=FindKeywords(fid,'minVoltage');
% str=fgetl(fid);
% pattern = '"XAxis" value="(-?[\d.]+)"';
% tokens = regexp(str, pattern, 'tokens');
% Xmin=str2num(tokens{1}{1});
% pattern = '"YAxis" value="(-?[\d.]+)"';
% str=fgetl(fid);
% tokens = regexp(str, pattern, 'tokens');
% Ymin=str2num(tokens{1}{1});
% 
% % Define scan amplitudes for X and Y axes
% ScanAmp_X= [Xmin Xmax];
% ScanAmp_Y= [Ymin Ymax];

% Find and extract the 2P imaging laser power
fseek(fid,0,'bof');
str=FindKeywords(fid,'UncagingLaser="Satsuma" UncagingLaserPower=');
% str=fgetl(fid);
pattern = 'UncagingLaser="Satsuma" UncagingLaserPower="(-?[\d.]+)"*';
tokens = regexp(str, pattern, 'tokens');
laserPower=str2num(tokens{1}{1});

% % Find and extract the number of lines per frame and pixels per line
% fseek(fid,0,'bof');
% str=FindKeywords(fid,'key="linesPerFrame" value');
% pattern = 'key="linesPerFrame" value="(-?[\d.]+)"';
% tokens = regexp(str, pattern, 'tokens');
% Xnum=str2num(tokens{1}{1});

% fseek(fid,0,'bof');
% str=FindKeywords(fid,'key="pixelsPerLine" value');
% pattern = 'key="pixelsPerLine" value="(-?[\d.]+)"';
% tokens = regexp(str, pattern, 'tokens');
% Ynum=str2num(tokens{1}{1});
% 
% % Find and extract the optical zoom, microns per pixel for X and Y axes, and bit depth
% fseek(fid,0,'bof');
% str=FindKeywords(fid,'opticalZoom');
% pattern = '"opticalZoom" value="(-?[\d.]+)"';
% tokens = regexp(str, pattern, 'tokens');
% ZoomX=str2num(tokens{1}{1});
% 
% fseek(fid,0,'bof');
% str=FindKeywords(fid,'micronsPerPixel');
% str=fgetl(fid);
% pattern = '"XAxis" value="(-?[\d.]+)"';
% tokens = regexp(str, pattern, 'tokens');
% microsPerPixelx=str2num(tokens{1}{1});
% pattern = '"YAxis" value="(-?[\d.]+)"';
% str=fgetl(fid);
% tokens = regexp(str, pattern, 'tokens');
% microsPerPixely=str2num(tokens{1}{1});




% Close the file
fclose(fid);

% Store the extracted information into a YAML structure
% yaml.SLM_Pixels_X=Xnum;
% yaml.SLM_Pixels_Y=Ynum;
% yaml.SLM_BitDepth=bitDepth;
% yaml.ScanAmp_X=ScanAmp_X;
% yaml.ScanAmp_Y=ScanAmp_Y;
% yaml.scan_Z=Zmicro;
% yaml.scan_ZInd=ZInd;
% yaml.scan_Zdescription=Zdescription;
% yaml.Zdepth_ETL=Zdepth_ETL;
% yaml.Zdepth_laserPower=laserPower;
% 
% yaml.FOVsize_OpticalZoom=ZoomX;
% yaml.FOVsize_PX=Xnum;
% yaml.FOVsize_PY=Ynum;
% yaml.FOVsize_UM=[microsPerPixelx*Xnum microsPerPixely*Ynum];
% yaml.umPerlPixelX=microsPerPixelx;
% yaml.umPerlPixely=microsPerPixely;
% 

function output=FindKeywords(filePtr,Keywords)
% This function searches for specific keywords in a file and returns the first line containing the keyword.

% Initialize output to -1, indicating no match found initially
output=-1;

% Loop through the file until the end is reached
while (feof(filePtr) ~= 1)
    % Read a line from the file
    junk=fgetl(filePtr);	% this will make the function execute much faster but should be optional
    
    % Check if the line contains the specified keyword
    if strfind(junk,Keywords)
        % If the keyword is found, store the line in the output and break out of the loop
        output=junk;
        break
    end
end


