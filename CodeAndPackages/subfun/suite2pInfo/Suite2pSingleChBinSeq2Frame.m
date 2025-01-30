function Data = Suite2pSingleChBinSeq2Frame(BinFile, Ly, Lx, nPlanes, FrameID)

diffID=diff(FrameID);
Checking=find(diffID<=0, 1);
if ~isempty(Checking)
   error('FrameID must be in ascending order');
   Data=[];
   return;
end


    % Open the binary file for reading
    fid = fopen(BinFile);

    % Set the number of bytes for data format (assuming uint16 data)
    Bytes = 2;

    % Get the number of frames to read
    nFrame = length(FrameID);

    % Initialize the Data matrix
    Data = zeros(Ly, Lx, nFrame);
    Jump = Ly * Lx * nPlanes * (FrameID(1) - 1);
            % Calculate the jump to the start of the current frame
        % Move the file pointer to the start of the current frame

    fseek(fid, Jump * Bytes, 'bof');

    % Loop through each specified frame
    for i = 1:length(FrameID)
        % Read data for the current frame and store it in Data matrix
        Data(:,:,i) = fread(fid, [Ly, Lx], 'uint16');
    end

    % Close the file
    fclose(fid);
    Data=double(Data);

end

% find(diff(FrameID)>1)



