function Data = Suite2pSingleChBin2Frame(BinFile, Ly, Lx, nPlanes, FrameID)

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
    if nPlanes==1
    Data = zeros(Ly, Lx, nFrame);
    elseif nPlanes>1
    else

    end

    if nPlanes==1
    % Loop through each specified frame
       Data = zeros(Ly, Lx, nFrame);
    for i = 1:length(FrameID)
        % Calculate the jump to the start of the current frame
        Jump = Ly * Lx * nPlanes * (FrameID(i) - 1);
        % Move the file pointer to the start of the current frame
        fseek(fid, Jump * Bytes, 'bof');
        ibase=(i-1)*nPlanes;
        % Read data for the current frame and store it in Data matrix
        Data(:,:,i) = fread(fid, [Ly, Lx], 'uint16');

    end
    elseif nPlanes>1
        Data = zeros(Ly, Lx, nFrame, nPlanes);

        for i = 1:length(FrameID)
        % Calculate the jump to the start of the current frame
          Jump = Ly * Lx * nPlanes * (FrameID(i) - 1);
        
        % Move the file pointer to the start of the current frame
        fseek(fid, Jump * Bytes, 'bof');
        % Read data for the current frame and store it in Data matrix
             for iplane=1:nPlanes
                  Data(:,:,i,iplane) = fread(fid, [Ly, Lx], 'uint16');
             end
        end
    else
         Data=[];
         return;
    end
    % Close the file
    fclose(fid);
end


% find(diff(FrameID)>1)



