function [Pos3D, Pos3DRaw, CaData, CaDataPlane, stat, yaml] = Suite2pMultiPlaneROIToXYZ(LoadPath, varargin)
    % Suite2pMultiPlaneROIToXYZ - Processes multi-plane ROI data from Suite2P
    % and converts it into 3D positions.
    %
    % Inputs:
    %    LoadPath - Path to the directory containing the Suite2P data.
    %    varargin - Optional, FileID as a string or number to specify the file.
    %
    % Outputs:
    %    Pos3D - 3D positions of cells identified by Suite2P.
    %    Pos3DRaw - Raw 3D positions of all cells before filtering.
    %    CaData - Calcium imaging data from Suite2P.
    %    CaDataPlane - Calcium imaging data separated by imaging planes.
    %    stat - Filtered cell statistics from Suite2P.
    %    yaml - Metadata extracted from the XML file.

    % Find the folder containing Suite2P data within the specified LoadPath
    Suite2PPath = findAllFolders(LoadPath, 'suite2p');
    
    % If only one Suite2P folder is found, use it
    if length(Suite2PPath) == 1
        Suite2PPath = Suite2PPath{1};
    elseif length(Suite2PPath) > 1
        for isuite=1:length(Suite2PPath)
            tempL(isuite)=length(Suite2PPath{isuite});
        end
        [~,i1]=min(tempL);
        Suite2PPath=Suite2PPath{i1};
    else
        disp(['No suite2p processed folder is found in ' LoadPath])
        return
    end

    % Determine the FileID based on input arguments
    if nargin == 1
        FileID = '001';
    else
        FileID = num2str(varargin{1});
        if length(FileID) == 1
            FileID = ['00' FileID];
        elseif length(FileID) == 2
            FileID = ['0' FileID];
        else
            % Do nothing if FileID is already in correct format
        end
    end

    % Locate the XML file that contains recording information
    xmlFile = dir([LoadPath '*TSeries*-' FileID '\TSeries-*-' FileID '.xml']);
    if isempty(xmlFile)
        % If no XML file is found, display a message and return empty outputs
        disp('Check Path and File Name, No .xml file is detected for recording information')
        SavePath = [];
        Pos3D = [];
        CaData = [];
        stat = [];
        yaml = [];
        confSet = [];
        CaDataPlane = [];
        return;
    end

    % Convert the XML file to YAML format
    xmlFile = [xmlFile.folder '\' xmlFile.name];
    yaml = xml2yaml(xmlFile);

    % Load the combined Suite2P data
    Suite2Temp = [Suite2PPath '\combined\Fall.mat'];
    CaData = load(Suite2Temp);

    % Get the number of planes from the Suite2P data
    PlaneN = double(CaData.ops.nplanes);

    % Get the list of plane folders within the Suite2P directory
    planFolder = dir([Suite2PPath 'plane*']);
    clear CaDataPlane;

    % Check if the number of plane folders matches the expected number of planes
    if length(planFolder) ~= PlaneN
        disp('Plane folder number does Not match plane #');
    else
        % Load the data for each plane
        for i = 1:PlaneN
            tempFolder = [planFolder(i).folder '\' planFolder(i).name '\Fall.mat'];
            if i == 1
                CaDataPlane(i) = load(tempFolder);
            else
                CaDataPlane = concatenateStructs(CaDataPlane, load(tempFolder));
            end
        end
    end

    % Get the dimensions of the imaging planes
    Lx = CaDataPlane(1).ops.Lx;
    Ly = CaDataPlane(1).ops.Ly;
    SLMtargetColor = [0.1 0.8 0.2];

    % Initialize variables for storing data
    CaDataNew = CaData;
    statTemp = {};
    Shown = [];
    CellPlaneID = [];

    % Collect statistics for each plane
    for i = 1:PlaneN
        statTemp = [statTemp CaDataPlane(i).stat];
        CellPlaneID = [CellPlaneID; zeros(length(CaDataPlane(i).stat), 1) + i];
        CaData.PlaneMeanImg(:, :, i) = CaDataPlane(i).ops.meanImg;
    end

    % Store the raw statistics and plane IDs
    statRaw = statTemp;
    CellPlaneIDRaw = CellPlaneID;

    % Filter out non-cell entries
    CellPlaneID(CaDataNew.iscell(:, 1) == 0) = [];
    stat = CaDataNew.stat;
    stat(CaDataNew.iscell(:, 1) == 0) = [];
    statTemp(CaDataNew.iscell(:, 1) == 0) = [];

    % Extract XY pixel positions for each cell
    xyPix = [];
    for iCell = 1:length(statTemp)
        xyPix(iCell, :) = statTemp{iCell}.med(end:-1:1);
    end

    % Store the processed statistics and positions
    CaData.statCell = statTemp;
    yaml.Zdepth_ETL = unique(round(yaml.Zdepth_ETL));
    CaData.CellPlaneID = CellPlaneID;
    Pos3D = [];
    Zmicro = yaml.Zdepth_ETL(CellPlaneID);
    Pos3D = [xyPix Zmicro(:)];

    % Extract raw XY pixel positions
    xyPixRaw = [];
    for iCell = 1:length(statRaw)
        xyPixRaw(iCell, :) = statRaw{iCell}.med(end:-1:1);
    end

    % Store the raw positions and plane IDs
    CaData.CellPlaneIDRaw = CellPlaneIDRaw;
    Pos3DRaw = [];
    Zmicro = yaml.Zdepth_ETL(CellPlaneIDRaw);
    Pos3DRaw = [xyPixRaw Zmicro(:)];
end
