function Data=AmpNormalizeDim(Data,NormDim,prcTh)

%% AmpNormalizeDim - Normalize data by amplitude along a specified dimension
% AmpNormalize - Normalize Data by Amplitude
%   Data = AmpNormalize(Data)
%   Data = AmpNormalize(Data, prcTh)
%   Data = AmpNormalize(Data, prcTh, NanVec)
%
% Input:
%   - Data: Input data to be normalized.
%   - varargin{1}: Percentile threshold for defining low and high limits for amplitude.
%   - varargin{2}: Optional, NanVec controls NaN replacement. NanVec(1)=1 sets values below lower limit to NaN,
%                  NanVec(2)=1 sets values above higher limit to NaN. Default is [0 0].
%
% Output:
%   - Data: Normalized data.
%
% Example Usage:
%   1. Data = AmpNormalize(Data)
%   2. Data = AmpNormalize(Data, [10 90])
%   3. Data = AmpNormalize(Data, [10 90], [1 0])

SubData = GetDataDim(Data,NormDim);

for iDim=1:length(SubData)
    SubData{iDim}=AmpNormalize(SubData{iDim},prcTh);
end

Data = ReconstructFromDataDim(SubData, NormDim);

end






function Data = GetDataDim(X, DimM)
    % GetDataDim - Extracts each submatrix along a specified dimension DimM
    %
    % Syntax: Data = GetDataDim(X, DimM)
    %
    % Inputs:
    %    X - Input multi-dimensional matrix (can be 4D or any N-D matrix)
    %    DimM - The dimension along which to extract submatrices
    %
    % Outputs:
    %    Data - A cell array where each element is a submatrix along DimM
    
    % Get the size of the input matrix
    matrixSize = size(X);
    
    % Determine the number of submatrices along the specified dimension
    numSubMatrices = matrixSize(DimM);
    
    % Initialize the output cell array
    Data = cell(1, numSubMatrices);
    
    % Prepare an indexing cell array
    idx = repmat({':'}, 1, ndims(X));
    
    % Loop through each submatrix along the specified dimension
    for i = 1:numSubMatrices
        idx{DimM} = i;             % Set the index for the DimM dimension
        Data{i} = squeeze(X(idx{:}));  % Extract the submatrix and remove singleton dimensions
    end
end





function X = ReconstructFromDataDim(Data, DimM)
    % ReconstructFromDataDim - Reconstructs a multi-dimensional matrix X
    % from submatrices stored in a cell array Data, arranged along dimension DimM.
    %
    % Syntax: X = ReconstructFromDataDim(Data, DimM)
    %
    % Inputs:
    %    Data - Cell array containing submatrices
    %    DimM - The dimension along which submatrices are to be arranged
    %
    % Outputs:
    %    X - Reconstructed multi-dimensional matrix
    
    % Determine the size of the submatrices
    submatrixSize = size(Data{1});
    
    % Determine the total number of submatrices along DimM
    numSubMatrices = length(Data);
    
    % Prepare the full size of the resulting matrix X
    fullSize = submatrixSize;
    fullSize(DimM) = numSubMatrices;  % Set the size along DimM
    
    % Initialize the output matrix X with zeros
    X = zeros(fullSize);
    
    % Prepare an indexing cell array
    idx = repmat({':'}, 1, length(fullSize));
    
    % Loop through each submatrix and place it in the correct position in X
    for i = 1:numSubMatrices
        idx{DimM} = i;  % Set the index for the DimM dimension
        X(idx{:}) = Data{i};  % Assign the submatrix to the appropriate slice in X
    end
end

