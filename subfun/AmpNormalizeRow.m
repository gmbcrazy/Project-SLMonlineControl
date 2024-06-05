function Data=AmpNormalizeRow(Data,varargin)

%% AmpNormalize - Normalize Data by Amplitude
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


for iRow=1:size(Data,1)
    Data(iRow,:)=AmpNormalize(Data(iRow,:),varargin);
end
