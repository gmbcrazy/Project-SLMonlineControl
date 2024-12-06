function Data=AmpNormalize(Data,varargin)

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

if nargin==1
   LowT=nanmin(Data(:));
   HighT=nanmax(Data(:));
   NanVec=[0 0];
elseif nargin==2
   prcTh=varargin{1};
   temp=prctile(Data(:),prcTh);

   LowT=temp(1);
   HighT=temp(2);
   NanVec=[0 0];

elseif nargin==3
   prcTh=varargin{1};
   NanVec=varargin{2};

   temp=prctile(Data(:),prcTh);

   LowT=temp(1);
   HighT=temp(2);
else

end
if NanVec(1)==1
Data(Data<LowT)=nan;
else
Data(Data<LowT)=LowT;
end
if NanVec(2)==1
Data(Data>HighT)=nan;
Data(Data>HighT)=HighT;
end

Data=(Data-LowT)/(HighT-LowT);

