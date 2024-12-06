function [corrPerm,corrTh]=ROIMatchPermutation(Mapping,ROIInfo,Palpha,PermutationParam)


FixedCellID=ROIInfo.FixedCellID;
MovingCellID=ROIInfo.MovingCellID;
FixRoiNAll=ROIInfo.FixRoiNAll;
MovingRoiNAll=ROIInfo.MovingRoiNAll;



if PermutationParam<=1

%% Extract distribution of all non-matching cell pair
% Loop through each Fixed cell to calculate correlation with each non-matching Moving cells
ii=1;
corrPerm=[];
for i=1:length(FixedCellID)
    [ismem,ismemI]=ismember(FixedCellID(i),Mapping(:,1));
    if ismem
       tempI=setdiff(MovingCellID,Mapping(ismemI,2));
    else
       tempI=MovingCellID;
    end
    temp1=FixRoiNAll(:,:,i);
    temp2=MovingRoiNAll(:,:,tempI);
    temp2=reshape(temp2,size(temp2,1)*size(temp2,2),length(tempI));
    RTemp=corr(temp1(:),temp2);
    corrPerm=[corrPerm;RTemp(:)];
end

else
    ShuffleNum=PermutationParam;
%% Extract distribution of non-matching cell pair using permutation method

% Generate random indices for Fixed and Moving cells for each permutation
P1=randi(length(FixedCellID),1,ShuffleNum);
P2=randi(length(MovingCellID),1,ShuffleNum);
corrPerm=[];
for i=1:ShuffleNum
    temp1=FixRoiNAll(:,:,P1(i));
    temp2=MovingRoiNAll(:,:,P2(i));
    corrPerm(i,1)=corr(temp1(:),temp2(:));
end


   
end


corrTh=prctile(corrPerm,(1-Palpha)*100);

