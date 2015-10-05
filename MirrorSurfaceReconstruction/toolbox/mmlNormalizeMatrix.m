function aNormalizedMatrix = mmlNormalizeMatrix(aMat)
%% This function is used to normalize each column of matrix 'aMat'.
% aMat should have dimension 3*n.

[rr cc] = size(aMat);
if rr > cc&&cc==3,
    aMat = aMat';
end

aNormalizedMatrix = aMat./repmat(sqrt(sum(aMat.^2,1)),3,1);

   