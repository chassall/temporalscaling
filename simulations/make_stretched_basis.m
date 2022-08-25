function stretchedBasis = make_stretched_basis(numRows,numCols,stretchMethod)
%MAKESTRETCHEDBASIS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    stretchMethod = 'box';
end

thisEye = eye(numCols);
stretchedBasis = imresize(thisEye,[numRows,numCols],stretchMethod);

end

