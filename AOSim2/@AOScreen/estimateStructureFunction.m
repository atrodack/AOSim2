function [SF,x] = estimateStructureFunction(S,NN,dim)
% [SF,x] = AOSCREEN.estimateStructureFunction(PIXEL_OFFSETS,[dim])
% 
% This AOScreen method estimates the PHASE structure function at lambdaRef
% by direct calculation of k^2 <(delta AOSCREEN)^2> with pixel offsets
% given by PIXEL_OFFSETS.  This is just in one direction, given by dim.  
% The result is the phase structure function at the specified pixel offsets
% and x is a vector of physical offsets.  
% Plot loglog(x,SF) to see the result.  The Fried length, r0, is where
% SF=6.88.

if(nargin<2)
    dim = 1;
end
if(dim~=2)
    dim = 1;
end

NMAX = S.size - 1;

NN(NN>NMAX(dim)) = []; % Can't have lags larger than the grid.
NN(NN<0) = []; % Negative lags are redundant.

SF = zeros(size(NN));

k2 = (2*pi/S.lambdaRef)^2;

for n=1:length(NN)
    if(dim==2)
        SF(n) = k2*mean(mean((S.grid_(:,1+NN(n):end)-S.grid_(:,1:end-NN(n))).^2));
    else
        SF(n) = k2*mean(mean((S.grid_(1+NN(n):end,:)-S.grid_(1:end-NN(n),:)).^2));
    end
end

SPACING = S.spacing;
x = NN*SPACING(dim);

