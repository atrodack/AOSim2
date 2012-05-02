function [pMinv,U,S,V,SI] = pseudoInv(M,sthresh)

% function [pMinv,U,S,V,SI] = pseudoInv(M,sthresh)
%
% Johanan L. Codona: 20070403
% 20081018 JLC: made thresh refer to max rather than abs.

[U,S,V] = svd(M);

si = 1./diag(S);
si(diag(S)<sthresh*S(1,1)) = 0;  % This "conditions" the matrix.
SI = diag(si);

% Now make it the right shape...
[sz1,sz2] = size(S);
if(sz2 > size(SI,1))
	SI = [SI;zeros(sz2-size(SI,1),size(SI,2))];
end
if(sz1 > size(SI,2))
	SI = [SI,zeros(size(SI,1),sz1-size(SI,2))];
end

% whos V D DI U

pMinv = V*SI*U';

return
