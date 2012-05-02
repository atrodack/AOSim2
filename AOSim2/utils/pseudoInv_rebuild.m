function pMinv = pseudoInv_rebuild(U,S,V,sthresh)

% function pMinv = pseudoInv_rebuild(U,S,V,sthresh)
%
% Johanan L. Codona: 20070403
% Johanan L. Codona: 20070424

% [U,S,V] = svd(M);

error('deprecated.  use the pseudoInv_rebuild2 version.');

s = diag(S);
SKEEP = (s<sthresh);
% SKEEP = (s<sthresh);

% si = 1./s(SKEEP);
si = 1./s;

si(~SKEEP) = 0;  % This "conditions" the matrix.
% SI = diag(si);

% Now make it the right shape...
% [sz1,sz2] = size(S);
% if(sz2 > size(SI,1))
% 	SI = [SI;zeros(sz2-size(SI,1),size(SI,2))];
% end
% if(sz1 > size(SI,2))
% 	SI = [SI,zeros(size(SI,1),sz1-size(SI,2))];
% end

% whos V D DI U

pMinv = V(:,1:length(si))*diag(si)*U(:,1:length(si))';

return
