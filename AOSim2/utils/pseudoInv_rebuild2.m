function pMinv = pseudoInv_rebuild2(U,s,V,sthresh)

% function pMinv = pseudoInv_rebuild2(U,s,V,sthresh)
%
% Johanan L. Codona: 20070403
% Johanan L. Codona: 20070424

% [U,S,V] = svd(M);

% s = diag(S);
SKEEP = (s>=sthresh);
fprintf('Keeping %d modes.\n',sum(SKEEP));
si = 1./s(SKEEP);
% si(SCUT) = 0;  % This "conditions" the matrix.
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

% error('broken');
pMinv = V(:,1:length(si))*diag(si)*U(:,1:length(si))';

return
