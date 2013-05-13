function pMinv = pseudoInv_rebuildN(U,s,V,Nmodes)

% function pMinv = pseudoInv_rebuildN(U,s,V,Nmodes)
%
% Johanan L. Codona: 20070403
% Johanan L. Codona: 20070424
% Johanan L. Codona: 20090804

if(Nmodes<=length(s))
	SKEEP = (1:length(s))<=Nmodes;
else
	warning('You asked for more modes than are possible.');
	SKEEP = true(size(s));
end
	
si = 1./s(SKEEP);

pMinv = V(:,1:length(si))*diag(si)*U(:,1:length(si))';

return
