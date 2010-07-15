function NORMS = normRows(vec)

%% NORMS = normRows(vec)
% 20080830: JLCodona

%% Only L2 norms for now.

NORMS = sqrt(sum(abs(vec).^2,2));

return;
