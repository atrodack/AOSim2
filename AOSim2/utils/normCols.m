function NORMS = normCols(vec)

%% NORMS = normCols(vec)
% 20080830: JLCodona

%% Only L2 norms for now.

NORMS = sqrt(sum(abs(vec).^2,1));

return;
