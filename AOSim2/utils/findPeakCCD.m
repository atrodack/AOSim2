function PEAK = findPeakCCD(ccd)

% PEAK = findPeakCCD(ccd)
% 
% JLCodona: 9/25/2007

[mx1,loc1] = max(ccd);
[mx2,loc2] = max(mx1);

PEAK = [loc1(loc2),loc2];
