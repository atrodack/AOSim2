function POINT = pickXYpoint(pnts)

% POINT = pickPoint()
% 
% 20080925: JLCodona: DMC49 Day.
% 20111016: Made report XY coords rather than pixels.

if(nargin==0)
	pnts=1;
end

POINT = ginput(pnts);

return;


