function POINT = pickPoint(pnts)

% POINT = pickPoint()
% 
% 20080925: JLCodona: DMC49 Day.

if(nargin==0)
	pnts=1;
end

POINT = round(fliplr(ginput(pnts)));

return;


