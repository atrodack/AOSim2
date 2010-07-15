function AOA = make(AOA)

% MAKE: Render the AOSegment grid.
% USAGE: S = make(S) or S.make or make(S)
%    >>> works on itself.
%
% This make is now a switch function based on pupilDef version.
% version 1 is the original definition matrix, see help AOSegment/make1.
% version 2 is TBD.
% 
% 20090412: JLCodona: Made into a version switch for AOSim2.

if(AOA.version == 1)
	AOA = make1(AOA);
else
	error('max AOSegment version is 1');
end

