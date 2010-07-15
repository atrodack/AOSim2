function A = setDomain(A,dom)
  
% SETDOMAIN: Set the domain for an AOGrid object. (NOTE: Does no FFTs.)
% 
% USAGE:
%  A = setDomain(A,dom)
%
% INPUTS:
%    A: An AOGrid object.
%    dom: either 'k' or 'x'.
%
% NOTE: This will probably screw things up unless you are starting clean.
% 
% Written by: Johanan L. Codona, Steward Observatory, CAAO
% July 12, 2002: First writing.
% 20090413 JLCodona: tweaked for AOSim2.

  switch dom
   case 'x'
	   A.domain_ = AOGrid.DOMAIN_SPACE;
   case 'k'
    A.domain_ = AOGrid.DOMAIN_FREQ;
    A.axis_ = AOGrid.AXIS_CORNER;
   otherwise
    error(['unknown domain: ' dom]);
  end
  
