function [array] = normalize(array)
  
% normalize: Normalize an array so its max value is unity.
% USAGE: [normed] = normalize(array)
%
% Johanan L. Codona, Steward Observatory, CAAO
% August 27, 2002 - Got tired of doing this manually!

% mx = max(abs(array(:)));
mx = max(array(:));
if(mx~=0)
	array = array/mx;
end
