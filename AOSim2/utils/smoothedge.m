function values = smoothedge(x,epsilon)
  
% SMOOTHEDGE: Make a smoothed edge transition from 0 to 1.
%
% USAGE: values = smoothedge(x,epsilon)
%
% x: array of inputs
% epsilon: scale over which transition is smoothed.
% values: values!
%
% Johanan L. Codona, CAAO, Steward Observatory, UA
% Sept. 6, 2002
  
  
  values = 1-.5*erfc(x/epsilon);
    
