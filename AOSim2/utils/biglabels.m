function biglabels(blurbx,blurby,fsize)
  
% BIGLABELS: Make a Big Bold statement.
%
% usage:  biglabels(blurbx,blurby,font size)
%
% Johanan L. Codona, Steward Observatory, CAAO
% Jan. 9, 2003
% Sep 16, 2005  Font size control
  
if(nargin<3)
    fsize = 12;
end

  xlabel(blurbx,'FontSize',fsize,'FontWeight', 'bold'); 
  ylabel(blurby,'FontSize',fsize,'FontWeight', 'bold'); 
