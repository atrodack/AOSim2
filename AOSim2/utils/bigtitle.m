function th=bigtitle(blurb,fsize)

% BIGTITLE: Make a Big Bold statement.
%
% usage: th=bigtitle(blurb)
%
% Johanan L. Codona, Steward Observatory, CAAO
% Jan. 9, 2003
% Sep 16, 2005  Font size control
% 20080828: JLC return a text handle.

if(nargin<2)
    fsize = 14;
end

th=title(blurb,'FontName','Arial','FontSize',fsize,'FontWeight', 'bold');

return;
