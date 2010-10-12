function show(AOG)
% SHOW: Make an image of the AOGrid.
%
% 20090408: JLCodona AOSim2

AOG.center;
[x,y] = AOG.coords;

if(AOG.isX)
%     if(isreal(AOG.grid_))
%         imagesc(x,y,AOG.grid());
%         axis square;
%         title([class(AOG) ' ' AOG.name ': axis:' AOG.axis ' domain:' AOG.domain ],'FontSize',14);
%         colorbar;
%         drawnow;
%     else
% 		%         plotCAmpl(AOG.grid(),1/2);
% 		%         axis square;
%         % colorbar;
%         drawnow;
% 	end
	
	AOG.plotC;
	title([class(AOG) ' ' AOG.name ': axis:' AOG.axis ' domain:' AOG.domain ],'FontSize',14);
    daspect([1 1 1]);
    %drawnow;
else
    imagesc(x,y,AOG.ndex,[-3 0]);
    % plotCAmpl(AOG.grid_,1/4);
    axis square;
    axis xy;
    title([class(AOG) ' ' AOG.name ': axis:' AOG.axis ' domain:' AOG.domain ],'FontSize',14);
    colorbar;
    daspect([1 1 1]);
    %drawnow;
end

end

