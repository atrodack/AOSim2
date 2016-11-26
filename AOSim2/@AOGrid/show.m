function show(AOG)
% SHOW: Make an image of the AOGrid.
%
% 20090408: JLCodona AOSim2

AOG.center;

if(AOG.isX)
    AOG.plotC;
    title(AOG.describe,'FontSize',14);
    daspect([1 1 1]);
else
    [x,y] = AOG.coords;
    imagesc(x,y,AOG.ndex,[-3 0]);
    axis square;
    axis xy;
    title(AOG.describe,'FontSize',14);
    colorbar;
    daspect([1 1 1]);
end

end

