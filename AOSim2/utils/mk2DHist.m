function [HISTIMG,XBINS,YBINS] = mk2DHist(DATAx,DATAy,pixelx,pixely)

% [HISTIMG,XBINS,YBINS] = mk2DHist(DATAx,DATAy,pixelx,pixely)
% 
% JLCodona 20120322

YLIMS = [min(DATAy) max(DATAy)];


NUMx = round(squeeze(DATAx)/pixelx);
XBINS = (min(NUMx):max(NUMx))*pixelx;

NUMx = NUMx - min(NUMx) + 1;
Nx = max(NUMx);

NUMy = round(squeeze(DATAy)/pixely);
YBINS = (min(NUMy):max(NUMy))*pixely;
NUMy = NUMy - min(NUMy) + 1;
Ny = max(NUMy);

if(Nx>1000 || Ny>1000)
    fprintf('Are you sure? The dims would be %dx%d\n',Nx,Ny);
    return;
end

HISTIMG = zeros(Ny,Nx);

for n=1:length(NUMx)
    HISTIMG(NUMy(n),NUMx(n)) = HISTIMG(NUMy(n),NUMx(n)) + 1;
end


