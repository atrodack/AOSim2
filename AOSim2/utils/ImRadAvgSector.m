function Rav = ImRadAvgSector(ARRAY,CENTER,RADII,dR,anglow,anghigh)

% Rav = ImRadAvgSector(ARRAY,CENTER,RADII,dR,anglow,anghigh)
%
% General purpose radial average tool for plain MATLAB arrays.
%
% Note that the sector limits are in degrees for convenience.
% 
% 20050920: Johanan L. Codona, Steward Observatory.

x = 1:size(ARRAY,2);
y = 1:size(ARRAY,1);

[X,Y] = meshgrid(x-CENTER(1),y-CENTER(2));

R = sqrt(X.^2 + Y.^2);
ANG = atan2(Y,X);

Rav = zeros(size(RADII));
d = dR/2;
% 
% ff = figure;
% colormap(gray);
% setAnim;

for nR=1:length(RADII)
    SELECT = R>=RADII(nR)-d/2 ...
        & R<=RADII(nR)+d/2 ...
        & ANG > anglow*pi/180 ...
        & ANG < anghigh*pi/180;

%     imagesc(SELECT); sqar; drawnow;
    
    Rav(nR) = mean(mean(ARRAY(SELECT)));
    fprintf('%d: %.1f %.1g\n', nR,RADII(nR),Rav(nR));
end

% close(ff);
