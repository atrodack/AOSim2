clear all;
clc;
close all;

nxy = 165;
spacing = 3;
nActs = 32;


G = AOGrid(nxy);
G.center;

[X,Y] = G.COORDS;
centerpix = G.AXIS_PIXEL;
cornerpix = [centerpix(1)-(nActs/2)*spacing,centerpix(2)-(nActs/2)*spacing;centerpix(1)+(nActs/2)*spacing,centerpix(2)-(nActs/2)*spacing;centerpix(1)-(nActs/2)*spacing,centerpix(2)+(nActs/2)*spacing;centerpix(1)+(nActs/2)*spacing,centerpix(2)+(nActs/2)*spacing]
edgepix = cornerpix(1,1):spacing:cornerpix(2,1);

if length(edgepix) > nActs
    removepix = length(edgepix) - nActs;
    edgepix = edgepix(1:end-removepix);
end

pixelpicx = zeros(nActs);
pixelpicy = zeros(nActs);

for ii = edgepix
    for jj = edgepix
        pixelpicx(ii-cornerpix(1,1)+1,jj-cornerpix(1,2)+1) = X(ii,jj);
        pixelpicy(ii-cornerpix(1,1)+1,jj-cornerpix(1,2)+1) = Y(ii,jj);
    end
end

G.show;
hold on
plot(pixelpicx,pixelpicy,'.');
hold off