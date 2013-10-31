function drawCircles(RADII,CENTER,scale,fmt)

%function drawCircles(RADII,CENTER,scale,[fmt])

if(nargin<4)
fmt = '--w';
end

TANGS = (0:90)/45*pi;

hold on
for RAD=RADII
    plot(CENTER(2)+RAD*scale*cos(TANGS),CENTER(1)+RAD*scale*sin(TANGS),fmt);
end
hold off
