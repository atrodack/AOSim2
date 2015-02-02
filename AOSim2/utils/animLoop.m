function animLoop(dcube, nloops, RNG, figNum,secsDelay)

% function animLoop(dcube,[nloops],[RNG],[figNum],[secsDelay])
%
% 20050718: JLC: added support for complex data cubes.

if(nargin<4)
    figNum = gcf;
end

if(nargin<2)
    nloops = 10
end

if(nargin<3)
    MIN = min(min(min(dcube)));
    MAX = max(max(max(dcube)));
else
    MIN = RNG(1);
    MAX = RNG(2);
end

if(MIN == MAX)
    MIN = MIN - 1;
    MAX = MAX + 1;
end

figure(figNum);

set(gcf,'DoubleBuffer','on');
set(gcf,'BackingStore','on');

%colormap(gray(256).^.5);

[NX NY NT] = size(dcube);

for iloop=1:nloops
    for iframe=1:NT
        if(isreal(dcube))
            imagesc(squeeze(dcube(:,:,iframe)),[MIN MAX]);
            colorbar;
        else
            %plotCAmpl(squeeze(dcube(:,:,iframe)),.5,[MIN MAX]);
            plotCAmpl(squeeze(dcube(:,:,iframe)),1/2);
        end
        %daspect([1 1 1]);
        title(sprintf('frame %d',iframe));
        
        drawnow;
        if(nargin>4)
            pause(secsDelay);
        end
    end
end
