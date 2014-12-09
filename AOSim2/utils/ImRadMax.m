function Rav = ImRadMax(ARRAY,CENTER,RADII,dR)

% Rmax = ImRadMax(ARRAY,CENTER,RADII,dR)
%
% Like radial average tool, but with max.
%
% 20140305: Johanan L. Codona, Steward Observatory.


x = 1:size(ARRAY,2);
y = 1:size(ARRAY,1);

[X,Y] = meshgrid(x-CENTER(2),y-CENTER(1));

R = sqrt(X.^2 + Y.^2);

Rav = zeros(size(RADII));

for nR=1:length(RADII)
	SEL = ...
		(R>=RADII(nR)-dR/2) & ...
		(R<=RADII(nR)+dR/2);
	Rav(nR) = max(ARRAY(vec1(SEL)));
    %imagesc(SEL); sqar; drawnow;
    %Rav(nR) = mean(mean(ARRAY(R>=RADII(nR)-dR/2&R<=RADII(nR)+dR/2)));
    %fprintf('%d: %.1f %.1g\n', nR,RADII(nR),Rav(nR));
end
