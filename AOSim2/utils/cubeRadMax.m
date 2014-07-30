function CUBE_ravg = cubeRadMax(CUBE,CENTER,RADII,WIDTH)

% CUBE_ravg = cubeRadAvg(CUBE,CENTER,RADII,WIDTH)


CUBE_ravg = zeros(length(RADII),size(CUBE,3));

for n=1:size(CUBE,3)
	modprint(n,10);
	Ravg = ImRadAvg(CUBE(:,:,n),CENTER,RADII,WIDTH);
	%plot(RADII,Ravg); drawnow;
	CUBE_ravg(:,n) = Ravg;
end

fprintf('\n');



