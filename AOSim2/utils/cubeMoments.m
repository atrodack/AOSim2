function [M1,M2,M3,M4,Mmax] = cubeMoments(BASENAME)

% [M1,M2,M3,M4,Minf] = cubeMoments(BASENAME)

FILES = dir([BASENAME '*.fits' ])
N = length(FILES);

% Find the mean in a first pass...

M1 = single(0);
M2 = single(0);
M3 = single(0);
M4 = single(0);
Mmax = single(0); 

for n=1:N
	modprint(n,10);
	CUBE = fits_read_image(FILES(n).name);

	M1 = M1 + CUBE;
	Mmax = max(CUBE,Mmax);
end
M1 = M1/N;

for n=1:N
	modprint(n,10);
	CUBE = fits_read_image(FILES(n).name);
    CUBE = CUBE - M1;
    
	M2 = M2 + CUBE.^2;
	M3 = M3 + CUBE.^3;
	M4 = M4 + CUBE.^4;
end

M2 = M2/N;
M3 = M3/N;
M4 = M4/N;

save([BASENAME '_WS_Moments.mat'],'M*','BASENAME','FILES');
