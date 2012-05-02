function CCD = downsampleCCD(ccd,nbin,mbin)

% function CCD = downsampleCCD(ccd,nbin,mbin)

N = floor(size(ccd,1)/nbin);
M = floor(size(ccd,2)/mbin);

CCD = zeros(N,M);

for n=1:N
	for m=1:M
	 CHUNK = ccd(nbin*(n-1)+(1:nbin),mbin*(m-1)+(1:mbin));
	 CCD(n,m) = sum(sum(CHUNK));
	end
end

return

