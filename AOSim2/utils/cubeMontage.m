function IMG = cubeMontage(CUBE,N1,N2)

% IMG = cubeMontage(CUBE,N1,N2)
% 
% JLCodona 20120427

% CUBE(:,:,N1*N2+1:end) = [];

N3 = size(CUBE,3);

if(nargin<3)
    N2 = ceil(N3/N1);
end

IMG = zeros(N1*size(CUBE,1),N2*size(CUBE,2));


CH1 = 1:size(CUBE,1);
CH2 = 1:size(CUBE,2);

for n=1:min(N1*N2,N3)
    start1 = 1+mod(n-1,N1)*size(CUBE,1);
    start2 = round((n-mod(n-1,N1)-1)/N1)*size(CUBE,2);

    IMG(start1+CH1,start2+CH2) = CUBE(:,:,n);
end
