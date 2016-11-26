function IMG1 = viewThrough(ATMO,IMG,pixel_size,RANGE)
% IMG1 = viewThrough(ATMO,IMG,pixel_size,RANGE)
% 
% JLCodona's version of the AOAtmo image warper.  
% This is a prototype for an AOAtmo visualization method.
% 
% 20151218 JLCodona.

ATMO.BEACON = [0 0 RANGE];
DOWNRANGE = max(0,RANGE-ATMO.z);

r0 = ATMO.totalFriedScale;

% STRIDE = max(1,round(r0/pixel_size));
STRIDE = max(1,round(ATMO.lambdaRef/r0*max(0,(RANGE-ATMO.z))/pixel_size/2))
STRIDE = 20

IMG = squeeze(IMG(:,:,1)); % Only monochrome at the moment.
CENTER = round(size(IMG)/2);
[X1,X2,~] = mkImageCoords(IMG,pixel_size,CENTER);

X1_ = X1(1:STRIDE:end,1:STRIDE:end);
X2_ = X2(1:STRIDE:end,1:STRIDE:end);

N1 = size(X1_,1);
N2 = size(X1_,2);

BINNING = round(r0./ATMO.dx)
% BINNING = 10

Q = downsampleCCD(ATMO.grid,BINNING,BINNING);

CTILTS = zeros([size(Q),N1,N2]);
clear Q

for n1=1:N1
    fprintf('%d/%d: ',n1,N1);
    for n2=1:N2
        fprintf('%d ',n2);
        ATMO.setBeacon([X2_(n1,n2), X1_(n1,n2), RANGE]);
        [g2,g1] = gradient(ATMO.grid,ATMO.dx);
        CTILTS(:,:,n1,n2) = downsampleCCD(g1+1i*g2,BINNING,BINNING)/BINNING^2;
        %imagesc(WAVEFRONTS(:,:,n1,n2));sqar;axis xy;axis off;drawnow;
    end
    fprintf('\n');
end

% plotComplex(squeeze(mean(mean(CTILTS))),2);

IMG1 = zeros(size(IMG));

for na1=1:size(CTILTS,1)
    fprintf('\n%d: ',na1);
    for na2=1:size(CTILTS,2)
        fprintf('%d ',na2);
        CMORPH = interp2(X2_,X1_,squeeze(CTILTS(na1,na2,:,:)),X2,X1);
        MORPH1 = min(max(1,round((X1 + real(CMORPH*DOWNRANGE))/pixel_size)+CENTER(1)),size(IMG,1));
        MORPH2 = min(max(1,round((X2 + imag(CMORPH*DOWNRANGE))/pixel_size)+CENTER(2)),size(IMG,2));
        
        %dIMG1 = zeros(size(IMG));
        for n=1:numel(IMG)
            IMG1(n) = IMG1(n) + IMG(MORPH1(n),MORPH2(n));
            %imagesc(IMG1);sqar;drawnow;
        end
        %imagesc(IMG1);sqar;drawnow;
    end
end

IMG1 = IMG1 / (size(CTILTS,1)*size(CTILTS,2));

fprintf('\n');
