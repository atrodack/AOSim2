% cd AOSim2/
make_the_MMT_AO_jlc
run_myMMT_justStartup
% addpath /data/mfitsio/
pwd
% JLC_OTFpipeline
% ls
% imagesc(PSF)
[PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
FoV = lambda/D*206265*12
lambda = F.lambda
FoV = lambda/D*206265*12
[PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
dth=FoV/256
[PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
imagesc(PSF)
imagesc(PSF);sqar;
imagesc(log10(normalize(PSF)),[-4 0]);sqar;
F.planewave*ATMO*A;
imagesc(log10(normalize(PSF)),[-4 0]);sqar;
[PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
imagesc(log10(normalize(PSF)),[-4 0]);sqar;
[X,Y,r] = mkImageCoords(PSF,1,size(PSF)/2);
DOPPLER = X;
LAG = r.^2;
plot(DOPPLER(:),LAG(:),'.','MarkerSize',1);

LDIMAGE(2048,600);
LDIMAGE=zeros(2048,600);
for n=1:numel(PSF)
    LDIMAGE(round(LAG(n)/1000)+1,round(DOPPLER(n)+300)) = PSF(n);
end
imagesc(LDIMAGE)
colorbar
imagesc(conv2(LDIMAGE,ones(5)))
imagesc(log10(conv2(LDIMAGE,ones(5))))
imagesc(log10(conv2(LDIMAGE,ones(5))));axis xy
imagesc(log10(conv2(LDIMAGE,ones(5))));axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(5))),[-3 0]);axis xy;colorbar;
LDIMAGE=zeros(2048,600);
for n=1:numel(PSF)
    LDIMAGE(round(LAG(n)/200)+1,round(DOPPLER(n)+300)) = PSF(n);
end
imagesc(log10(conv2(normalize(LDIMAGE),ones(5))),[-3 0]);axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(5))),[-2 0]);axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(1))),[-2 0]);axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(1))),[-5 0]);axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-5 0]);axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-4 0]);axis xy;colorbar;
LDIMAGE=zeros(2048,600);
for n=1:numel(PSF)
    LDIMAGE(round(LAG(n)/50)+1,round(DOPPLER(n)+300)) = PSF(n);
end
imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-4 0]);axis xy;colorbar;
imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-3 0]);axis xy;colorbar;
ATMO.time
make_the_MMT_AO_jlc
make_the_GMT_AO
run_myGMT_justStartup
F.planewave*ATMO*A;
[PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
FoV = lambda/D*206265*12
D
dth=FoV/1024
[PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
DOPPLER = X;
[X,Y,r] = mkImageCoords(PSF,1,size(PSF)/2);
DOPPLER = X;
LAG = r.^2;
LDIMAGE(round(LAG(n)/50)+1,round(DOPPLER(n)+300)) = PSF(n);
plot(DOPPLER(:),LAG(:),'.','MarkerSize',1);
2.5e6/1024
2.5e6/3000
for n=1:numel(PSF)
    LDIMAGE=zeros(2048,2*1110);
    LDIMAGE=zeros(1024,2*1110);
    for n=1:numel(PSF)
        LDIMAGE(round(LAG(n)/1000)+1,round(DOPPLER(n)+1050)) = PSF(n);
    end
    imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-3 0]);axis xy;colorbar;
    imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-2 0]);axis xy;colorbar;
    pwd
    dt = .01
    for n=1:100
        ATMO.time = -0.5 + n*dt;
        F.planewave*ATMO*A;
        [PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
        LDIMAGE=zeros(1024,2*1110);
        for nn=1:numel(PSF)
            LDIMAGE(round(LAG(nn)/1000)+1,round(DOPPLER(nn)+1050)) = PSF(nn);
        end
        imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-2 0]);axis xy;colorbar;
        biglabels('Doppler (biased)','lag');
        bigtitle(sprintf('sim time = %.3f',ATMO.time);
        for n=1:100
            ATMO.time = -0.5 + n*dt;
            F.planewave*ATMO*A;
            [PSF,thx,thy] = F.mkPSF(FoV,dth);PSF(isnan(PSF)) = 0;
            LDIMAGE=zeros(1024,2*1110);
            for nn=1:numel(PSF)
                LDIMAGE(round(LAG(nn)/1000)+1,round(DOPPLER(nn)+1050)) = PSF(nn);
            end
            imagesc(log10(conv2(normalize(LDIMAGE),ones(3))),[-2 0]);axis xy;colorbar;
            biglabels('Doppler (biased)','lag');
            bigtitle(sprintf('sim time = %.3f',ATMO.time));
            saveJPEG(sprintf('/tmp/FRAME_%04d.jpg',n));
        end
