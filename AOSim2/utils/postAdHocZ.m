% Build an explicit set of low-order patches to the reconstructor.

SCALE = 1e-6;
Nmax = 2; % Highest radial order Zernike

Nsegs = length(RECON.A.segList);
slopes = RECON.WFS.slopes;
TEST_ACTS = zeros(length(slopes),2*Nsegs);
clear slopes;

nTEST_ = 0;
for n=1:Nsegs
    SEG = (RECON.DM.actuators(:,4)==n);
    nTEST_ = nTEST_ + 1;
    TEST_ACTS(SEG,nTEST_) = SCALE * demean(RECON.DM.actuators(SEG,1));
    nTEST_ = nTEST_ + 1;
    TEST_ACTS(SEG,nTEST_) = SCALE * demean(RECON.DM.actuators(SEG,2));
end

subplot(2,1,1);
imagesc(TEST_ACTS);

TEST_SLOPES = RECON.processTEST_ACTS(TEST_ACTS);

subplot(2,1,2);
imagesc(TEST_SLOPES);

