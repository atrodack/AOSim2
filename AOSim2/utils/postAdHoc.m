% Build an explicit set of low-order patches to the reconstructor.

SCALE = 1e-6;

save RECONSTRUCTOR_BACKUP RECON

Nsegs = length(RECON.A.segList);
slopes = RECON.WFS.slopes;
TESTVECTORS = zeros(length(slopes),2*Nsegs);
clear slopes;

ntest = 0;
for n=1:Nsegs
    SEG = (RECON.DM.actuators(:,4)==n);
    ntest = ntest + 1;
    TESTVECTORS(SEG,ntest) = SCALE * demean(RECON.DM.actuators(SEG,1));
    ntest = ntest + 1;
    TESTVECTORS(SEG,ntest) = SCALE * demean(RECON.DM.actuators(SEG,2));
end

subplot(2,1,1);
imagesc(TESTVECTORS);

TESTSLOPES = RECON.processTestVectors(TESTVECTORS);

subplot(2,1,2);
imagesc(TESTSLOPES);

ITSLOPES = pseudoInv(TESTSLOPES,1e-4);
Rfix = TESTVECTORS*ITSLOPES;
RECON.RECONSTRUCTOR = RECON.RECONSTRUCTOR + 0.9*Rfix; % 0.9 is an example.
