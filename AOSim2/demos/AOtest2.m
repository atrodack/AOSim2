STAR = [0 0 1e6];
ATMO.BEACON = STAR;
ATMO.layers{1}.Wind = [5 0];
ATMO.layers{2}.Wind = [1 -1]*30;

Fwfs = AOField(A);
F.FFTSize = 2048*[1 1]
Fwfs.lambda = RECON.lambda;

F = AOField(A);
F.lambda = AOField.KBAND;

DM.setActs(0);
% DM.addRippleActs(.3*[1 1],500e-9,0);
% DM.addRippleActs(.5*[-1.2 .4],300e-9,pi/2);
% DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);
% DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);
% DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);
% DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);
% DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);

touch(DM);

colormap(gray);

A.trueUp;
% A.segList{3}.piston = 300e-9;
% A.segList{5}.tiptilt = 0.5/206265 * [1 1];
A.touch;

mask = (A.grid>0.5);

Ipeak = 0;

f = F.grid_(:);
N1=2;N2=2;
for n=1:1000
	ATMO.time = n/1000;
	% WFS light...
	WFS.sense(Fwfs.planewave*ATMO*A*DM);
	DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
	DM.removeMean;
	
	% Science light...
	f_last = f;
	F.planewave*ATMO*A*DM;
	f = F.grid_(:);
	
	subplot(N1,N2,1);
	%F.plotPSF(0.5,[-3 0],0.01);
	%imagesc(sqrt(F.mkPSF(0.5,0.01)));sqar;
	FOV = 0.2;
	RNG = FOV * [-1 1];
	PSF = F.mkPSF(FOV,FOV/256);
	Ipeak = max(Ipeak,max(PSF(:)));
	imagesc(RNG,RNG,sqrt(PSF/Ipeak));sqar;
	%axis off;
    %title(sprintf('sqrt(PSF) (\\gamma=2) t=%.3f',ATMO.time));
    title(sprintf('H-band PSF (\\gamma=2) t=%.3f',ATMO.time));
    xlabel('arcsecs');
    ylabel('arcsecs');
    %ylabel('new phase');
    
	subplot(N1,N2,2);
	A.show;
	WFS.quiver(1);
% 	xlim([-1 1]*D/2*1.2);
%     ylim([-1 1]*D/2*1.2);
    
    % 	subplot(N1,N2,3);
    % 	% 	imagesc(F.interferometer(1));sqar;
    % 	plot(angle(f_last),angle(f),'.','MarkerSize',4);
    % 	axis square;
    % 	xlim([-1 1]*pi);
    % 	ylim([-1 1]*pi);
    % 	xlabel('old phase');
    % 	ylabel('new phase');
    % 	grid;
    
    subplot(N1,N2,4);
    % 	plot(F.grid_(:),'.','MarkerSize',1);
    % 	xlim([-1 1]*1);
    % 	ylim([-1 1]*1);
    % 	axis square;
    % 	title('pupil field complex amplitude');
    % 	hist(angle(F.grid_(:)),linspace(-pi,pi,64));
    xScale = linspace(-pi,pi,64);
    % bar(xScale,histc(angle(F.grid_(:)),xScale,'k').
   
    g=F.grid_;
	binDat = histc(angle(g(mask)),xScale);
	[vals,indx] = max(binDat);
	phase0 = xScale(indx);
	binDat = histc(angle(exp(-1i*phase0)*g(mask)),xScale);
	%bar(xScale,binDat);
	plot(xScale,binDat,'k.');
    %ylim([0 2500]);
    ylim([0 8000]);
	title('Science wavelength phase histogram');
	xlabel('pupil phase');
	ylabel('frequency');
	
	
    subplot(N1,N2,3);
    % 	plot(F.grid_(:),'.','MarkerSize',4);
	% 	xlim([-1 1]*1);
	% 	ylim([-1 1]*1);
    % 	axis square;
    % 	F.show;
    % 	title('pupil field complex amplitude');
    imagesc(sqrt(F.interferometer(1)),[0 2.5]);sqar;axis off;
    title('Science Band Interferometer');
    
    
% 	drawnow;
	saveJPEG(sprintf('/tmp/FRAME_%04d.jpg',n),160);
	
end
