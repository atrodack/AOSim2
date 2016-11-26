% A demonstration of using AOSim2 to watch the evolution of various
% statistics beyond a Kolmogorov phase screen.
% Use the Kuiper 61" to view the field.
% 
% 20150221 JLCodona

SAVE_IMGS = false;

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.

D = 1.54;
secondary = 14.5/100;

SPACING = 0.01;            % fine spacing makes nice pupil images but is really overkill.
aa = SPACING;              % for antialiasing.
% aa = 0.04;
spider = 0.0254;
% spider = 0.01;

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
   0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

% Since this demo only uses one AOSegment, I will not use the AOAperture wrapper.  
% I will only use AOSegment.  This is fine for simple pupils.

A = AOSegment;
A.spacing(SPACING);
A.name = 'Kuiper 61inch Primary';
A.pupils = PUPIL_DEFN;
A.make;

clf;
colormap(gray);

% A.show;
% input 'Press ENTER to Continue...'
% pause(3);

%% Make a Kolmogorov phase screen.
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
%TURBULENCE.lambdaRef = AOField.VBAND; %This is the default.

TURBULENCE.spacing(.02);
TURBULENCE.setR0(r0); 
TURBULENCE.make;

% TURBULENCE.show;
% input 'Press ENTER to Continue...'
% pause(3);

%% Make an AOField object.

F = AOField(A);
F.resize(1024); % make it big to study the field before the pupil.
F.FFTSize = 1024; % Used to compute PSFs, etc.
F.lambda = lambda;

[x,y] = F.coords;

F.planewave*A;
% F.show;

% input 'Continue...'

% This adds a reference wave to the field and computes the intensity.
% imagesc(x,y,F.interferometer(1),[0 3]);
% sqar;
% axis xy;
% drawnow;
% input 'Continue...'

THld = F.lambda/D * 206265; % Lambda/D in arcsecs.

% 
F.planewave*A; % Just go through the pupil.
[PSF,thx,thy] = F.mkPSF(5,THld/4);
PSFmax = max(PSF(:)); % Save for normalizing.

PSF = PSF/PSFmax; % make the brightest value =1.

imagesc(thx,thy,log10(PSF),[-4 0]); 
axis square;
axis xy;
colorbar;

% input 'Continue...'

PIXEL_RANGE = 513 + (-400:400);
xx = x(PIXEL_RANGE);
yy = y(PIXEL_RANGE);

PUPIL_PIXELS = abs(x)<0.6*D;

%% Now propagate from the screen to a large distance.  
% No aperture because we are studying the evolution of the field.

N1=2; N2=2;

% MAX_RANGE selected so that the Fresnel scale is 1/10 of the grid.
% GRID_SIZE = max(F.extent);
% ZMAX = (GRID_SIZE/50).^2/F.lambda;

% RANGES = [1:10,20:10:100,200:100:1000,2000:1000:1e4];
% RANGES = logspace(1,log10(ZMAX),100);
RANGES = logspace(1,6,100);
SELECT = 513+(-15:15); % points for the phasor plot.

for z=RANGES
    fprintf('Range: %.1f m\n',z);
    % Assume the phase screen is in the z= plane.
    
    % Start the field over each time rather than repeatedly propagating.  
    % It minimizes numerical issues..
    
    F.planewave*TURBULENCE;
    %F.propagate(z,0); % don't filter high angles.
    F.propagate(z);
    
    PSI = F.subGrid(PIXEL_RANGE,PIXEL_RANGE);
    IRRADIANCE = PSI .* conj(PSI); 
    
    subplot(N1,N2,1);
    imagesc(xx,yy,IRRADIANCE,[0 3]);
    daspect([1 1 1]);
    axis xy;
    colorbar;
    setFoV(1.0); % set the size of the plot window.
    
    drawCircles([D secondary]/2,[0 0],1,'r-')

    if(z<1000)
        title(sprintf('Irradiance z=%.1f m',z));
    else
        title(sprintf('Irradiance z=%.2f km',z/1000));
    end
    subplot(N1,N2,2);

    %hold on;
    plot(F.subGrid(SELECT,SELECT),'k.','MarkerSize',1);
    %hold off;
    %drawCircles(1,[0 0],1,'r--');
    grid;
    xlim([-1 1]*3);
    ylim([-1 1]*3);
    daspect([1 1 1]);
    title('Small patch of Field');
    xlabel('real');
    ylabel('imag');
    
    F*A;
    subplot(N1,N2,3);
    % F.show;
    PSI = F.subGrid(PUPIL_PIXELS,PUPIL_PIXELS);
    plotComplex(PSI,1);
    axis xy;
    axis off;
    title('Complex Pupil Field');
    
    subplot(N1,N2,4);
    [PSF,thx,thy] = F.mkPSF(3,THld/5);
    imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
    daspect([1 1 1]);
    axis xy;
    colorbar;

    title('PSF');
    
    drawnow; 

    if(SAVE_IMGS)
        savePNG(sprintf('WPRM_z%.0fm.png',z),150);
    end
end

%% Do it again.... (processing saved data from the earlier loop would be 
% faster, but I can't count on how much memory you have.)

SELECT = 128:8:1024-128; % Large region for phasor plot.

% select PSD binning
[XX,YY] = meshgrid(xx,yy);
RR = sqrt(XX.^2+YY.^2);
% lbinning = log10(1); % octaves
lbinning = 1/20; 
RR_ = 10.^(round(log10(RR)/lbinning)*lbinning);
RLIST = unique(sort(RR_(:)));
KAPPA = RLIST/xx(end)*pi/F.dx;
warning('off','MATLAB:Axes:NegativeDataInLogAxis');

RANGES = logspace(1,6,50);
for z=RANGES
    fprintf('[pass 2] Range: %.1f m\n',z);

    % Assume the phase screen is in the z= plane.
    
    % Start the field over each time rather than repeatedly propagating.  
    % It minimizes numerical issues..
    
    F.planewave*TURBULENCE;
    F.propagate(z,0); % don't filter high angles.
    
    subplot(N1,N2,1);
    % IRRADIANCE = F.mag2;
    PSI = F.subGrid(PIXEL_RANGE,PIXEL_RANGE);
    IRRADIANCE = PSI .* conj(PSI); 
    
    imagesc(xx,yy,IRRADIANCE,[0 3]);
    daspect([1 1 1]);
    axis xy;
    colorbar;

    if(z<1000)
        title(sprintf('Irradiance z=%.1f m',z));
    else
        title(sprintf('Irradiance z=%.0f km',z/1000));
    end
    
    % Complex scatterplot
    subplot(N1,N2,2);
    %hold on;
    %plot(F.subGrid(SELECT,SELECT),'k.','MarkerSize',1);
    plot(PSI(1:8:end,1:8:end),'k.','MarkerSize',1);
    %hold off;
    %drawCircles(1,[0 0],1,'r--');
    grid;
    xlim([-1 1]*3);
    ylim([-1 1]*3);
    daspect([1 1 1]);
    title('Complex Field');
    xlabel('real');
    ylabel('imag');

    %% PDF
    subplot(N1,N2,3);

    %[C,B] = hist(angle(PSI(:)),360);
    [C,B] = hist(IRRADIANCE(:),linspace(0,10,1000));
    C = C/sum(C); % make histogram be indep of binning choices.
    %plot(B,C);
    semilogy(B,C);
    %semilogx(B,C);
    %loglog(B,C);
    grid;
    xlim([B(1) B(end)]);
    ylim(10.^[-6 0]);
    %title('Phase Probability');
    title('Irradiance Probability');
    xlabel('Irradiance');
    ylabel('Probability');
    
    %% Irradiance Power Spectrum
    subplot(N1,N2,4);
    
    PSD = fft2(IRRADIANCE);
    PSD = fftshift(PSD.*conj(PSD));

    PSD_ = zeros(size(RLIST));
    PSD_sigma = zeros(size(RLIST));
    for nk=1:length(RLIST)
        PSD_(nk) = mean(PSD(RR_(:)==RLIST(nk)));
        PSD_sigma(nk) = std(PSD(RR_(:)==RLIST(nk)));
    end
    
    %loglog(KAPPA,PSD_);
    %loglog(KAPPA,PSD_,'k-',KAPPA,PSD_+PSD_sigma,'y--',KAPPA,PSD_-PSD_sigma,'y--');
    loglog(KAPPA,PSD_,'k-','LineWidth',3);
    hold on;
    loglog(KAPPA,PSD_+PSD_sigma,'k--',KAPPA,PSD_-PSD_sigma,'k--');
    hold off;
    grid;
    
    xlim([0.1 1e3]);
    ylim([1 1e10]);
    
    title('Irradiance Power Spectrum');
    xlabel('spatial freq');
    ylabel('power spectrum');
        
    drawnow; 
end

warning('on','MATLAB:Axes:NegativeDataInLogAxis');
