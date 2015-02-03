clc;
close all;
clear all;


%% Assess the flag situation
useatmo = input('Set useatmo flag: ');

if useatmo == true
    DMpokescreen = false;
else
    DMpokescreen = true;
    SHWFScheck = false;
end

%% Initialize Entrance Pupil
LAMBDA = AOField.JBAND;
% load('dOTFDM.mat');
% dOTFDM = DM;
% clear DM;
make_the_MMT_AO_jlc;

thld = LAMBDA/D*206265;
FOV = 100 * thld;
PLATE_SCALE = thld/3;


%% SHWFS
WFS.qscale = 1;

% Can change SHWFScheck to false if desired
if useatmo == true
    SHWFScheck = true;
end

%% Reconstructor
if SHWFScheck == true
    % RECON.program(D,6*sqrt(2)); % Use Fourier modes. OWD is ~6 lambda/D for programming.
    RECON.zprogram(D,9);  % program using Zernikes.
    % RECON.dhprogram(D,10); % program using disk harmonics.
    
    RECON.rebuild(500);
    RECON.Nmodes;
    DM.nActs;
    fprintf('\n\n\n');
end

%% Field Initialization
F = AOField(A);
F.FFTSize = 1024;
F2 = AOField(A);
F2.FFTSize = 1024;
Fwfs = AOField(A);
Fwfs.lambda = F.lambda; 
Fwfs.FFTSize = 1024;
FSHwfs = AOField(A);
FSHwfs.lambda = RECON.lambda; 
FSHwfs.FFTSize = 1024;

%% dOTF Object
W = AOdOTF(A,FOV,PLATE_SCALE);
[x,y] = W.coords;
% test multiple finger positions/widths
numtests = 1;
xpos = [-2,-1,0,1,2];
ypos = [-2.5,-1.5,0.5,1.5,2.5];
width = [0.125,0.25,0.33,0.5,0.67];

%% Atmosphere Creation
% Use Atmosphere for aberration Flag
% useatmo = false;


% Atmosphere
if useatmo == true;
    ATMO = AOAtmo(A);
    WFlow = AOScreen(1024,0.17,500e-9);
    WFlow.name = 'Lower altitude turbulence';
    WFhigh = AOScreen(2048,0.20,500e-9);
    WFhigh.name = 'High altitude turbulence';
    ATMO.addLayer(WFlow,1000);
    ATMO.addLayer(WFhigh,8000);
    ATMO.layers{1}.Wind = [3 1];
    ATMO.layers{2}.Wind = [1 -1]*20;
    r0 = ATMO.totalFriedScale;
    th_scat = AOField.VBAND/r0*206265;
    fprintf('The total r0 is %f cm.\n',100*ATMO.totalFriedScale);
    fprintf('The seeing is %.2f arcsecs.\n',th_scat);
    % Turning this off is like using dynamic refocus.
    ATMO.GEOMETRY = false;
    
    % Guide star selection
    STAR = [0 0 1e10];
    GUIDE_STAR = STAR; % pick one.
    SCIENCE_OBJECT = STAR; % pick one.
    ATMO.BEACON = GUIDE_STAR;
end

%% DM
A.trueUp;
DM.setActs(0);
% DMpokescreen flag
% DMpokescreen = true;

%DMpokescreen flag check
if DMpokescreen == true
    clf;
    A.show;
    DM.plotActuators(1);
    colormap(gray);
    pause(2.5);
%     pokeact = [218,221];
    pokeact = randi(336,1,randi(15,1,1));
    for ii = 1:length(pokeact)
        DM.addPoke(pokeact(ii),1e-6);
        fprintf('Actuator # %d poked\n',pokeact(ii));
    end
end

touch(DM);
touch(A);

%% dOTFDM
dOTFDM = DM.copy;
touch(dOTFDM);


%% Some Grid Objects
k = 2*pi/F.lambda;
g = AOGrid(length(A.grid));
g.constant(1);












%% do the dOTF
if useatmo == true
    ps = ATMO;
else
    ps = DM;
end

W.calibrateWFS(2.5,1,0.125,Fwfs.planewave,ps);
% store figure handles
h = figure(1);

%% Close the Loop
if useatmo == true
    %% Atmo Screen
    for n = 1:1
        %             get time based on FPS value of 550
        t = n/550;
        %                         update Atmo time property
        %             ATMO.time = t-1.5;
        W.sense(Fwfs.planewave * dOTFDM * A,ATMO);
        OPL = W.Phase/k;
        g.grid(OPL);
        pistonvec = g.interpGrid(DM.actuators(:,1),DM.actuators(:,2));
        dOTFDM.bumpActs(pistonvec);
        dOTFDM.removeMean;
        dOTFDM.render;
        
        
        if SHWFScheck == true
            WFS.sense(FSHwfs.planewave*ATMO*A*DM);
            gain = 1;
            DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
            DM.removeMean;
        end
        
        % Plotting Options
% % %              Plot PSFs, OTFs, magnitude of dOTF, retrieved phase of dOTF
%         W.show;
%         
% % %                  Plot Computed vs. Incident Phase
%         clf;
%         subplot(1,2,1)
%         imagesc(x,y,W.Phase .* A.combined.grid_);
%         colormap(jet);
%         colorbar;
%         axis xy;
%         bigtitle(sprintf('Computed Phase at t = %0.4f',t),10);
%         daspect([1,1,1]);
%         subplot(1,2,2)
%         Fwfs.show;
%         imagesc(x,y,ATMO);
%         colorbar;
%         bigtitle('Incident Phase',12);
%         colormap(jet);
%         axis xy;
%         daspect([1,1,1]);
%         
% % %                         Plot Residual Phase
%         clf;
%         subplot(1,2,1)
%         imagesc(x,y,angle(Fwfs.grid_).*A.combined.grid_);
%         colormap(jet);
%         colorbar;
%         axis xy;
%         daspect([1,1,1]);
%         subplot(1,2,2)
%         imagesc(x,y,angle(Fwfs.grid_).*A.combined.grid_ - W.Phase);
%         colormap(jet);
%         colorbar;
%         axis xy;
%         daspect([1,1,1]);
%         
% % %                         Plot PSFs
%         clf;
%         gain = 1;
%         compphase = W.Phase .* A.combined.grid_;
%         actualphase = angle(Fwfs.grid_).*A.combined.grid_;
%         residphase = (actualphase - (gain*compphase));
%         actualphase = padarray(actualphase,[512,512]);
%         residphase = padarray(residphase,[512,512]);
%         img = abs(fftshift(fft2(fftshift(actualphase)))).^2;
%         img2 = abs(fftshift(fft2(fftshift(residphase)))).^2;
%         subplot(1,2,1)
%         imagesc(img);
%         axis xy;
%         daspect([1,1,1]);
%         colormap(gray);
%         colorbar;
%         bigtitle(sprintf('t = %0.4f',t),12);
%         subplot(1,2,2)
%         imagesc(img2);
%         axis xy;
%         daspect([1,1,1]);
%         colormap(gray);
%         colorbar;
%         bigtitle(sprintf('t = %0.4f',t),12);
        
        
        
        clf;
        N1 = 1; N2 = 2; N3 = 0;
        if SHWFScheck== true
            N3 = 1;
        end
        
        subplot(N1+N3,N2+N3,2);
        F2.planewave * ATMO * dOTFDM * A;
        psfplot = F2.mkPSF(FOV,PLATE_SCALE);
        [sizex,sizey] = size(psfplot);
        plotwin = floor((sizex/2))- round(floor(sizex/2)/3) : floor(sizex/2) + round(floor(sizex/2)/3);
        F2.touch;
        imagesc(psfplot(plotwin,plotwin));
        daspect([1,1,1]);
        axis xy;
        colormap(jet);
        bigtitle('Corrected PSF',20);
        
        
        subplot(N1+N3,N2+N3,1);
        k = (2*pi)/LAMBDA;
        imagesc(dOTFDM.grid .* A.grid);
        daspect([1,1,1]);
        axis xy;
        colormap(jet);
        colorbar;
        bigtitle('Calculated OPL',20);
        
        if SHWFScheck == true
            subplot(2,3,4)
            imagesc(DM.grid .* A.grid);
            daspect([1,1,1]);
            axis xy;
            colormap(jet);
            colorbar;
            bigtitle(sprintf('DM Shape from SHWFS at t = %0.3f',n),20);
            
            subplot(2,3,5)
            F.planewave*ATMO*A*DM;
            PSFplot = F.mkPSF(FOV,PLATE_SCALE);
            F.touch;
            imagesc(PSFplot(plotwin,plotwin))
            daspect([1,1,1]);
            axis xy;
            colormap(jet);
            bigtitle('AO loop PSF',20);
            
            subplot(2,3,6)
            F.planewave * ATMO * A;
            uncorrpsf = F.mkPSF(FOV,PLATE_SCALE);
            F.touch;
            imagesc(uncorrpsf(plotwin,plotwin));
            daspect([1,1,1]);
            axis xy;
            colormap(jet);
            bigtitle('Uncorrected PSF',20);
            
            
        end
        
        
        
        
        
        %% Other Commands
        drawnow;
        
%                     M(n) = getframe(h);
        
        %% Store before looping
        %             Store Stuff
        psf0{n} = W.PSF0;
        psf1{n} = W.PSF1;
        otf0{n} = W.OTF0;
        otf1{n} = W.OTF1;
        dotf{n} = W.dOTF;
        
        %           Clear dOTF properties in class
        %             W.cleardOTF;
    end
    
    
    
elseif DMpokescreen == true
    %% Use DM poke Screen
    for n = 1:1
        %             get time based on FPS value of 550
        t = n/550;
        W.sense(Fwfs.planewave * g,DM);
        OPL = W.Phase/k;
        wave = exp((1*2*pi*1i/F.lambda)*OPL);
        g.setgrid_(1*wave);
        
        %             if SHWFScheck == true
        %                 WFS.sense(FSHwfs.planewave*ATMO*A*DM);
        %                 gain = 1;
        %                 DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
        %                 DM.removeMean;
        %                 OPL2 = DM.grid;
        %                 wave2 = exp((1*2*pi*1i/F.lambda)*OPL2);
        %                 g2.setgrid_(wave2);
        %             end
        
        
        %% Plot it
        clf;
        N1 = 1; N2 = 2;
        
        subplot(N2,N2,2);
        F2.planewave * DM * A * g;
        psfplot = F2.mkPSF(FOV,PLATE_SCALE);
        [sizex,sizey] = size(psfplot);
        plotwin = floor((sizex/2))- round(floor(sizex/2)/3) : floor(sizex/2) + round(floor(sizex/2)/3);
        F2.touch;
        imagesc(psfplot(plotwin,plotwin));
        daspect([1,1,1]);
        axis xy;
        colormap(jet);
        bigtitle('Corrected PSF',20);
        
        
        subplot(N2,N2,1);
        imagesc(OPL.*A.grid);
        daspect([1,1,1]);
        axis xy;
        colormap(jet);
        colorbar;
        bigtitle('Calculated OPL',20);
        
        subplot(N2,N2,3)
        imagesc(DM.grid .* A.grid);
        daspect([1,1,1]);
        axis xy;
        colormap(jet);
        colorbar;
        bigtitle(sprintf('DM Shape at t = %0.3f',n),20);
        
        subplot(N2,N2,4)
        F.planewave * A * DM;
        PSFplot = F.mkPSF(FOV,PLATE_SCALE);
        F.touch;
        imagesc(PSFplot(plotwin,plotwin))
        % F.show;
        daspect([1,1,1]);
        axis xy;
        colormap(jet);
        bigtitle('Uncorrected PSF',20);
        
        
        %% Command it
        drawnow;
        %             M(n) = movie(h);
        
        %% Store before looping
        %             Store Stuff
        psf0{n} = W.PSF0;
        psf1{n} = W.PSF1;
        otf0{n} = W.OTF0;
        otf1{n} = W.OTF1;
        dotf{n} = W.dOTF;
        
        
        
        
    end
    
    
else
    %% No Phase Screen
    W.mkOTF(Fwfs);
    
    W.mkdOTF;
    W.mkMask;
    W.show;
    %         Store Stuff
    psf0 = W.PSF0;
    psf1 = W.PSF1;
    otf0 = W.OTF0;
    otf1 = W.OTF1;
    dotf = W.dOTF;
    
    W.cleardOTF;
end



%% Make a Movie
% Uncomment to generate movie file
% movie2avi(M, 'sorta_working')