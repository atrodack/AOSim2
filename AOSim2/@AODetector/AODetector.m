classdef AODetector < handle
    properties(GetAccess='public', SetAccess='public')
        numPixels=NaN;
        readNoise;
        darkCurrent;
        fullWell;
        gain;
%         spacing;
        quantumEfficiency;
        pixelSize;  %PMH added to get old LBT script to work 100512
    end
    methods
        % constructor
        function DET = AODetector() 

            %DET = DET@AOGrid(npix);
            DET.numPixels=NaN;
            DET.readNoise=0;
            DET.darkCurrent=0;
            DET.fullWell = 0;
            DET.gain = 0;
            DET.pixelSize = 0;
            DET.quantumEfficiency = 1;
            
        
        end
        
        function DET = SetNumPixels(DET,Nx,Ny)
            
            if nargin < 3
                Ny = Nx;
            end
            DET.numPixels = [Ny,Nx];
        end
        function psf = CreatePSF(DET,f,useNoise)
            Ny = DET.numPixels(1);
            Nx = DET.numPixels(2);
            if nargin < 3
                useNoise = false;
            end
            psf = abs(fftshift2d(fft2(f,Ny,Nx))).^2 ;
            if useNoise == 1 
                % Poisson Shot Noise
                Sum0 = sum(psf(:));
                N0 = sum(sum(abs(f)));
                psf = psf*(1e-12*N0/Sum0);
                psf = 1e12*imnoise(psf,'poisson');
                % Read Noise
                psf = psf + randn(size(psf))*DET.readNoise;
                % Dark Current? psf = psf + randn(size(psf)) *
                % DET.darkCurrent*integrationTime; Where do i get the
                % integration time?
            end
        end
    end
    
end