function PS = make(PS,L0,fixLF)
% PS = make(PS)
% 
% New, extensible turbulent screen maker.  
% This method computes a new random realization of a screen.
% It uses the parameters in the AOScreen class to determine the model.
% 
% Supported turbulence models (specified by static class constants):
% 
% PS = AOScreen(N);
% PS.TURBULENCE_MODEL = model_number.
% 
% The model numbers are 
% AOScreen.KOLMOGOROV    -- pure power law turbulence spectrum.
% AOScreen.TATARSKI      -- KOLMOGOROV with an inner scale.
% AOScreen.VON_KARMAN    -- TATARSKI modified to include an outer scale (L0).
% AOScreen.MODIFIED_ATMO -- VON_KARMAN modified to include the Hill inner
%                           scale enhancement.
% All models now have a variable power law PHI_n index ALPHA.  Set this in
% the AOScreen before make()'ing it.

% Start with the envelope, the random part is common to all models.

if(nargin>1)
   fprintf('WARNING: Setting the outer scale and fixLF flag in the make args is deprecated.\n');
   fprintf('WARNING: Set these params in the AOScreen object first and then call make.\n');
   
   PS.L0 = L0;
   PS.LF_FRACTAL_PATCH = fixLF;
   
end

[KX,KY] = PS.KCOORDS;
KR2 = KX.^2 + KY.^2;

kappa0 = 2*pi/PS.L0;
kappam = 5.92/PS.inner_scale;

switch PS.TURBULENCE_MODEL
    case AOScreen.KOLMOGOROV
        PSD = 0.033 * PS.Cn2 * KR2.^(-PS.ALPHA/2);
        
    case AOScreen.TATARSKI
        PSD = 0.033 * PS.Cn2 * KR2.^(-PS.ALPHA/2);
        PSD = PSD .* exp(-KR2./kappam^2); % Inner scale
        
    case AOScreen.VON_KARMAN
        PSD = 0.033 * PS.Cn2 .* (kappa0^2+KR2).^(-PS.ALPHA/2); % Outer scale
        PSD = PSD .* exp(-KR2./kappam^2); % Inner scale
                
    case AOScreen.MODIFIED_ATMO
        kappal = 3.3/PS.inner_scale;
        K = sqrt(KR2);
        PSD = 0.033 * PS.Cn2 .* (kappa0^2+KR2).^(-PS.ALPHA/2); % Outer scale
        PSD = PSD .* (1 + 1.802*(K/kappal) - 0.254*(K/kappal).^(7/6)); % Hill bump.
        PSD = PSD .* exp(-KR2./kappal^2); % Modified Inner scale
        
    otherwise % Use the Von Karman model as a default.
        PSD = 0.033 * PS.Cn2 .* (kappa0^2+KR2).^(-PS.ALPHA/2); % Outer scale
        PSD = PSD .* exp(-KR2./kappam^2); % Inner scale
        
end

PSD(PS.FAXIS_PIXEL(1),PS.FAXIS_PIXEL(2)) = 0;  % Kill DC for zero mean.

PS.zero.addNoise(1);
PS.grid(PS.grid.*sqrt(PSD));

PS.grid(PS.fft);
PS.grid(PS.grid*(24.6 / sqrt(prod(PS.spacing)) / sqrt(prod(PS.size))  ));

if(PS.LOW_FREQ_FIX) % Patch up the real part using the imag part.
    
    % TBD
    fprintf('NOTICE: The Low Freq Patch is not included in this version YET.\n');
    
end

PS.real; % Throw away the imaginary part.  

PS.touched = false;

