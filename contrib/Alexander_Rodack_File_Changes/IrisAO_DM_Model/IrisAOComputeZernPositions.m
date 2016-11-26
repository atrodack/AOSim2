function [ PTTpos ] = IrisAOComputeZernPositions( wavelength, Zernike_Modes, Zernike_Coefficient_waves )
%[ PTTpos ] = IrisAOComputeZernPositions( wavelength, Zernike_Numbers, Zernike_Coefficient_waves )
%
% Help Section
% Inputs:
%wavelength =  the wavelength. Must be specified in meters!
%            
%Zernike_Modes = a vector of any number of integers between 1 and 23,
%                  corresponding to the single number indexing of Zernikes
%                  (not Noll, just follow the standard Pyramid)
%
%Zernike_Coefficient_waves = Zernike Coefficients in waves.  The function
%                            will error out if the number of inputs in the
%                            Zernike_Modes vector is different than the
%                            number of inputs in this vector. This means 
%                            that each Mode should get its own coefficient. 
%                            The order of coefficients should be the same 
%                            as the order for the modes.
%
% Output:
%PTTpos is a 37x3 matrix that contains the computed piston, tip, and tilt
%values for each segemnt.  Column 1 is the piston values (in units of
%meters), Column 2 is the tip values (in units of radians), and
%Column 3 is the tilt values (again in units of radians). 
%
%

lambda = wavelength*10^6;
ZernModes = Zernike_Modes;
ZernCoeffs = Zernike_Coefficient_waves;

if length(ZernModes)~=length(ZernCoeffs)
    error('Must be a Zernike Mode for each Zernike Coefficient');
end

NumSegments = 37;

% Initalize PTTpos
PTTpos = zeros(NumSegments,3);

% Initialize Coefficient Array (needs to be 21x1)
ZernikeCoefficientArray = zeros(21,1);

% Convert Coefficients from waves to distances
ZernCoeffs = ZernCoeffs*lambda;

% Add in Coefficients
ZernikeCoefficientArray(ZernModes) = ZernCoeffs;

% Load in Zernike Basis Data and Get PTT Positions
PTTPositionArray = GetPTT489ZernikePostions(ZernikeCoefficientArray);

% Shrink to the size of PTT111
PTTPositionArray = PTTPositionArray(1:37,:);
PTTpos = PTTPositionArray;
PTTpos(:,1) = PTTpos(:,1)/2; %make the vector fit how AOSim2 needs to see it
PTTpos(:,1) = PTTpos(:,1) .* 10^-6; %make the values in units of meters
PTTpos(:,2) = PTTpos(:,2) .* 10^-3; %make the values in units of radians
PTTpos(:,3) = PTTpos(:,3) .* 10^-3; %make the values in units of radians


end

