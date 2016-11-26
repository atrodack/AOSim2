function PTTPositionArray = GetPTT489ZernikePostions(ZernikeCoefficientArray)
%==========================================================================
%	GetPTT489ZernikePostions.m
%
%	Copyright© 2012, Iris AO, Inc.
%
%	Authors: Michael A. Helmbrecht, Carl J. Kempf
%
%	Origination Date: 2/24/2012
%
%	Description:	
%       Calculates the PTT values for a PTT489 DM array using Zernike basis
%       functions up to 5th order. Zernike coefficients are in microns rms.
%       The Zernike basis functions follow the ANSI-Z80.28-2004 standard.
%
%	Function format: 
%		PTTPostionArray = GetPTT489ZernikePostions(ZernikeCoefficientArray)
%
%	Function Arguments:
%		ZernikeCoefficientArray: A [21x1] vector containing Zernike
%		coefficients in µm rms.
%
%   Funtion Outputs:
%       PTTPositionArray: A [169x3] array of segment PTT values.
%
%	
%   Usage: PTTPositionArray = GetPTT489ZernikePostions(ZernikeCoefficientArray);
%
%   Version: 0.0.0.2
%
%	NOTES:
%
%   REVISIONS:
%   o  v0.0.0.1, 2/27/12, Michael Helmbrecht
%   Loads the basis set to the calling function the first time it is
%   called. Subsequent calls get the basis set from the calling function.
%   This is ~15X faster than reading from disk everytime as v0.0.0.0 did.
%   Execution times on a 2005 vintage laptop are 65 µs once the basis set
%   has been loaded to disk.
%
%   o  v0.0.0.2, 2/27/12, Michael Helmbrecht
%   Similar concept as v0.0.0.1 but loads one large [3x169x21] array that
%   contains the basis sets for piston, tip, and tilt. This executes in ~20
%   µs.
%
%   BUGS:
% 
%==========================================================================

%--------------------------------------------------------------------------
% Intiialization
%--------------------------------------------------------------------------
NumModes = 21;

% Load the modes into the calling workspace is this is the first time
% through. Otherwise get them from the calling workspace.
try
    PTT489_PTT_Modes = evalin('caller', 'PTT489_PTT_Modes');
catch
    load('PTT489_PTT_Modes')
    assignin('caller', 'PTT489_PTT_Modes', PTT489_PTT_Modes);
end


% Check the ZernikeCoefficientArray dimension -- we want a column vector
ZSize = size(ZernikeCoefficientArray);
if ZSize(2) == 21
    ZernikeCoefficientArray = ZernikeCoefficientArray';
end
    
%--------------------------------------------------------------------------
% PTT Calculations
%--------------------------------------------------------------------------
PTTPositionArray(:,1) = squeeze(PTT489_PTT_Modes(1,:,:))*ZernikeCoefficientArray;
PTTPositionArray(:,2) = squeeze(PTT489_PTT_Modes(2,:,:))*ZernikeCoefficientArray;
PTTPositionArray(:,3) = squeeze(PTT489_PTT_Modes(3,:,:))*ZernikeCoefficientArray;

