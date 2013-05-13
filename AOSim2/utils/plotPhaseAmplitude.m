function plotPhaseAmplitude(I, PHASE, RANGE)

%function plotPhaseAmplitude(I,PHASE,RANGE)

% whos I PHASE RANGE

if(nargin<3)
    minI = min(I(:));
    maxI = max(I(:));
else
    minI = RANGE(1);
    maxI = RANGE(2);
end

I = (max(min(I,maxI),minI)-minI)/(maxI-minI); % should now be 0<I<1

PHASE = PHASE - 1.1;  % better alignment of colors.

RGB(:,:,1) = I .* (1+cos(PHASE-pi/2))/2;
RGB(:,:,2) = I .* (1+cos(2*PHASE-pi/3))/2;
RGB(:,:,3) = I .* (1+cos(PHASE+pi/2))/2;

imagesc(RGB,[0 1]);
sqar;



