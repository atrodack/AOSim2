% Make a HV 5/7 atmosphere.
% 

dx = 0.04;
PS_SIZE = 2048;

if(~exist('WIND'))
    WIND = 30;
end

if(~exist('Cn2_Ground'))
    Cn2_Ground = 1e-14;
end

LAMBDA = AOField.VBAND;

HEIGHTS = [0,logspace(0,4.5,25)];
THICKS = diff(HEIGHTS);
THICKS(end+1) = THICKS(end);

CN2 = HufnagelValley(HEIGHTS,WIND,Cn2_Ground);

ATMO = AOAtmo(512);
ATMO.name = 'Hufnagel-Valley Atmosphere';
ATMO.spacing(dx);
ATMO.lambdaRef = LAMBDA;

for n=1:length(HEIGHTS)
    ps = AOScreen(PS_SIZE);
    ps.name = sprintf('Layer %d: %.0f m',n,HEIGHTS(n));
    ps.spacing(dx);
    ps.lambdaRef = LAMBDA;
    ps.setCn2(CN2(n),THICKS(n));
    %ps.interpolate_method = 'cubic';
    ps.interpolate_method = 'linear';
    
    ATMO.addLayer(ps,HEIGHTS(n));
    
    ATMO.layers{n}.Wind = [10 15] + [0 WIND]*HEIGHTS(n)/10000; % Jet Stream.
end

ATMO.show;

% r00 = HVATMO.totalFriedScaleStar
% NAPER = round(D/r00);
% HVATMO.spacing(r00);
% HVATMO.resize(NAPER);

ATMO.useGeometry(false);

fprintf('\nr0 for the beacon is %f m.\n',ATMO.totalFriedScale);