function CN2 = HufnagelValley(ALTITUDES,WIND,Cn2_ground)

% CN2 = HufnagelValley(ALTITUDES,WIND,[Cn2_ground])
% ALTITUDES are in m.


GROUND_LAYER    = 0.1; % km
LOW_LAYER       = 1.5; % km

if(nargin<3)
    Cn2_ground = 1.7e-14;
end

Cn2_low         = 2.7e-16;
Cn2_Tropo       = 8.148e-26;

if(nargin<2)
    WIND = 50; % m/s Default rms wind in the jet stream.
end

Z = ALTITUDES/1000;

CN2 = Cn2_ground*exp(-Z/GROUND_LAYER) ...
    + Cn2_low*exp(-Z/LOW_LAYER) ...
    + Cn2_Tropo * WIND.^2 * Z.^10 .* exp(-Z);



