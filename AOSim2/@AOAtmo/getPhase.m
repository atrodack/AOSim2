function OPD = getPhase(ATMO,X,Y,z,XYZb,includeGeometry)

% function OPD = getPhase(ATMO,X,Y,z,BEACONxyz,includeGeometry)
% 
% Part of AOSim's AOAtmo.
% 
% JLCodona 20090221

global env;

N = length(ATMO.layers);
% fprintf('There are %d phase screens.\n', N);

OPD = zeros(size(X));
for n=1:N
% 	fprintf('AtmoLayer %d (%s): z=%g spacing=[%g %g]\n', ...
% 		n, ATMO.layers{n}.name, ATMO.layers{n}.z, ATMO.layers{n}.dXY);
	
	SZn = size(ATMO.layers{n}.phase);
	XYn = ATMO.layers{n}.XY;
	dXYn = ATMO.layers{n}.dXY;
	Zn = ATMO.layers{n}.z;
	Zb = XYZb(3);
	WIND = ATMO.layers{n}.Wind;
	
	if(Zn>Zb)
% 		fprintf('AtmoLayer %d (%s): z=%g spacing=[%g %g]\n', ...
% 			n, ATMO.layers{n}.name, ATMO.layers{n}.z, ATMO.layers{n}.dXY);
		fprintf('Screen %d is beyond the beacon.\n',n);
		continue;
	end
	
	% Construct the layer coordinates on-the-fly for interpolation.
	x = XYn(1) + (0:SZn(1)-1)*dXYn(1);
	y = XYn(2) + (0:SZn(2)-1)*dXYn(2);
	[Xn,Yn] = meshgrid(x,y);
	
	beta = (Zn-z)/(XYZb(3)-z); % fraction of the way to the beacon.
	
	Xinterp = X*(1-beta) + XYZb(1)*beta + WIND(1)*env.time;
	Yinterp = Y*(1-beta) + XYZb(2)*beta + WIND(2)*env.time;
	
% 	dOPD = interp2(Xn,Yn,ATMO.layers{n}.phase,Xinterp,Yinterp,'linear');
	dOPD = qinterp2(Xn,Yn,ATMO.layers{n}.phase,Xinterp,Yinterp);
% 	dOPD(isnan(dOPD)) = 0;
	OPD = OPD + dOPD;
end

if(includeGeometry)
	dL = sqrt((X-XYZb(1)).^2+(Y-XYZb(2)).^2+(z-XYZb(3)).^2);
	dL = dL - mean(dL(:)); % this is arbitrary, I just don't want too much wrapping going on.
	OPD = OPD + dL;
end

OPD = OPD*env.k;

return;
