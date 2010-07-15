% number of actuators, radius in degrees from center (seen from sphere
% center)

D = 8.4;
LambdaRef = 0.85e-6;

% Nr = 14;
meta_r = (1:Nr)';
N1 = 6; %this needs to be an integer ~2pi to make it have isotropic density.

meta = [round(N1*meta_r) meta_r]; 

Nact = sum(meta(:,1)) + 1;
lact_ = sqrt(pi*(D/2)^2/Nact);
R0=scaleR0(0.15,500e-9,LambdaRef);

fprintf('This %d ring DM has %d actuators with a mean spacing of %.3fm. ',Nr, Nact,lact_);
fprintf('This is %.3f r0 at %.2f micron.\n', lact_/R0,LambdaRef*1e6);

Actuators = [];

for arow=1:size(meta,1)
	r = meta(arow,2);
	N = meta(arow,1);
	for n=1:N
		th = (n-1)/N*2*pi;
		Actuators(end+1,:) = r*[cos(th) sin(th)];
	end
end

Actuators(end+1,:) = [0 0];

Act = Actuators / Nr * D / 2; 
plot(Act(:,1),Act(:,2),'o'); sqar;
drawCircles(D/2,[0 0],1,'k-');
setFoV(D/2);
saveJPEG(sprintf('GMT_DM%d_%dRings_jlc_proposal.jpg',Nr,Nact),120);

voronoi(Act(:,1),Act(:,2)); sqar;
drawCircles(D/2,[0 0],1,'k-');
setFoV(D/2);
saveJPEG(sprintf('GMT_DM%d_%dRings_voronoi_jlc_proposal.jpg',Nr,Nact),120);
