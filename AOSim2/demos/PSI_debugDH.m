amp = Fwfs.lambda/4;

NTH = 4;
NZ = 500;

% ABERCUBE = [];
% SLOPES = [];
% DMCUBE = [];
% COEFS = [];

ABERCUBE = zeros([Fwfs.size NZ]);
SLOPES = zeros(2*RECON.WFS.nSubAps,NZ);
DMCUBE = zeros([Fwfs.size NZ]);
COEFS = zeros([2 NZ]);

ncase = 1;
for n=1:7
	for m=0:n
		for th=(1:NTH)*2*pi/NTH
			ABER.zero;
			ABER.addDiskHarmonic(n, m,amp*cos(th));
			ABER.addDiskHarmonic(n,-m,amp*sin(th));
			ABERCUBE(:,:,ncase) = ABER.grid;
			
			WFS.sense(Fwfs.planewave*A*ABER);
			SLOPES(:,ncase) = WFS.slopes;
			
			% 	Fwfs.show;
			WFS.quiver;
			bigtitle(sprintf('Z(%d,%d)',n,m));
			setFoV(4);
			drawnow;
			
			DM.setActs(RECON.RECONSTRUCTOR * WFS.slopes);
			DMCUBE(:,:,ncase) = DM.grid;
			COEFS(:,ncase) = [n,m];
			
			ncase = ncase + 1;
		end
	end
end

ABERCUBE(:,:,ncase:end) = [];
DMCUBE(:,:,ncase:end) = [];
COEFS(:,ncase:end) = [];
SLOPES(:,ncase:end) = [];
