amp = Fwfs.lambda/4;

NZ = 408;

% ABERCUBE = [];
% SLOPES = [];
% DMCUBE = [];
% COEFS = [];

ABERCUBE = zeros([Fwfs.size NZ]);
SLOPES = zeros(220,NZ);
DMCUBE = zeros([Fwfs.size NZ]);
COEFS = zeros([2 NZ]);

ncase = 1;
for n=1:8
	for m=mod(n,2):2:n
		for th=(0:16)*2*pi/16
			ABER.zero;
			ABER.addZernike(n, m,amp*cos(th));
			ABER.addZernike(n,-m,amp*sin(th));
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
