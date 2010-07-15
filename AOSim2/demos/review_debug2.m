N1=2; N2=2;
ZSCALE = 2;

[x,y] = A.coords;
a = A.grid;
anan = a;
anan(a<0.8) = nan;

for n=1:size(ABERCUBE,3)
	clf;
	
	subplot(N1,N2,1);
	imagesc(x,y,anan.*ABERCUBE(:,:,n)/RECON.lambda*2*pi,[-1 1]);
	sqar;colorbar;
	bigtitle(sprintf('Zernike(%d,%d): %d',COEFS(:,n),n));
	
	subplot(N1,N2,2);
	imagesc(x,y,anan.*DMCUBE(:,:,n)/RECON.lambda*2*pi,[-1 1]);
	sqar;colorbar;
	
	%subplot(N1,N2,[3 4]);
	subplot(N1,N2,3);
	imagesc(x,y,anan.*(ABERCUBE(:,:,n)-DMCUBE(:,:,n))/RECON.lambda*2*pi,0.2*[-1 1]);
	sqar;colorbar;
	
	subplot(N1,N2,4);
	surf(x,y,anan.*(ABERCUBE(:,:,n)-DMCUBE(:,:,n))/RECON.lambda*2*pi,...
		'LineStyle','none');
	lt = light();
	set(lt,'Position',[1 1 0.5]);
	zlim([-1 1]*ZSCALE);
    daspect([1 1 ZSCALE]);
	
	
	drawnow;

% 	input 'next? '
end
