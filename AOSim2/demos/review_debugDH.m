N1=1; N2=2;
ZSCALE = 2;

[x,y] = A.coords;
a = A.grid;
anan = a;
anan(a<0.8) = nan;

for n=1:size(COEFS,2)

    subplot(N1,N2,1);
    surf(x,y,demean(anan.*(ABERCUBE(:,:,n))*(2*pi/Fwfs.lambda)),'LineStyle','none');
    zlim([-1 1]*ZSCALE);
    daspect([1 1 ZSCALE]);
    lt=light();
    set(lt,'Position',[1 1 .5]);
    bigtitle(sprintf('Input was Disk Harmonic(%d,%d): %d',COEFS(:,n),n));

    subplot(N1,N2,2);
    surf(x,y,demean(anan.*(ABERCUBE(:,:,n)-DMCUBE(:,:,n))*(2*pi/Fwfs.lambda)),'LineStyle','none');
    zlim([-1 1]*ZSCALE);
    daspect([1 1 ZSCALE]);
    lt=light();
    set(lt,'Position',[1 1 .5]);
    bigtitle(sprintf('Input-Reconstructed'));
    drawnow;
    
%     input 'next> ';
end
