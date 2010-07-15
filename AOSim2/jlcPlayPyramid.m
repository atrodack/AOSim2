[HALO,thx,thy] = Fwfs.mkHALO(1,0.005);
[THX,THY]=meshgrid(thx,thy);
HALO(isnan(HALO)) = 0;

% CH = 200 + (-63:64);
CH = 200 + (-47:48);
RNG=-15:3:15;

% Iwfs00 = abs(fftshift2d(fft2(fftshift2d(HALO.*(THX>thx(nx)).*(THY>thy(ny)))))).^2;
Iwfs00 = abs(fftshift2d(fft2(fftshift2d(HALO.*(THX>=0).*(THY>=0))))).^2;
Iwfs01 = abs(fftshift2d(fft2(fftshift2d(HALO.*(THX>=0).*(THY<0))))).^2;
Iwfs10 = abs(fftshift2d(fft2(fftshift2d(HALO.*(THX<0).*(THY>=0))))).^2;
Iwfs11 = abs(fftshift2d(fft2(fftshift2d(HALO.*(THX<0).*(THY<0))))).^2;

Iwfs00 = Iwfs00(CH,CH);
Iwfs01 = Iwfs01(CH,CH);
Iwfs10 = Iwfs10(CH,CH);
Iwfs11 = Iwfs11(CH,CH);

Ipyr = [Iwfs00,Iwfs01;Iwfs10,Iwfs11];

Ipyramid(:,:,n) = Ipyr;

% figure(2);
% clf;
% colormap(gray);
% imagesc(Ipyr.^(1/2));sqar;colorbar;axis xy;
%  title(sprintf('Pyramid (No wobble) (\\lambda=%.2g\\mu{}m, \\gamma=2) t=%.3f',...
%         Fwfs.lambda*1e6,t));
% drawnow;
% % Iwfs = 0;
% % for nx=200+RNG
% %     for ny=200+RNG
% %         Iwfs = Iwfs + abs(fftshift2d(fft2(fftshift2d(HALO.*(THX<thx(nx)).*(THY>thy(ny)))))).^2;
% %     end
% % end
% saveJPEG(sprintf('/tmp/examplePyr_%04d.jpg',n));

