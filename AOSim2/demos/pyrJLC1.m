% OFFSET = 10;

CH=1024+(-64:64);
% Q = circshift(ifft2(ifftshift2d(PSI.*(KX<10*dK))),[1 1]*128);Q=Q(1:256,1:256);subplot(1,2,1);plotCAmpl(PSI(CH,CH));subplot(1,2,2);imagesc(abs(Q));sqar;axis xy;
% Q = circshift(ifft2(ifftshift2d(PSI.*(KX<10*dK))),[1 1]*128);Q=Q(1:256,1:256);subplot(1,2,1);plotCAmpl(F.subGrid(CH,CH));subplot(1,2,2);imagesc(abs(Q));sqar;axis xy;
% Q = circshift(ifft2(ifftshift2d(PSI.*(KX<10*dK))),[1 1]*128);Q=Q(1:256,1:256);subplot(1,2,1);plotCAmpl(F.grid);subplot(1,2,2);imagesc(abs(Q));sqar;axis xy;
% Q = circshift(ifft2(ifftshift2d(PSI.*(KX<10*dK))),[1 1]*128);Q=Q(1:256,1:256);subplot(1,2,1);[gx,gy]=gradient(F.grid);plotCAmpl(gx);subplot(1,2,2);imagesc(abs(Q));sqar;axis xy;
% Q = circshift(ifft2(ifftshift2d(PSI.*(KX<10*dK))),[1 1]*128);Q=Q(1:256,1:256);subplot(1,2,1);[gx,gy]=gradient(F.grid);plotCAmpl(gx,1/4);subplot(1,2,2);imagesc(abs(Q));sqar;axis xy;
% Q = circshift(ifft2(ifftshift2d(PSI.*(KX<10*dK))),[1 1]*128);Q=Q(1:256,1:256);subplot(1,2,1);[gx,gy]=gradient(F.grid);imagesc(abs(gx));colorbar;sqar;subplot(1,2,2);imagesc(abs(Q));sqar;axis xy;

% PSI = F.fft;
 
Q = circshift(ifft2(ifftshift2d(F.fft.*(KX<OFFSET*dK))),[1 1]*128);
Q1 = Q(1:256,1:256);
Q = circshift(ifft2(ifftshift2d(F.fft.*(KX>OFFSET*dK))),[1 1]*128);
Q2 = Q(1:256,1:256);

I1 = abs(Q1).^2;
I2 = abs(Q2).^2;

subplot(2,2,1);

[gx,gy] = gradient(F.grid);
plotCAmpl(gx,1/4)
% imagesc(-abs(gx),[-.2 0]);
% colorbar;
% sqar;
axis xy

% F.show;

subplot(2,2,2);

imagesc(I1);
sqar;
axis xy;

subplot(2,2,3);

imagesc(I1-I2);
sqar;
axis xy;

subplot(2,2,4);

imagesc(I2);
sqar;
axis xy;
