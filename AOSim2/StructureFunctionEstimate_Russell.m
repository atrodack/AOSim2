load data/NewGMTPupil.mat

r0_1= .1858;
r0_2 = .3089;
L0_large = 10^6;
L0_small = 30;

WFlow = AOScreen(2048,r0_1,500e-9);     
WFhigh = AOScreen(2040,r0_2,500e-9);  

WFlow.setOuterScale(L0_large);
WFhigh.setOuterScale(L0_small);

dphase_low = 0;
dphase_high = 0;
s_real = .02:.04:25;

for i = 1:100
    
    WFlow.make;
    WFhigh.make;

    [dPhase_meanSquare,s,dPhase_meanSquareSigma,dPhase_] = WFlow.SFestimate(A,1000);
    dphase_low = dphase_low + interp1(s,dPhase_meanSquare,s_real);

    [dPhase_meanSquareH,sH,dPhase_meanSquareSigmaH,dPhase_H] = WFhigh.SFestimate(A,1000);
    dphase_high = dphase_high + interp1(sH,dPhase_meanSquareH,s_real);

end
dphase_low = dphase_low ./i;
dphase_high = dphase_high./i;
%%
figure;
loglog(s_real,dphase_low);
hold on;
plot(s_real([1 end]),[1 1]*6.88,'k--');
plot([1 1]*r0_1,[0.1 100],'k-');
plot(s_real,6.88*(s_real/r0_1).^(5/3),'g--');
%%%added by RPK
legend('MeanSquare','6.88','r0','6.88 * (s/r0)^5^/^3','location','best')
title(['r0: ', num2str(r0_1),' outerscale: ',num2str(L0_large)])

figure;
loglog(s_real,dphase_high);
hold on;
plot(s_real([1 end]),[1 1]*6.88,'k--');
plot([1 1]*r0_2,[0.1 100],'k-');
plot(s_real,6.88*(s_real/r0_2).^(5/3),'g--');
%%%added by RPK
legend('MeanSquare','6.88','r0','6.88 * (s/r0)^5^/^3','location','best')
title(['r0: ', num2str(r0_2),' outerscale: ',num2str(L0_small)])
