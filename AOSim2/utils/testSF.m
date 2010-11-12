function [SF1,SF2] = testSF(S)

% [SF1,SF2] = testSF(S)
% This is a utility to validate the structure function of an AOScreen.
% 

SPACING = S.spacing;
SIZE = S.size;
lambda = S.lambdaRef;
k2 = (2*pi/lambda)^2;

N1 = min(100,SIZE(1)-1);
N2 = min(100,SIZE(2)-1);

SF1 = zeros(N1,2);
SF2 = zeros(N2,2);
% SF1 = zeros(N1,3);
% SF2 = zeros(N2,3);

SF1(:,1) = (1:N1)*SPACING(1);
SF2(:,1) = (1:N2)*SPACING(2);

g = S.grid_(:,1:8:end);
for n=1:N1
    dSurface2 = (g(1:SIZE(1)-n,:)-g(1+n:SIZE(1),:)).^2;
    SF1(n,2) = k2*mean(dSurface2(:));
    %SF1(n,3) = k2*std(dSurface2(:))/sqrt(numel(dSurface2));
end

g = S.grid_(1:8:end,:);
for n=1:N2
    dSurface2 = (g(:,1:SIZE(2)-n)-g(:,1+n:SIZE(2))).^2;
    SF2(n,2) = k2*mean(dSurface2(:));
    %SF2(n,3) = k2*std(dSurface2(:))/sqrt(numel(dSurface2));
end

loglog(SF1(:,1),SF1(:,2),'k-','LineWidth',2);
hold on;

% for n=1:N1
%     loglog(SF1(n,1)*[1 1],SF1(n,2)+SF1(:,3)*[-1 1],'k-');
% end
% 
loglog(SF2(:,1),SF2(:,2),'b-','LineWidth',2);
% for n=1:N2
%     loglog(SF2(n,1)*[1 1],SF2(n,2)+SF2(:,3)*[-1 1],'b-');
% end

% Plot r0 and 6.88 which is where the structure function should be at r0.
plot([1 1]*S.r0,[0.1 100],'r--',...
    SF1([1 end],1),6.88*[1 1],'r--');

plot(SF2(:,1),6.88*(SF2(:,1)/S.r0).^(5/3),'r--');

hold off;

grid;

% loglog(SF1(:,1),SF1(:,2),...
%     SF1(:,1),SF1(:,2),...
%     SF2(:,1),6.88*(SF2(:,1)/S.r0).^(5/3),'r--');
