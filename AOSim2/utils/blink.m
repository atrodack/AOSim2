function blink(IMG1,IMG2,tau,N,RANGE)

% blink(IMG1,IMG2,tau,N,RANGE)

if(nargin<5)
    MIN = min(min(IMG1(:)),min(IMG2(:)));
    MAX = max(max(IMG1(:)),max(IMG2(:)));
    RANGE = [MIN MAX];
end

if(nargin<4)
    N = 100
end

if(nargin<3 || isempty(tau))
    tau = 1/8
end

for n=1:N
    imagesc(IMG1,RANGE);sqar;colorbar;drawnow;
    pause(tau);
    imagesc(IMG2,RANGE);sqar;colorbar;drawnow;
    pause(tau);
end
