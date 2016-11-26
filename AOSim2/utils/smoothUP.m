function OUTPUT = smoothUP(X,dx)

% OUTPUT = smoothUP(X,dx)
% Slow turn on at 0;

OUTPUT = (erf(X*(pi/dx))+1)/2;

