function G = gaussian(X,sigma)
% G = gaussian(X,sigma)
% Compute a unit Gaussian at the points X with a sigma.
% sigma defaults to 1.
% 20150321 JLCodona

if(nargin<2)
    G = exp(-0.5*(X).^2);
else
    G = exp(-0.5*(X/sigma).^2);
end

