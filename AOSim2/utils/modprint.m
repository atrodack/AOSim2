function modprint(n,N)

% modprint(n,N)
% 
% JLCodona. I kept needing this, so...


fprintf('%5d ',n);
if(mod(n,N)==0) 
    fprintf('\n'); 
end


