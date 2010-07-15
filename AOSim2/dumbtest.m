Ntest = 50

thxy = randn(Ntest,2)/206265;
[X,Y] = COORDS(A);

for n=1:Ntest
    fprintf('%d ',n);
    if(mod(n,20)==0)
        fprintf('\n');
    end
    MYSTAR = [thxy(n,:) 1]*1e10;
    opl = ATMO.OPL_(X,Y,0,MYSTAR);
    
end
