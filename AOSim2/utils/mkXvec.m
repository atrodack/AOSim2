function x = mkXvec(sz,dx)

%  x = mkXvec(sz,dx)

x = 1:sz;
x = fftshift(x);
x = x-x(1);

x = x*dx;
