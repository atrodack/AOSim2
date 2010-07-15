%% @file
% test file

%%
% test function returns one
% @param car input variable
% @return one

function m = doxytest(car)
  n=car;

  % normal comment
  m=n;
  

%%
% test function2 returns nothing
% @param philbert input variable
%
function subfunct(philbert)
  n=philbert/2;			% end of line comment
  m=n;

%%
% last test funciton
% @param g input param 1
% @param a input param 2
% @param d input param 3
% @return a,b value of param d

%
 function [a, b]= last_Func(g, a,d)
  a=d;
  subfunct(a);
  b=doxytest(g)
  
  