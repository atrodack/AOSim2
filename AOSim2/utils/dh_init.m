% dh_init: Initialize disk harmonic global variables
%
% Norman Mark Milton                             May 6, 2008
%
global dh_bjprime_zero_tab
% load ([mydatadir() '/dh/' 'besseljprimezeros200.dat']);
load ('data/besseljprimezeros200.mat');
dh_bjprime_zero_tab  = besseljprimezeros200;
clear besseljprimezeros200;
