function n=dims(Q)

% function n=dims(Q)
%  Compute the number of non-singleton dimensions in an object.
% John Codona: 20050826
% 20060206: (JLC) New smarter version. 

n = sum(size(Q)>1);
