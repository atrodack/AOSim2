function CUBE = cubeConv2(CUBE,KERNEL,MODE)

% CUBE = cubeConv2(CUBE,KERNEL,[MODE])
% MODE defaults to 'same'
% see doc conv2.

if(nargin<3)
    MODE = 'same';
end

whos KERNEL
MODE
for n=1:size(CUBE,3)
    CUBE(:,:,n) = conv2(CUBE(:,:,n),KERNEL,MODE);
end
