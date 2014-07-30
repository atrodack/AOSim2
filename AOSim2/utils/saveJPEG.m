function saveJPEG(filename,rez)

%function saveJPEG(filename,[rez])
% 
% Save the current figure as a JPEG image.
% 
% rez: pixels/inch (72 is canonical)
% 
% Johanan Codona

if(nargin<2)
    rez = 72;
end

resolution = sprintf('-r%d',rez);
print(resolution,'-djpeg',filename);
