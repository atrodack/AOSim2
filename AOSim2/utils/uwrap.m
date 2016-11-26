function uwphase = uwrap(phase,ALGO)

% function uwphase = uwrap(phase,[algorithm])
% Unwrap phase using fancy external programs.
% Available algorithms are 
% 1 gold,
% 2 unwt,
% 3 unmg,
% 4 flyn,
% 5 fmg,
% 6 lpno,
% 7 mcut,
% 8 pcg,
% 9 qual.
%
% 20110628 JLC generalized this to use all of the unwrapping algorithms.

% PROGRAMS=uwdiff flyn fmg gold jump lpno mcut pcg qual unmg unwt

PROGRAMS = {'gold','unwt','unmg','flyn','fmg','lpno','mcut','pcg','qual'};
CODEDIR = '/home/jlc/LAB'; 

uwphase = [];

if(nargin<2)
    ALGO = 'gold';
else
    if(ischar(ALGO))
        if(strcmp(ALGO,'help'))
            fprintf('Algorithms available: (use number or name)\n');
            for n=1:length(PROGRAMS)
                fprintf('\t%2d: %s\n',n,PROGRAMS{n});
            end
            return;
        end
    else
        ALGO = PROGRAMS{ALGO};
    end
end

% if(nargin>1)
%     USE_FFT = true;
% else
%     USE_FFT = false;
% end

[nx,ny] = size(phase);

% if(USE_FFT)
    NX = 2^floor(log2(nx)+1)+1;
    NY = 2^floor(log2(ny)+1)+1;
%     NX = 2^floor(log2(nx)+1);
%     NY = 2^floor(log2(ny)+1);
    phase = padarray(phase,[(NX-nx) (NY-ny)],'post');
% else
%     NX = nx;
%     NY = ny;
% end

%TEMPDIR = '/mnt/ramdisk';
TEMPDIR = '/tmp';
RAND = floor(1000*rand);

RAW = sprintf('%s/_uwrap_%dx%d_%010d.phase',TEMPDIR,NX,NY,RAND);
RESULT = sprintf('%s/_uwrap_%dx%d_%010d.uwphase',TEMPDIR,NX,NY,RAND);

fid = fopen(RAW,'wb');
fwrite(fid, phase(:), 'float');
fclose(fid);

% fprintf('/home/jlc/LAB/gold -input pengaha1.phase -format float -xsize 513 -ysize 513 -dipole no -output pengahaout1')
% system('/home/jlc/LAB/gold -input pengaha1.phase -format float -xsize 513 -ysize 513 -dipole no -output pengahaout1')


% if(USE_FFT)
%     CMD = sprintf('/home/jlc/LAB/unwt -input "%s" -format float -xsize %d -ysize %d -output "%s"  > uwrap.out 2> uwrap.err',RAW,NX,NY,RESULT);
% else
%     CMD = sprintf('/home/jlc/LAB/gold -input "%s" -format float -xsize %d -ysize %d -output "%s" > uwrap.out 2> uwrap.err',RAW,NX,NY,RESULT);
% end

CMD = sprintf('%s/%s -input "%s" -format float -xsize %d -ysize %d -output "%s"  > /tmp/uwrap.out 2> /tmp/uwrap.err',...
    CODEDIR,ALGO,RAW,NX,NY,RESULT);

% fprintf('CMD: %s\n',CMD);
[ret,cmdout] = unix(CMD);

% if(ret ~= 0)
%     fprintf('COMMAND: %s\nOUTPUT: %s\n',CMD,cmdout);
%     !cat /tmp/uwrap.err
%     %error('The system call failed.');
% end

[fid,errmsg] = fopen(RESULT,'r');
if(fid<0)
    fprintf('ERROR: %s --> %s\n',RESULT,errmsg);
    fprintf('Hint: make sure the ~/LAB/* routines exist from Gigglia & Pritt.\n');
end

uwphase = fread(fid,NX*NY,'float');
uwphase = reshape(uwphase,NX,NY);
uwphase = uwphase(1:nx,1:ny);
fclose(fid);

delete(RAW);
delete(RESULT);

return;
