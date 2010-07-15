function fitswritesimple(data,filename)
%FITSWRITE saves a Matlab matrix as a FITS image
%
%Usage: fitswrite(data, filename)
%
%Known deficiencies
%
%1. Doesn't support arrays of imaginary numbers.
%2. Only handles simple 2 dimensional arrays of data.
%
%Author: R. G. Abraham, Institute of Astronomy, Cambridge University
%        abraham@ast.cam.ac.uk

%retrieved from: http://physics.oulu.fi/ftp/personal/zieger/KittPeak/fitswrite.m
%UPDATES: 3/11/10 Phil Hinz, 
%--added check for MAC computer
%--changed output data from float to double
%--added fliplr and rot90 to data manipulation; interchanged nrow and ncol

                             


[nrow,ncol]=size(data);
header_cards = [make_card('SIMPLE','T');       ...
      make_card('BITPIX',-64);       ...
      make_card('NAXIS',2);          ...
      make_card('NAXIS1',ncol);      ...
      make_card('NAXIS2',nrow);      ...
      make_card('BSCALE',1.0);       ...
      make_card('BZERO',0.0);        ...
      make_card('CREATOR','Matlab'); ...
      make_card('DATE',date);        ...
      make_card('END')];

header_record = make_header_record(header_cards);
%[ncards,dummy]=size(header_cards);
%fprintf(header_record(1,:));

fid=fopen(filename,'w');
fwrite(fid,header_record','char');

% try to figure out if we need to swap bytes. This is
% imperfect as I don't know the endian-ness of each
% architecture, so I'm only looking for ones I know for 
% sure are big-endian.
friend = computer;
if strmatch(friend,'PCWIN')
   bswap = 'b';
elseif strmatch(friend,'LNX86')
   bswap = 'b';   
elseif strmatch(friend,'ALPHA')
   bswap = 'b';
elseif strmatch(friend,'MACI')
   bswap = 'b';
else
   bswap = 'b';
end

%PMH: flip FITS file
data=fliplr(data);
data=rot90(data);


fwrite(fid,data,'double',bswap);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function card=make_card(keyword,value)
%MAKE_CARD turns a set of strings into a valid FITS card

%Make keyword field 8 bytes long
lk=length(keyword);
if (lk > 8) & (nargin>1)
	error('Keyword must be less than or equal to 8 characters!')
elseif (lk < 8 )
	keyword=[keyword,setstr(ones(1,8-lk)*32)];
end;

%Deal with both straight keyword and keyword/value pair
if (nargin==1)
	%Keyword without a value
	card=keyword;	
else
	%Key/value pair has an equal sign and space at bytes 9 and 10
	card=[keyword,'= '];

	%Now output the value. The FITS standard wants things to start 
	%in different columns depending on what type of data the
	%value holds, according to the following rules:
	%
	%  Logical: T or F in column 30
	%
	%  Character string: A beginning quote in column 11 and an
	%  ending quote between columns 20 and 80.
	%
	%  Real part of an integer or floating point number: right 
	%  justified, ending in column 30.
	%
	%  Imaginary part: right justified, ending in
	%  column 50, starting after column 30 (NB. I won't bother 
	%  to support an imaginary part in this M-file, and will 
	%  let some radio astronomer who needs it add it if they want).

	if isstr(value)
  	    %Test for logical. If logical it goes in column 30 
		if (length(value)==1) & (strmatch(upper(value),'T') | strmatch(upper(value),'F'))
 			card=[card,setstr(ones(1,19)*32),value];	
		else	
			%Value must be a character string. Pad if less than 8
			%characters long.
			lv=length(value);
		    if (lv > 70)
		   error('Value must be less than 70 characters long!')
		    elseif (lv < 10 )
 	 	   value=[value,setstr(ones(1,8-lv)*32)];
		    end;
			card=[card,'''',value,''''];
		end;	
	else
		%Value must be a number. Convert to a string. Maximum
		%precision is set to 10 digits
		value=num2str(value,10);
		lv=length(value);
	
		%Write this out so it will end on column 30
		card=[card,setstr(ones(1,20-lv)*32),value];	
	end;
end;

%Now pad the output to make it exactly 80 bytes long
card=[card,setstr(ones(1,80-length(card))*32)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hrec=make_header_record(card_matrix)

[nrow,ncol] = size(card_matrix);
n_blanks = 36 - rem(nrow,36);
blank_line = setstr(ones(1,80)*32);
hrec = [card_matrix; repmat(blank_line,n_blanks,1)];
