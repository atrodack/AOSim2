% Function Name:
%    fits_write_image
%
% Description: Writes an image and a header to a FITS file.
%
% Usage:
%    fits_write_image(FILENAME, IMAGE);
%    fits_write_image(FILENAME, IMAGE, HEADER);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%    IMAGE: An image array representing the filename.
%    HEADER: Optional. A structure array containing keywords to be modified.
%
% Returns:
%    Nothing.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function fits_write_image(FILENAME, IMAGE, HEADER);

% Author: Damian Eads <eads@lanl.gov>
%
% This software and ancillary information (herein called ``Software'')
% called MFITSIO is made available under the terms described here. The
% SOFTWARE has been approved for release with associated LA-CC number:
% LA-CC-02-085.
% 
% This SOFTWARE has been authored by an employee or employees of the
% University of California, operator of the Los Alamos National Laboratory
% under Contract No. W-7405-ENG-36 with the U.S. Department of Energy. 
% The U.S. Government has rights to use, reproduce, and distribute this
% SOFTWARE.  The public may copy, distribute, prepare derivative works and
% publicly display this SOFTWARE without charge, provided that this Notice
% and any statement of authorship are reproduced on all copies.  Neither
% the Government nor the University makes any warranty, express or
% implied, or assumes any liability or responsibility for the use of this
% SOFTWARE.  If SOFTWARE is modified to produce derivative works, such
% modified SOFTWARE should be clearly marked, so as not to confuse it with
% the version available from LANL.

