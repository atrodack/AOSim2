% Function Name:
%    fits_read_header
%
% Description: Reads a FITS file and stores all header information in
% a structure array.
%
% Usage:
%    HEADER = fits_read_header(FILENAME);
%
% Arguments:
%    FILENAME: A character array representing the filename.
%
% Returns:
%    HEADER:   The fits header as a structure array.
%
% Type 'mfitsio_license' to display the MFITSIO licensing agreement.

function HEADER = fits_read_header(FILENAME);

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
%
% Consult the COPYING file for more details on licensing.

