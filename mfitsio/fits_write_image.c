/**
 * File:        fits_write_image.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_write_image.
 * Created:     December 2, 2002
 *
 * This software and ancillary information (herein called ``Software'')
 * called MFITSIO is made available under the terms described here. The
 * SOFTWARE has been approved for release with associated LA-CC number:
 * LA-CC-02-085.
 *
 * This SOFTWARE has been authored by an employee or employees of the
 * University of California, operator of the Los Alamos National Laboratory
 * under Contract No. W-7405-ENG-36 with the U.S. Department of Energy. 
 * The U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE.  The public may copy, distribute, prepare derivative works and
 * publicly display this SOFTWARE without charge, provided that this Notice
 * and any statement of authorship are reproduced on all copies.  Neither
 * the Government nor the University makes any warranty, express or
 * implied, or assumes any liability or responsibility for the use of this
 * SOFTWARE.  If SOFTWARE is modified to produce derivative works, such
 * modified SOFTWARE should be clearly marked, so as not to confuse it with
 * the version available from LANL.
 *
 * Consult the COPYING file for more details on licensing.
 */

#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <malloc.h>

#include "mex.h"
#include "matrix.h"
#include "mfitsio.h"

/**
 * Executes the fits_write_image function. Type "help fits_write_image" for
 * more information.
 *
 * @param nlhs             The total number of output arguments.
 * @param plhs             The output arguments.
 * @param nrhs             The total number of input arguments.
 * @param prhs             The input arguments.
 */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int buflen, status;
  char *filename;
  mfitsio_header *header;
  mfitsio_info *info;
  mxArray *headerm = 0;
  mxArray *img = 0;

  /* Check proper input and output */
  if (nrhs < 2)
    mexErrMsgTxt ("At least two inputs required.");
  else if (nlhs > 0)
    mexErrMsgTxt ("Too many output arguments.");
  else if (!mxIsChar (prhs[0]))
    mexErrMsgTxt ("Input must be a string (filename).");

  if (mxGetM (prhs[0]) != 1)
    mexErrMsgTxt ("Input must be a row vector string (filename).");

  buflen = (mxGetM (prhs[0]) * mxGetN (prhs[0])) + 1;

  filename = (char*) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);
  img = (mxArray *) prhs[1];

  /**  if (mxGetClassID(img) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("The image should be a M x N x P x ... double array. Use "
                 "double(IMG) to cast.");
		 }**/
  if (nrhs == 3)
  {
    headerm = (mxArray *) prhs[2];
    if (mxGetClassID (headerm) != mxSTRUCT_CLASS)
    {
      mexErrMsgTxt ("The header must be a 1x1 struct array.");
    }
  }
  else
  {
    int dims[] = { 1, 1 };
    int ndims = 2;
    int nfields = 1;
    const char **fields;
    fields = (const char **) malloc (sizeof (char *) * nfields);
    fields[0] = "_nil_";
    headerm = mxCreateStructArray (ndims, dims, nfields, fields);
  }
  mfitsio_write_image (filename, headerm, img);

  return;
}

