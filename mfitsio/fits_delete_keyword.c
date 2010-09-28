/**
 * File:        fits_delete_keyword.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The gateway program for fits_delete_keyword.
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
 * Executes the fits_delete_keyword function. Type "help fits_delete_keyword"
 * for more information.
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
  char *filename, *keyword;
  mfitsio_header *header;

  /* Check proper input and output */
  if (nrhs != 2)
  {
    mexErrMsgTxt ("Two inputs required.");
  }
  else if (nlhs > 1)
  {
    mexErrMsgTxt ("Too many output arguments.");
  }
  else if (!mxIsChar (prhs[0]) || !mxIsChar (prhs[1]))
  {
    mexErrMsgTxt ("Inputs must be strings (filename, keyword).");
  }
  else if (mxGetM (prhs[0]) != 1 || mxGetM (prhs[1]) != 1)
  {
    mexErrMsgTxt ("Inputs must be a row vector strings (filename, keyword).");
  }
  buflen = (mxGetM (prhs[0]) * mxGetN (prhs[0])) + 1;
  filename = (char *) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[0], filename, buflen);

  buflen = (mxGetM (prhs[1]) * mxGetN (prhs[1])) + 1;
  keyword = (char *) mxCalloc (buflen, sizeof (char));
  status = mxGetString (prhs[1], keyword, buflen);


  mfitsio_delete_keyword (filename, keyword);
  return;
}
