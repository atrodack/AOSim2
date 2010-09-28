/**
 * File:        mfitsio.h
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

#ifndef MFITSIO_H

/** Ensures prototypes are processed only once. **/
#define MFITSIO_H

/** Macro specifies which function to execute when a warning is reported. */
#define MFITSIO_WARN mexWarnMsgTxt

/** Macro specifies which function to execute when an error is reported. */
#define MFITSIO_ERR mexErrMsgTxt

/** Macro specifies which function to execute when writing standard text. */
#define MFITSIO_PRINTF mexPrintf

/** Macro specifies which function to execute when freeing memory. */
#define MFITSIO_FREE free

/** Macro specifies which function to execute when allocating memory using
    the malloc invocation syntax. */
#define MFITSIO_MALLOC malloc

/** Macro specifies which function to execute when allocating memory using
    the calloc invocation syntax. */
#define MFITSIO_CALLOC calloc

/** Macro specifies which status flag should be used for functions which
    only read to a FITS image. */

#define MFITSIO_READONLY READONLY

/** Macro specifies which status flag should be used for functions which
    only read and write to a FITS image. */

#define MFITSIO_READWRITE READWRITE

#include "mex.h"
#include "matrix.h"
#include <fitsio.h>

/** This structure holds keyword information extracted from a FITS file
    header. */

typedef struct
{
  char *key;
  char *value;
  int type;
}
mfitsio_record;

/** This structure stores an entire extracted header. */

typedef struct
{
  mfitsio_record **records;
  int length;
}
mfitsio_header;

/** This structure stores information about a FITS data array or a MATLAB
    image. */

typedef struct
{
  int naxis;
  int *naxes;
  int bitpix;
  int check;
}
mfitsio_info;

void mfitsio_free_record (mfitsio_record * record);
void mfitsio_free_header (mfitsio_header * header);
void mfitsio_free_info (mfitsio_info * info);
mfitsio_record *mfitsio_parse_card (const char *cardtext);
int mfitsio_ignore_card (const char *s);
mfitsio_header *mfitsio_read_header (const char *filename);
mxArray *mfitsio_adapt_frecord (const mfitsio_record * record);
mxArray *mfitsio_adapt_fheader (const mfitsio_header * record);
mfitsio_info *mfitsio_read_info (const char * filename);
mxArray *mfitsio_read_image (const char *filename,
			      const mfitsio_info * info);
mxArray *mfitsio_read_image_impl (const char *filename,
				  const mfitsio_info * info,
				  long *fpixels, long *lpixels,
				  long *lnaxes, long *inc);
void mfitsio_write_image (const char *filename,
			   const mxArray * header, mxArray * img);
void mfitsio_write_image_impl (const char *filename, const mxArray * header,
			       mxArray * img, long *fpixels, long *lpixels,
			       long *lnaxes, long *inc, mfitsio_info *info);
void mfitsio_write_header (const char *filename, const mxArray * header);
void mfitsio_delete_keyword (const char *filename, const char *keyword);
mfitsio_info *mfitsio_calc_info (const mxArray * img);
void mfitsio_write_info (fitsfile * fptr, const mfitsio_info * info);
int mfitsio_is_scalar (const mxArray * array);
int mfitsio_forbidden (const char *s);
mxArray *mfitsio_get_mfield (const mxArray * array, int index,
			     const char *field_name);
double mfitsio_get_mscalar (const mxArray * array);
void mfitsio_check_coordinate (const mxArray *crd1, const mxArray *crd2,
			       const int naxis, const int *naxes);
void mfitsio_check_size_info(const mxArray *crd1);
long *mfitsio_convert_dbl2long (const double *dbl, const int size);
long *mfitsio_convert_int2long (const int *in, const int size);
long *mfitsio_create_ones_vector (const int size);
int *mfitsio_get_region_size(const long *crd1, const long *crd2,
			     const int naxis);
mxArray *mfitsio_get_lpixels(const mxArray *spixels, const mxArray *img, 
			     int naxis);
void mfitsio_write_image_execute(fitsfile *file,
				 mxArray *img, long group, long naxis,
				 long *naxes, long *fpixel, long *lpixel,
				 int *status);
mxArray *mfitsio_read_image_execute(const mfitsio_info *info, int *diff,
				    fitsfile *fptr,
				    int group, int naxis, long *lnaxes,
				    long *fpixel, long *lpixel,
				    long *inc, int *anynul, int *status);
#endif /** MFITSIO_H **/
