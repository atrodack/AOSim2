/**
 * File:        mfitsio.c
 * Author:      Damian Ryan Eads <eads@lanl.gov>
 * Description: The main "work horse" functions for the MATLAB-fitsio interface
 *              are implemented in this file.
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

/** Include necessary headers from the C standard library. */
#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <malloc.h>

#if !defined(_MSC_VER)
#include <unistd.h>
#endif

/** Include mex related headers. */
#include "mex.h"
#include "matrix.h"

/** Include local headers. */
#include "mfitsio.h"

#define F_OK 0

/*************************************************************************
 *  Allocation and deallocation functions for mfitsio data structures. *
 *************************************************************************/

/**
 * Deallocates a mfitsio_record structure.
 *
 * @param  record       The record you wish to deallocate.
 */

void
mfitsio_free_record (mfitsio_record * record)
{
  MFITSIO_FREE (record->key);
  MFITSIO_FREE (record->value);
}

/**
 * Deallocates a mfitsio_header structure.
 *
 * @param  header       The record you wish to deallocate.
 */

void
mfitsio_free_header (mfitsio_header * header)
{
  int i;
  for (i = 0; i < header->length; i++)
  {
    mfitsio_free_record (header->records[i]);
  }
  MFITSIO_FREE (header->records);
}

/**
 * Deallocates a mfitsio_info structure.
 *
 * @param  info         The info you wish to deallocate.
 */

void
mfitsio_free_info (mfitsio_info * info)
{
  MFITSIO_FREE (info->naxes);
  MFITSIO_FREE (info);
}

/**
 * Parses a fits header "card". Comments are ignored. A mfitsio_record
 * structure is returned. When this structure is no longer needed, it should
 * be deallocated.
 *
 * @param  cardtext     The text of the fits header card.
 *
 * @return              A mfitsio_record structure containing the key
 *                      and value which was parsed from the card. NULL is
 *                      returned if the field is to be ignored or if the card
 *                      is a comment line.
 */

mfitsio_record *
mfitsio_parse_card (const char *cardtext)
{
  int i = 0, j = 0, k = 0, len_card = 0, len_key = 0, len_value = 0;
  char value = 0;
  mfitsio_record *record = 0;

  /** If the field is to be ignored, then do not bother parsing it. */
  if (!mfitsio_ignore_card (cardtext))
  {
    /** Allocate memory for the record. */
    record = (mfitsio_record *) MFITSIO_MALLOC (sizeof (mfitsio_record));
    record->key = 0;
    record->value = 0;

    /** Acquire the length of the card string. */
    len_card = strlen (cardtext);

    /** Search for the terminating space or equal sign. This will allow us
     to determine the memory needed for the key string. */
    for (i = 0; i < len_card; i++)
    {
      /** When a space is encountered, break out. */
      if (cardtext[i] == ' ')
      {
	break;
      }
      /** When an equal sign is encountered, break out. Set the marker
          variable k to true so that we don't further remove any more
          equal signs for this card. */
      else if (cardtext[i] == '=')
      {
	k = 1;
	break;
      }
    }

    /** Allocate memory for the key. */
    len_key = i + 1;
    record->key = (char *) MFITSIO_MALLOC (len_key * sizeof (char));

    /** Copy the key text into the buffer. */
    memcpy (record->key, cardtext, i);

    /** Place a null terminator. */
    record->key[i] = 0;

    /** Start after the key text. */
    j = i + 1;

    /** If we did not remove the equal sign separating the key and value,
	do so now. */
    if (k == 0)
    {
      for (; j < len_card; j++)
      {
	if (cardtext[j] == '=')
	{
	  j++;
	  break;
	}
      }
    }

    /** Trim surrounding spaces. */
    for (; j < len_card; j++)
    {
      if (cardtext[j] != ' ')
      {
	break;
      }
    }

    /** Determine how much memory is needed for the value. */
    len_value = (len_card - j + 1);

    /** Allocate memory for the value of the keyword. */
    record->value = (char *) MFITSIO_MALLOC (len_value * sizeof (char));
    memcpy (record->value, cardtext + j, len_value - 1);

    /** Place a null terminator. */
    record->value[len_value - 1] = 0;

    /** Remove the comment from the value string simply by placing
	null terminators.
        Note: This is a kludge. In the future, only the memory needed
	should be allocated. */

    for (j = len_value - 1; j >= 0; j--)
    {
      value = record->value[j];
      if (value == '\'')
      {
	break;
      }
      else if (value == '/' || value == ' ')
      {
	record->value[j] = '\0';
      };
    }
  }

  /** Return the result of the parsing. NULL is returned if the card is to
      be ignored. */

  return record;
}

/**
 * Reads a header from a fits file. Before the result of this function can be
 * used by MATLAB, the header must be adapted to a MATLAB struct array using
 * mfitsio_adapt_fheader function.
 *
 * @param filename      The filename of the fits file to read.
 *
 * @return              A mfitsio_header structure acquired by reading
 *                      the header in the fits file.
 */

mfitsio_header *
mfitsio_read_header (const char *filename)
{
  fitsfile *fptr = 0;
  mfitsio_header *header = 0, *new_header = 0;
  char card[FLEN_CARD];
  int status = 0, nkeys = 0, jkeys = 0, i = 0, j = 0;
  mfitsio_record *record = 0;

   /** Open the fits file. */
  fits_open_file (&fptr, filename, MFITSIO_READONLY, &status);
  if (status)
  {
    MFITSIO_PRINTF ("1  An error occurred: %d\n", status);
    MFITSIO_ERR ("Unable to open FITS file.");
  }

   /** Retrieve the number of header keywords stored. */
  fits_get_hdrspace (fptr, &nkeys, NULL, &status);
  if (status)
  {
    MFITSIO_PRINTF ("2  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not access header.");
  }

   /** Allocate memory for the header structure. */
  header = (mfitsio_header *) MFITSIO_MALLOC (sizeof (mfitsio_header));

   /** Set the length of the header to the number of keys. */
  header->length = nkeys;

   /** Allocate memory for the records. */
  header->records =
    (mfitsio_record **) MFITSIO_MALLOC (nkeys * sizeof (mfitsio_record *));

   /** For each line in the header, read it and parse it. Place the result
       in the header structure. */
  for (i = 1; i <= nkeys; i++)
  {
    fits_read_record (fptr, i, card, &status);
    if (status)
    {
      MFITSIO_PRINTF ("3  An error occurred: %d\n", status);
      MFITSIO_ERR ("Could not access header.");
    }
    record = mfitsio_parse_card (card);
    header->records[i - 1] = record;

     /** If the line is not ignored, increase the number of read lines. */
    if (record != 0)
    {
      jkeys++;
    }
  }

   /** If there exists fields we wish to ignore, filter them out into
       a new header structure. **/

  if (jkeys != nkeys)
  {
     /** Allocate memory for the new header. */
    new_header =
      (mfitsio_header *) MFITSIO_MALLOC (sizeof (mfitsio_header));
    new_header->records =
      (mfitsio_record **) MFITSIO_MALLOC (jkeys *
					    sizeof (mfitsio_record *));

     /** Store the unignored records in the new structure. */
    new_header->length = jkeys;
    j = 0;
    for (i = 0; i < nkeys; i++)
    {
      if (header->records[i] != 0)
      {
	new_header->records[j++] = header->records[i];
      }
    }
     /** Deallocate the old header structure. */
    MFITSIO_FREE (header->records);
    MFITSIO_FREE (header);

     /** Set the header pointer to the new header. */
    header = new_header;
  }

  /** Finally, close the file. */
  fits_close_file (fptr, &status);
  if (status)
  {
    MFITSIO_PRINTF ("4  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not close file.");
  }

  return (status != 0) ? 0 : header;
}

/**
 * Adapts a mfitsio_record structure into a MATLAB variable. The syntax of the
 * text of the value determines the type of the MATLAB variable.
 *
 * @param record        The record to adapt.
 *
 * @return              A MATLAB variable.
 */

mxArray *
mfitsio_adapt_frecord (const mfitsio_record * record)
{
  int len = 0, i, j;
  char *temp = 0;
  char f = 0;

   /** Variables to store arguments for calls to MATLAB functions. */
  mxArray *retval = 0;
  mxArray *iargs[] = { 0 };
  mxArray *oargs[] = { 0 };

   /** Acquire the length of the keyword text. */
  len = strlen (record->value);

   /** If the field is not empty, parse it. **/
  if (len)
  {

     /** Acquire the first character of the keyword text. */
    f = record->value[0];

     /** If a string is enclosed within single quotes, create a MATLAB
	 string with the single quotes removed. */
    if (f == '\'')
    {
       /** Allocate memory for the temporary string. */
      temp = (char *) MFITSIO_MALLOC ((len - 1) * sizeof (char));

       /** Copy the text past the single quote. */

      memcpy (temp, record->value + 1, len - 2);

       /** Assign a null terminator at the spot of the terminating quotation.*/
      temp[len - 2] = 0;

       /** This variable determines where the next null terminator should
	   be placed. -1 means don't place another null terminator. */
      i = -1;
      for (j = strlen (temp) - 1; j >= 0; j--)
      {
	if (temp[j] != ' ')
	{
	  i = j + 1;
	  break;
	}
      }
      if (i != -1)
      {
	temp[i] = '\0';
      }

       /** Create a MATLAB string to hold the text. */
      retval = mxCreateString (temp);
    }
     /** If the value is a logical, create a MATLAB logical variable. */
    else if (f == 'T' || f == 'F' || f == 't' || f == 'f')
    {
      retval = mxCreateLogicalScalar (f == 'T' || f == 't');
    }
     /** If the value is a number, create an integer or a double. */
    else if ((f >= '0' && f <= '9') || f == '-' || f == '.')
    {
       /** If a . exists in the string, then the variable is assumed to be
	   a double. */
      if (strchr (record->value, '.'))
      {
	retval = mxCreateDoubleScalar (atof (record->value));
      }
       /** Otherwise it is an integer. A MATLAB int32 object is created. **/
      else
      {
	 /** The value is parsed. */
	retval = mxCreateDoubleScalar (atof (record->value));
	iargs[0] = retval;

	 /** And then converted to an integer. */
	mexCallMATLAB (1, oargs, 1, iargs, "int32");
	retval = oargs[0];
      }
    }
     /** If the syntax is incorrect, we will convert the text to a string. */
    else
    {
      retval = mxCreateString (strdup (record->value));
      i = -1;
      for (j = strlen (record->value) - 1; j >= 0; j--)
      {
	if (record->value[j] != ' ')
	{
	  i = j + 1;
	  break;
	}
      }
      if (i != -1)
      {
	record->value[i] = '\0';
      }
    }
  }
   /** If the text was blank, a double value of zero will be returned. */
  else
  {
    retval = mxCreateDoubleScalar (0);
  }
  return retval;
}

/**
 * Adapts a mfitsio_header structure into a MATLAB struct array. 
 *
 * @param header        The header to adapt.
 *
 * @return              A MATLAB struct object corresponding to the parsed
 *                      mfitsio_header.
 */

mxArray *
mfitsio_adapt_fheader (const mfitsio_header * header)
{
  const char **fnames = 0;
  int buflen, status = 0, i;
  char *temp = 0;
  mxArray *retval = 0;

  if (header == 0)
  {
    MFITSIO_ERR ("33 Header unreadable.");
  }
  /** The number of field names is equal to the length of the header. */
  fnames = (const char **)MFITSIO_CALLOC (header->length, sizeof (*fnames));

  /** Go through each record in the header and copy the text of the keys
      into the array. */
  for (i = 0; i < header->length; i++)
  {
    /** Acquire the buffer length. */
    buflen = strlen (header->records[i]->key);

    /** Allocate memory for the buffer. */
    temp = (char *) MFITSIO_MALLOC (sizeof (char) * buflen + 1);

    /** Copy the field value into the new buffer. */
    memcpy (temp, header->records[i]->key, buflen);

    /** Place a null terminator. */
    temp[buflen] = 0;

    /** Store the buffer in the array of keys. */
    fnames[i] = temp;
  }

  /** Create a MATLAB struct array to store the header information. */
  retval = mxCreateStructMatrix (1, 1, header->length, fnames);

  /** For each field of the MATLAB struct, set the value. Since the record
      text is a character string, its syntax must be analyzed before assigning
      a type. cfitiosm_adapt_frecord will do all the dirty work for us. */

  for (i = 0; i < header->length; i++)
  {
    mxSetFieldByNumber (retval, 0, i,
			mfitsio_adapt_frecord (header->records[i]));
  }
  return retval;
}

/**
 * The dimensions and bitpix are acquired from the FITS file header.
 *
 * @param filename      The filename of the FITS file.
 *
 * @return              A mfitsio_info struct containing the information.
 */

mfitsio_info *
mfitsio_read_info (const char * filename)
{
  mfitsio_info *retval = 0;
  fitsfile *fptr = 0;
  long *naxes;
  int status = 0, i;

  /** Open the file. */
  fits_open_file (&fptr, filename, MFITSIO_READONLY, &status);
  if (status)
  {
    MFITSIO_PRINTF ("5  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not open file for reading.");
  }

  /** Allocate memory for the info structure. */
  retval = (mfitsio_info *) MFITSIO_MALLOC (sizeof (mfitsio_info));

  /** Store the important fields we really care about in an info structure. */
  fits_get_img_type(fptr, &(retval->bitpix), &status);
  if (status)
  {
    MFITSIO_PRINTF ("6  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not retrieve BITPIX keyword.");
  }

  /** Get the number of dimensions. */
  fits_get_img_dim(fptr, &(retval->naxis), &status);
  if (status)
  {
    MFITSIO_PRINTF ("7  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not retrieve NAXIS, total number of dimensions.");
  }

  /** Create a temporary array to store the image dimensions as a long. */
  naxes = (long *) MFITSIO_MALLOC (sizeof (long) * retval->naxis);

  /** Get the image dimensions. */
  fits_get_img_size(fptr, retval->naxis, naxes, &status);
  if (status)
  {
    MFITSIO_PRINTF ("8  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not open file for reading.");
  }
  retval->naxes = (int *) MFITSIO_MALLOC (retval->naxis * sizeof (int));

  /** Retrieve the values for each dimension. */
  for (i = 1; i <= retval->naxis; i++)
  {
    retval->naxes[i - 1] = (int) naxes[i - 1];
  }

  /** Free the temporary array. */
  MFITSIO_FREE ( naxes );

  /** Finally, close the file. */
  fits_close_file (fptr, &status);
  if (status)
  {
    MFITSIO_PRINTF ("4  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not close file.");
  }

  /** Return the result of the extraction. */
  return retval;
}

/**
 * Read a FITS image and store the result in a MATLAB array.
 *
 * @param filename      The filename of the FITS image.
 * @param info          The dimensions and bitpix of the image.
 *
 * @return              A mfitsio_info struct containing the information.
 */

mxArray *
mfitsio_read_image (const char *filename, const mfitsio_info * info) {
  mxArray *image = 0;
  long *fpixels = 0, *lpixels = 0, *lnaxes = 0, *inc = 0;
  int i = 0;
  fpixels = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  lpixels = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  lnaxes = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  inc = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  /** Set all the parameters for each dimension. */
  for (i = 0; i < info->naxis; i++)
    {
      lpixels[i] = info->naxes[i];
      fpixels[i] = 1;
      inc[i] = 1;
      lnaxes[i] = info->naxes[i];
    }
  image = mfitsio_read_image_impl(filename, info, fpixels,
				  lpixels, lnaxes, inc);
  MFITSIO_FREE (fpixels);
  MFITSIO_FREE (lpixels);
  MFITSIO_FREE (inc);
  MFITSIO_FREE (lnaxes);
  return image;
}

/**
 * Read a FITS image and store the result in a MATLAB array.
 *
 * @param filename      The filename of the FITS image.
 * @param info          The dimensions and bitpix of the image.
 *
 * @return              A mfitsio_info struct containing the information.
 */

mxArray *
mfitsio_read_image_impl (const char *filename, const mfitsio_info * info,
			 long *fpixels, long *lpixels,
			 long *lnaxes, long *inc)
{
  fitsfile *fptr = 0;
  mxArray *img = 0;
  double *pr = 0;
  int status = 0, i = 0;
  int *diff = 0;

  diff = mfitsio_get_region_size(fpixels, lpixels, info->naxis);
  /** Allocate memory to store the starting and ending coordinates. */

  /** Open the file. */
  fits_open_file (&fptr, filename, MFITSIO_READONLY, &status);
  if (status)
  {
    MFITSIO_PRINTF ("9  An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not open file for reading.");
  }

  /** Allocate memory for the image. */
  img = mfitsio_read_image_execute(info, diff, fptr, 0, info->naxis, lnaxes,
				   fpixels, lpixels, inc, 0, &status);

  MFITSIO_FREE (diff);

  /** Close the FITS file. */
  fits_close_file (fptr, &status);
  if (status)
  {
    MFITSIO_PRINTF ("12 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not close file.");
  }

  /** Return the resultant image. */
  return img;
}

mxArray *mfitsio_read_image_execute(const mfitsio_info *info, int *diff,
				    fitsfile *fptr,
				    int group, int naxis, long *lnaxes,
				    long *fpixel, long *lpixel,
				    long *inc, int *anynul, int *status) {
  mxArray *img;
  void *data;
  switch (info->bitpix)
    {
    case BYTE_IMG:
      img = mxCreateNumericArray (info->naxis, diff, mxUINT8_CLASS, mxREAL);
      data = mxGetData(img);
      fits_read_subset_byt (fptr, group, info->naxis, lnaxes, fpixel, lpixel,
			    inc, 0, (unsigned char*)data, anynul, status);
      break;
    case SHORT_IMG:
      img = mxCreateNumericArray (info->naxis, diff, mxINT16_CLASS, mxREAL);
      data = mxGetData(img);
      fits_read_subset_sht (fptr, 0, info->naxis, lnaxes, fpixel, lpixel,
			    inc, 0, (short*)data, 0, status);
      break;
    case LONG_IMG:
      img = mxCreateNumericArray (info->naxis, diff, mxINT32_CLASS, mxREAL);
      data = mxGetData(img);
      fits_read_subset_int (fptr, 0, info->naxis, lnaxes, fpixel, lpixel,
			    inc, 0, (int*)data, 0, status);
      break;
    case FLOAT_IMG:
      img = mxCreateNumericArray (info->naxis, diff, mxSINGLE_CLASS, mxREAL);
      data = mxGetData(img);
      fits_read_subset_flt (fptr, 0, info->naxis, lnaxes, fpixel, lpixel,
			    inc, 0, (float*)data, 0, status);
      break;
    case DOUBLE_IMG:
      img = mxCreateNumericArray (info->naxis, diff, mxDOUBLE_CLASS, mxREAL);
      data = mxGetData(img);
      fits_read_subset_dbl (fptr, 0, info->naxis, lnaxes, fpixel, lpixel,
			    inc, 0, (double*)data, 0, status);
      break;
    default:
      MFITSIO_PRINTF ("52 Unsupported input BITPIX: %d\n", info->bitpix);
      MFITSIO_ERR ("Cannot continue.");
    }  
  return img;
}



/**
 * Write a MATLAB array/image to a FITS file.
 *
 * @param filename      The filename of the FITS file.
 * @param info          The MATLAB struct array representing the keywords
 *                      to store.
 * @param img           The MATLAB array/image to store in the FITS file.
 */

void
mfitsio_write_image (const char *filename, const mxArray * header,
		      mxArray * img)
{

  mfitsio_info *info = 0;
  long *fpixels = 0, *lpixels = 0, *lnaxes = 0, *inc = 0;
  int i;
  /** Calculate dimension and bitpix information from the MATLAB image. */
  info = mfitsio_calc_info (img);

  /** Allocate memory for parameters to the fits_write_subset function. */
  fpixels = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  lpixels = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  lnaxes = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  inc = (long *) MFITSIO_MALLOC (sizeof (long) * info->naxis);
  
  /** Set all the parameters for each dimension. */
  for (i = 0; i < info->naxis; i++)
    {
      lpixels[i] = info->naxes[i];
      fpixels[i] = 1;
      inc[i] = 1;
      lnaxes[i] = info->naxes[i];
    }

  mfitsio_write_image_impl(filename, header, img, fpixels,
			   lpixels, lnaxes, inc, info);
  MFITSIO_FREE (fpixels);
  MFITSIO_FREE (lpixels);
  MFITSIO_FREE (inc);
  MFITSIO_FREE (lnaxes);
  mfitsio_free_info(info);
}

/**
 * Write a MATLAB array/image to a FITS file.
 *
 * @param filename      The filename of the FITS file.
 * @param info          The MATLAB struct array representing the keywords
 *                      to store.
 * @param img           The MATLAB array/image to store in the FITS file.
 * @param fpixels       The starting pixel coordinates.
 * @param lpixels       The ending pixel coordinates.
 * @param lnaxes        The
 */

void
mfitsio_write_image_impl (const char *filename, const mxArray * header,
			  mxArray * img, long *fpixels, long *lpixels,
			  long *lnaxes, long *inc, mfitsio_info *info)
{

  fitsfile *fptr = 0;
  mxArray *newImg = 0;
  int status = 0, nkeys, i, j;
  double *pr = 0;

  /** Write the header to the FITS file. */
  mfitsio_write_header (filename, header);

  /** Open the file for writing. */
  fits_open_file (&fptr, filename, MFITSIO_READWRITE, &status);
  if (status)
  {
    MFITSIO_PRINTF ("13 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not open file for writing.");
  }

  /** Transpose the planes. */
  newImg = img;
  /**  if (info->naxis > 0)
  {
    mxArray *iargs[] = { 0 };
    mxArray *oargs[] = { 0 };
    iargs[0] = img;
    mexCallMATLAB (1, oargs, 1, iargs, "fits_transpose_all");
    newImg = oargs[0];
    }**/

  /** for (i = 0; i < info->naxis; i++)
      {**/
    /** MFITSIO_PRINTF("fp: %d lp: %d\n", fpixels[i], lpixels[i]);**/
  /**}**/

  /** Get the data location of the resultant image. */

  /** Resize the old image based on the new bitpix and axes. */
  if (info->check == 0) {
    fits_resize_img (fptr, info->bitpix, info->naxis, lnaxes, &status);
  }
  if (status)
  {
    MFITSIO_PRINTF ("14 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not resize FITS file image.");
  }

  /**
  MFITSIO_PRINTF("BITPIX: %d, NAXIS: %d NAXIS1: %d NAXIS2: %d \n",
                 info->bitpix, info->naxis, lnaxes[0], lnaxes[1] );**/

  /** Write the image to the FITS file.*/
  /**  fits_write_subset_dbl (fptr, 0, info->naxis, lnaxes, fpixels, lpixels,
       pr, &status);**/
  mfitsio_write_image_execute(fptr, img, 0, info->naxis, lnaxes, fpixels,
			      lpixels, &status);

  if (status)
  {
    MFITSIO_PRINTF ("15 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not write the image to the FITS file.");
  }

  /** Close the FITS file. */
  fits_close_file (fptr, &status);
  if (status)
  {
    MFITSIO_PRINTF ("39 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not close file.");
  }
}

void mfitsio_write_image_execute(fitsfile *file,
				 mxArray *img, long group, long naxis,
				 long *naxes, long *fpixel, long *lpixel,
				 int *status) {
  mxClassID classId;
  void *data;
  data = mxGetData(img);
  classId = mxGetClassID (img);
  switch (classId)
    {
    case mxSINGLE_CLASS:
      fits_write_subset_flt(file, group, naxis, naxes, fpixel, lpixel,
			    (float*)data, status);
      break;
    case mxDOUBLE_CLASS:
      /** If double, store the value as a double. */
      fits_write_subset_dbl(file, group, naxis, naxes, fpixel, lpixel,
			    (double*)data, status);
      break;
    case mxUINT8_CLASS:
      /** If an unsigned 8 bit int, store the value as such. */
      fits_write_subset_byt(file, group, naxis, naxes, fpixel, lpixel,
			    (unsigned char*)data, status);
      break;
    case mxUINT16_CLASS:
      /** If an unsigned 16 bit int, store the value as such. */
      fits_write_subset_usht(file, group, naxis, naxes, fpixel, lpixel,
			     (unsigned short int*)data, status);
      break;
    case mxINT16_CLASS:
      /** If a signed 16 bit int, store the value as such. */
      fits_write_subset_sht(file, group, naxis, naxes, fpixel, lpixel,
			    (short int*)data, status);
      break;
    case mxUINT32_CLASS:
      /** If an unsigned 32 bit int, store the value as
	  an unsigned short. */
      fits_write_subset_uint(file, group, naxis, naxes, fpixel, lpixel,
			    (unsigned int*)data, status);
      break;
    case mxINT32_CLASS:
      /** If a signed 32 bit int, store the value as a signed short. */
      fits_write_subset_int(file, group, naxis, naxes, fpixel, lpixel,
			    (int*)data, status);
      break;
    case mxLOGICAL_CLASS:
      /** If a logical, store the value as such. */
      break;
    case mxUINT64_CLASS:
      /** If an unsigned 64 bit int, store the value as an unsigned long.*/
#ifdef HAVE_LONGLONG
      fits_write_subset_lnglng(file, group, naxis, naxes, fpixel, lpixel,
			    (LONGLONG*)data);
			    
      break;
#endif
    case mxINT64_CLASS:
      /** If a signed 64 bit int, store the value as a signed long. */
#ifdef HAVE_LONGLONG
      fits_write_subset_lnglng(file, group, naxis, naxes, fpixel, lpixel,
			    (LONGLONG*)data);
      break;
#endif
    default:
      /** Otherwise, report an error.*/
      MFITSIO_PRINTF
	("51 Unsupported data type for image.\n");
      MFITSIO_ERR ("Cannot continue.");
      break;
    }
}

/**
 * Delete a record from a header in a FITS file.
 *
 * @param filename      The filename of the FITS file.
 * @param keyword       The keyword of the record to delete.
 */

void
mfitsio_delete_keyword (const char *filename, const char *keyword)
{
  fitsfile *fptr = 0;
  int status = 0;
  mfitsio_record *record = 0;
  char *tmp = 0;

  fits_open_file (&fptr, filename, MFITSIO_READWRITE, &status);
  if (status)
  {
    MFITSIO_PRINTF ("16 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not open file for writing.");
  }
  fits_movabs_hdu (fptr, 1, 0, &status);
  if (status)
  {
    MFITSIO_PRINTF ("17 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not access header.");
  }
  if (!mfitsio_forbidden (keyword))
  {
    tmp = strdup (keyword);
    fits_delete_key (fptr, tmp, &status);
    free (tmp);
  }
  else
  {
    MFITSIO_PRINTF ("19 It is forbidden to delete the keyword '%s'\n", keyword);
    MFITSIO_ERR ("Could not delete keyword.");
  }

  if (status)
  {
    MFITSIO_PRINTF ("18 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not delete keyword.");
  }
}

/**
 * Write header information stored in a MATLAB struct array to a FITS file.
 *
 * @param filename      The filename of the FITS file.
 * @param header        The MATLAB struct array containing the information.
 */

void
mfitsio_write_header (const char *filename, const mxArray * header)
{
  fitsfile *fptr = 0;
  int numFields = 0, i = 0, status = 0, singleton = 0, buflen = 0;
  mxClassID classId;
  char *fieldName = 0, *stringval = 0;
  float floatval = 0.0f;
  double dblval = 0.0;
  unsigned char ucharval = 0;
  char charval = 0;
  short shortval = 0;
  unsigned short ushortval = 0;
  int intval = 0;
  unsigned int uintval = 0;
  long longval = 0;
  unsigned long ulongval = 0;
  mxArray *field = 0;

  int naxis = 2;
  int bitpix = 8;
  long naxes[] = { 1, 1, 1 };

  /** If the file does not exist, create it.
      Note to self: access is a function in unistd.h. */
  if (access (filename, F_OK))
  {
    fits_create_file (&fptr, filename, &status);
    if (status)
    {
      MFITSIO_PRINTF ("20 An error occurred: %d\n", status);
      MFITSIO_ERR ("Could not create new FITS file.");
    }
    /** Create an image within the FITS file. */
    fits_create_img (fptr, bitpix, naxis, naxes, &status);
    if (status)
    {
      MFITSIO_PRINTF ("21 An error occurred: %d\n", status);
      MFITSIO_ERR ("Could not create an image within the FITS file.");
    }
    /** Place warning text within the FITS file so to warn other programs
	that the long card convention is being used. */
    fits_write_key_longwarn (fptr, &status);
    if (status)
    {
      MFITSIO_PRINTF ("22 An error occurred: %d\n", status);
      MFITSIO_ERR ("Could not write information to header.");
    }
  }
  else
  {
    /** Open the file for write. */
    fits_open_file (&fptr, filename, MFITSIO_READWRITE, &status);
  }
  if (status)
  {
    MFITSIO_PRINTF ("23 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not open file for writing.");
  }
  /** Go to the HDU. */
  fits_movabs_hdu (fptr, 1, 0, &status);
  if (status)
  {
    MFITSIO_PRINTF ("24 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not access header.");
  }
  /** Retrieve the number of fields from the MATLAB struct array. */
  numFields = mxGetNumberOfFields (header);

  /** For each field in the MATLAB struct array store it in the FITS
      file *only* if it is good thing to do so! */
  for (i = 0; i < numFields; i++)
  {
    /** Retrieve the next field name. */
    fieldName = (char *) mxGetFieldNameByNumber (header, i);

    /** We really don't want users modifying BITPIX and NAXIS keywords
        themselves. This information is automatically calculated whenever
	an image is modified. Whether the field is modifiable is determined
        by mfitsio_forbidden. */
    if (!mfitsio_forbidden (fieldName))
    {
      /** Retrieve the field and class information from the MATLAB variable.*/
      field = mxGetField (header, 0, fieldName);
      classId = mxGetClassID (field);

      /** We do not want to deal with singleton arrays unless they're
          strings.*/
      singleton = mfitsio_is_scalar (field);

      if (singleton || classId == mxCHAR_CLASS)
      {

	/** Depending on the MATLAB class, store the appropriate data type. */
	switch (classId)
	{
	case mxSINGLE_CLASS:
	  /** If single store the value as a float. */
	  floatval = (float) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TFLOAT, fieldName, &floatval, 0, &status);
	  break;
	case mxDOUBLE_CLASS:
	  /** If double store the value as a double. */
	  dblval = (double) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TDOUBLE, fieldName, &dblval, 0, &status);
	  break;
	case mxCHAR_CLASS:
	  /** If a character string, store the value as a string. */
	  buflen = (mxGetM (field) * mxGetN (field)) + 1;
	  stringval = (char *) MFITSIO_CALLOC (buflen, sizeof (char));

	  status = mxGetString (field, stringval, buflen);
	  fits_update_key_longstr (fptr, fieldName, stringval, 0, &status);
	  break;
	case mxUINT8_CLASS:
	  /** If an unsigned 8 bit int, store the value as such. */
	  ucharval = (unsigned char) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TBYTE, fieldName, &ucharval, 0, &status);
	  break;
	case mxINT8_CLASS:
	  /** If a signed 8 bit int, store the value as such. */
	  charval = (char) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TBYTE, fieldName, &charval, 0, &status);
	  break;
	case mxUINT16_CLASS:
	  /** If an unsigned 16 bit int, store the value as
              an unsigned short. */
	  ushortval = (unsigned short) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TUSHORT, fieldName, &ushortval, 0, &status);
	  break;
	case mxINT16_CLASS:
	  /** If a signed 16 bit int, store the value as a signed short. */
	  shortval = (short) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TSHORT, fieldName, &shortval, 0, &status);
	  break;
	case mxUINT32_CLASS:
	  /** If an unsigned 32 bit int, store the value as an unsigned long.*/
	  ulongval = (unsigned int) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TUINT, fieldName, &ulongval, 0, &status);
	  break;
	case mxINT32_CLASS:
	  /** If a signed 32 bit int, store the value as a signed long. */
	  longval = (int) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TINT, fieldName, &longval, 0, &status);
	  break;
	case mxLOGICAL_CLASS:
	  /** If a logical, store the value as such. */
	  intval = (int) mfitsio_get_mscalar (field);
	  fits_update_key (fptr, TLOGICAL, fieldName, &intval, 0, &status);
	  break;
	default:
	  /** Otherwise, report an error.*/
	  MFITSIO_PRINTF
	    ("25 Unsupported data type for header field \'%s\'.\n", fieldName);
	  MFITSIO_ERR ("Cannot continue.");
	  break;
	}
      }
      else
      {
	/** If there exists a non-singleton dimension of a non-char array,
	    report an error. */
	MFITSIO_PRINTF
	  ("26 There exists a non-singleton dimension in the value of \'"
	   "%s\'. A scalar is required. The field will be ignored.\n");
	MFITSIO_WARN ("A header field will not be written.");
      }
      /** If there was an error while writing a key, report it. */
      if (status)
      {
	MFITSIO_PRINTF
	  ("27 An error occurred: %d\n. Could not write key for field "
	   "\'%s\'.", status, fieldName);
	MFITSIO_ERR ("Could not write key.");
      }
    }
  }
  /** Close the file. */
  fits_close_file (fptr, &status);
  if (status)
  {
    MFITSIO_PRINTF ("28 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not create an image within the FITS file.");
  }
}

/**
 * Write bitpix and dimension information to a FITS header.
 *
 * @param filename      The filename of the FITS file.
 * @param info          The mfitsio_info struct containing the dimension
 *                      and bitpix information.
 */

void
mfitsio_write_info (fitsfile * fptr, const mfitsio_info * info)
{
  char naxisf[8];
  int status = 0, i = 0, val = 0;
  /** Because stupid fitsio avoids the use of const, I am forced to copy
      each value into a temporary variable to please the compiler. */

  /** Store the NAXIS and BITPIX information in the FITS file. */
  val = info->naxis;
  fits_update_key (fptr, TINT, "NAXIS", &val, "", &status);
  if (status)
  {
    MFITSIO_PRINTF ("29 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not store NAXIS information.");
  }
  val = info->bitpix;
  fits_update_key (fptr, TINT, "BITPIX", &val, "", &status);
  if (status)
  {
    MFITSIO_PRINTF ("30 An error occurred: %d\n", status);
    MFITSIO_ERR ("Could not store BITPIX information.");
  }
  /** For each dimension, update the corresponding key value. */
  for (i = 0; i < info->naxis; i++)
  {
    snprintf (naxisf, 7, "NAXIS%d", i);
    fits_update_key (fptr, TINT, naxisf, (info->naxes + i), "", &status);
    if (status)
    {
      MFITSIO_PRINTF ("31 An error occurred: %d\n", status);
      MFITSIO_ERR ("Could not store NAXISX information.");
    }
  }
}

/**
 * Calculate bitpix and dimension information from a MATLAB image.
 *
 * @param img         The MATLAB array where the image is stored.
 *
 * @return            A mfitsio_info struct containing the bitpix
 *                    and dimension information.
 */

mfitsio_info *
mfitsio_calc_info (const mxArray * img)
{
  mfitsio_info *retval = 0;
  const int *naxesm;
  int i;
  mxClassID classId;
  /** Allocate memory for the info structure. */
  retval = (mfitsio_info *) MFITSIO_MALLOC (sizeof (mfitsio_info));

  /** Retrieve the axes information from the image. */
  retval->naxis = mxGetNumberOfDimensions (img);
  retval->naxes = (int *) MFITSIO_MALLOC (sizeof (int) * retval->naxis);
  retval->check = 0;

  /** Copy the dimension information into the info structure. */
  naxesm = mxGetDimensions (img);
  for (i = 0; i < retval->naxis; i++)
  {
    retval->naxes[i] = naxesm[i];
  }

  /** Depending on the class of the image, assign different bitpixes. */
  classId = mxGetClassID (img);
  switch (classId)
  {
  case mxSINGLE_CLASS:
    retval->bitpix = FLOAT_IMG;
    break;
  case mxDOUBLE_CLASS:
    retval->bitpix = DOUBLE_IMG;
    break;
  case mxUINT8_CLASS:
    retval->bitpix = BYTE_IMG;
    break;
  case mxINT8_CLASS:
    retval->bitpix = BYTE_IMG;
    break;
  case mxUINT16_CLASS:
    retval->bitpix = SHORT_IMG;
    break;
  case mxINT16_CLASS:
    retval->bitpix = SHORT_IMG;
    break;
  case mxUINT32_CLASS:
    retval->bitpix = LONG_IMG;
    break;
  case mxINT32_CLASS:
    retval->bitpix = LONG_IMG;
    break;
  case mxCHAR_CLASS:
  case mxUINT64_CLASS:
  case mxINT64_CLASS:
  case mxLOGICAL_CLASS:
  default:
    MFITSIO_PRINTF ("32 Unsupported data type for header image.\n");
    MFITSIO_ERR ("Cannot continue.");
    break;
  }
  return retval;
}

/**
 * Determine whether the MATLAB array is singleton.
 *
 * @param array       The MATLAB variable to check for singleton dimensions.
 *
 * @return            True is the array is a scalar.
 */

int
mfitsio_is_scalar (const mxArray * array)
{
  int ndim = 0, nsd = 1, i = 0;
  const int *dims;
  /** Get the total number of dimensions. */
  ndim = mxGetNumberOfDimensions (array);

  /** Get the dimensions array. */
  dims = mxGetDimensions (array);

  /** Check for a non-singleton dimension. */
  for (i = 0; i < ndim; i++)
  {

    /** If a non-singleton dimension is found, change the return status. */
    if (dims[i] != 1)
    {
      nsd = 0;
      break;
    }
  }
  return nsd;
}

/**
 * This function determines whether modifying a certain header keyword is
 * allowed. The reason for checking is that we do not want users modifying
 * NAXIS, BITPIX, and other important headers as they can be calculated from
 * a MATLAB image. Thus, the situation where the dimension and bitpix keywords
 * differ from the stored image is avoided.
 *
 * @param array       The MATLAB variable to check for singleton dimensions.
 *
 * @return            True is the array is a scalar.
 */

int
mfitsio_forbidden (const char *s)
{
  return !strncmp (s, "NAXIS", 5) || !strcmp (s, "SIMPLE")
    || !strcmp (s, "_nil_") || !strcmp (s, "BITPIX");
}

/**
 * This function determines whether modifying a certain header keyword is
 * allowed. The reason for checking is that we do not want users modifying
 * NAXIS, BITPIX, and other important headers as they can be calculated from
 * a MATLAB image. Thus, the situation where the dimension and bitpix keywords
 * differ from the stored image is avoided.
 *
 * @param array       The MATLAB variable to check for singleton dimensions.
 *
 * @return            True is the array is a scalar.
 */

int
mfitsio_ignore_card (const char *s)
{
  return !strncmp (s, "COMMENT  ", 9) || !strncmp (s, "LONGSTRN=", 9)
    || !strncmp (s, "EXTEND  =", 9) || !strncmp (s, "SIMPLE  =", 9);
}

/**
 * A wrapper function for mxGetField. Whenever a field cannot be retrieved,
 * an error is reported and execution aborts.
 *
 * @param array       The MATLAB struct array.
 * @param index       The index for the record.
 * @param field_name  The key for the record.
 *
 * @return            The value of the field requested.
 */

mxArray *
mfitsio_get_mfield (const mxArray * array, int index, const char *field_name)
{
  mxArray *field;
  if (array == 0)
  {
    MFITSIO_PRINTF ("34 Header does not exist.\n");
    MFITSIO_ERR ("Cannot continue.");
  }
  field = mxGetField (array, index, field_name);
  if (field == 0)
  {
    MFITSIO_PRINTF ("35 Field '%s' not found in MATLAB struct array.\n",
		     field_name);
    MFITSIO_ERR ("Cannot continue.");
  }
  return field;
}

/**
 * A wrapper function for mxGetScalar. Whenever a scalar cannot be retrieved,
 * zero is returned.
 *
 * @param array       The MATLAB struct array.
 * @param index       The index for the record.
 * @param field_name  The key for the record.
 *
 * @return            The numeric value as a double.
 */

double
mfitsio_get_mscalar (const mxArray * array)
{
  return (array) ? mxGetScalar (array) : 0;
}

/**
 * Checks for validity of starting and ending pixel coordinates. In particular,
 * pixels must be within the range defined by the NAXES keywords and a region
 * along a particular axis profile must not begin at the ending coordinate.
 *
 * @param      crd1   The starting coordinate.
 * @param      crd2   The ending coordinate.
 * @param      naxis  The number of axes.
 * @param      naxes  The size of each axis.
 */

void mfitsio_check_coordinate(const mxArray *crd1, const mxArray *crd2,
			      const int naxis, const int *naxes) {
  int i;
  double *data1 = 0, *data2 = 0;
  if (mxGetNumberOfDimensions(crd1) > 2) {
    MFITSIO_PRINTF ("40 Too many  in your starting pixel coordinate."
		    " It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetNumberOfDimensions(crd2) > 2) {
    MFITSIO_PRINTF ("41 Too many dimensions in your ending pixel coordinate. "
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetM(crd1) != 1) {
    MFITSIO_PRINTF ("42 Too many rows in your starting pixel coordinate. "
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetM(crd2) != 1) {
    MFITSIO_PRINTF ("43 Too many rows in your ending pixel coordinate. "
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetN(crd1)<naxis) {
    MFITSIO_PRINTF ("44 Not enough values in row vector defining"
		    " starting pixel coordinate. "
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetN(crd1)>naxis) {
    MFITSIO_PRINTF ("45 Too many values in row vector defining"
		    " starting pixel coordinate. "
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetN(crd2)<naxis) {
    MFITSIO_PRINTF ("46 Not enough values in row vector defining"
		    " end pixel coordinate."
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  if (mxGetN(crd2)>naxis) {
    MFITSIO_PRINTF ("47 Too many values in row vector defining"
		    " starting pixel coordinate."
		    "It must be a 1x%d vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  data1 = (double*) mxGetData(crd1);
  data2 = (double*) mxGetData(crd2);
  for (i = 0; i < naxis; i++ ) {
    if (data1[ i ] < 1 || data1[ i ] > naxes[ i ]) {
      MFITSIO_PRINTF("36 Invalid component at index %d of starting"
		     " pixel: %d\n", i + 1, (int)data1[ i ]);
      MFITSIO_ERR("Starting coordinate out of range.");
    }
    if (data2[ i ] < 1 || data2[ i ] > naxes[ i ]) {
      MFITSIO_PRINTF("37 Invalid component at index %d of ending"
		     " pixel: %d\n", i + 1, (int)data2[ i ]);
      MFITSIO_ERR("Ending coordinate out of range.");
    }
    if (data2[ i ] < data1[ i ]) {
      MFITSIO_PRINTF("38 Component at index %d of starting pixel must be"
                     "less than or equal to component %d of ending pixel.\n");
      MFITSIO_ERR("Invalid subset (hyper)region.");
    }
  }
}

/**
 * Checks for validity of size information. The NAXIS and NAXES information
 * must be non-negative.
 *
 * @param      size_unit   The size unit as a row vector.
 */

void mfitsio_check_size_info(const mxArray *size_unit) {

  int i;
  double *data = 0;
  int naxis;
  if (mxGetNumberOfDimensions(size_unit) > 2) {
    MFITSIO_PRINTF ("48 Too many dimensions defining your axes."
		    " The size unit must be a 1xN vector.\n");
    MFITSIO_ERR ("Invalid coordinate.");
  }

  if (mxGetM(size_unit) != 1) {
    MFITSIO_PRINTF ("49 Too many rows defining your axes. "
		    " The size unit must be a 1xN vector.\n", naxis);
    MFITSIO_ERR ("Invalid coordinate.");
  }
  data = (double*) mxGetData(size_unit);
  naxis = mxGetN(size_unit);
  for (i = 0; i < naxis; i++ ) {
    if (data[ i ] < 0) {
      MFITSIO_PRINTF("50 Invalid component at index %d of size unit"
		     ": %d\n", i + 1, (int)data[ i ]);
      MFITSIO_ERR("Size unit out of range.");
    }
  }
}


/**
 * Converts a double array to a long array.
 *
 * @param     dbl   The array to convert.
 * @param     size  The size of the array to convert.
 *
 * @return          The converted array.
 */

long *mfitsio_convert_dbl2long(const double *dbl, const int size) {
  int i;
  long *l = (long*)MFITSIO_MALLOC(sizeof(long)*size);
  for (i = 0; i < size; i++ ) {
    l[ i ] = (long)dbl[ i ];
  }
  return l;
}

/**
 * Converts an int array to a long array.
 *
 * @param     in    The array to convert.
 * @param     size  The size of the array to convert.
 *
 * @return          The converted array.
 */

long *mfitsio_convert_int2long(const int *in, const int size) {
  int i;
  long *l = (long*)MFITSIO_MALLOC(sizeof(long)*size);
  for (i = 0; i < size; i++ ) {
    l[ i ] = (long)in[ i ];
  }
  return l;
}

/**
 * Creates a one vector of a particular size.
 *
 * @param     size  The size of the ones array to create.
 *
 * @return          The converted array.
 */

long *mfitsio_create_ones_vector(const int size) {
  int i;
  long *l = (long*)MFITSIO_MALLOC(sizeof(long)*size);
  for (i = 0; i < size; i++ ) {
    l[ i ] = 1;
  }
  return l;
}

/**
 * Computes the size of a region defined by two coordinates.
 *
 * @param      crd1   The starting coordinate.
 * @param      crd2   The ending coordinate.
 * @param      naxis  The number of axes.
 *
 * @return            The size of the region as an array.
 */

int *mfitsio_get_region_size(const long *crd1, const long *crd2,
			     const int naxis) {
  int i;
  int *diff = 0;
  diff = (int*)MFITSIO_MALLOC(sizeof(int)*naxis);
  for (i = 0; i < naxis; i++ ) {
    diff[ i ] = (int)(crd2[ i ] - crd1[ i ]) + 1;
  }
  return diff;
}

/**
 * Computes the ending pixel value based on the starting pixel.
 *
 * @param      spixels The starting pixel.
 * @param      img     The image.
 * @param      naxis   The number of axes.
 *
 * @return             The size of the region as an array.
 */

mxArray *mfitsio_get_lpixels(const mxArray *spixels,
			     const mxArray *img, 
			     int naxis) {
  int i;
  int dim[] = { 1, 1 };
  mxArray *array;
  double *sdata;
  double *ldata;
  const int *imgdims;
  int nimgdims;
  dim[ 1 ] = naxis;
  imgdims = mxGetDimensions(img);
  array = mxCreateNumericArray (2, dim, mxDOUBLE_CLASS, mxREAL);  
  sdata = (double*) mxGetData(spixels);
  ldata = (double*) mxGetData(array);
  nimgdims = mxGetNumberOfDimensions(img);
  for (i = 0; i < nimgdims; i++ ) {
    ldata[ i ] = sdata[ i ] + (double)imgdims[ i ] - 1;
  }
  for (; i < naxis; i++ ) {
    ldata[ i ] = sdata[ i ];
  }
  return array;
}

char *mfitsio_logical_m2f(const mxArray *img) {
  const int *imgdims;
  mxLogical *logicals;
  char *retval = 0;
  long prod = 1, j;
  int nimgdims, i;
  logicals = mxGetLogicals(img);
  imgdims = mxGetDimensions(img);
  nimgdims = mxGetNumberOfDimensions(img);
  for (i = 0; i < nimgdims; i++ ) {
    prod *= ((long)imgdims[ i ]);
  }
  retval = (char*)MFITSIO_MALLOC(sizeof(char) * prod);
  for (j = 0; j < prod; j++ ) {
    retval[ j ] = logicals[ j ] ? 'T' : 'F';
  }
  return retval;
}
