#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

/* this version also gets keywords for shifts fits */

int
getglobalpsfsizes(char *psfname, int *correct_distortion, int *pwidth, int *pheight, int *order,
                  int *chipvariation, int *chiporder, int *ncoeffs,
                  int *numstar, int *nxchip, int *nychip, int *xchipsampling,
                  int *ychipsampling, int *big_gap, double *hxsize,
                  double *hysize, float *poserror, double *scalefactor, int *sorder,
                  int *schipvariation, int *schiporder, int *sncoeffs)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis, hdunum, hdutype;
    long anaxes[3] = { 1, 1, 1 };
    char comment[64];
    long nrows;
    double dnull;
    float fnull;
    int anynull, ncol;
    int noglobalshift = { 0 };
    long frow, felem, nelem;

    /*
     * open input images 
     */
    fits_open_file(&afptr, psfname, READONLY, &status);
    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * read dimensions 
     */
    fits_get_img_dim(afptr, &anaxis, &status);
    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis < 2 || anaxis > 3) {
        fflush(stdout);
        fprintf(stderr,
                "Error: expected 2 or 3 dimensions in PSF fits file\n");
        exit(2);
    }

    fits_get_img_size(afptr, anaxis, anaxes, &status);
    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis == 2) {
        anaxes[2] = 1;
        order[0] = 0;
    } else {
        if (fits_read_key(afptr, TINT, "ORDER", order, comment, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword ORDER\n");
            return (status);
        }
    }

    pwidth[0] = (int) anaxes[0];
    pheight[0] = (int) anaxes[1];
    ncoeffs[0] = (int) anaxes[2];

    /*
     * get remaining keywords 
     */
    if (fits_read_key(afptr, TINT, "CORRECTD", correct_distortion, comment, &status)) {
      fflush(stdout);
      fits_report_error(stderr, status);  /* print error message */
      fprintf(stderr," keyword CORRECTD\n");
      fprintf(stderr," distortion correction keyword\n");
      return (status);
    }
    if (fits_read_key
        (afptr, TDOUBLE, "PIXSCALE", scalefactor, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        fprintf(stderr, " keyword PIXSCALE\n");
        exit(2);
    }
    if (fits_read_key
        (afptr, TINT, "CHIPVAR", chipvariation, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHIPVAR\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "CHORDER", chiporder, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHORDER\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "NXCHIP", nxchip, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword NXCHIP\n");
        *nxchip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NYCHIP", nychip, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword NYCHIP\n");
        *nychip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "XSAMPL", xchipsampling, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword XSAMPL\n");
        *xchipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "YSAMPL", ychipsampling, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword YSAMPL\n");
        *ychipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "BIG_GAP", big_gap, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword BIG_GAP\n");
        *big_gap = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HXSIZE", hxsize, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword HXSIZE\n");
        *hxsize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HYSIZE", hysize, comment, &status)) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword HYSIZE\n");
        *hysize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TFLOAT, "TOL", poserror, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword TOL\n");
        *poserror = 0;
        status = 0;
    }

    if (fits_read_key(afptr, TINT, "SORDER", sorder, comment, &status)) {
        fflush(stdout);
        printf
            (" warning from getglobalpsfsizes: globalshifts keyword SORDER has not been found \n");
        //fits_report_error(stderr, status); /* print error message */
        //printf(" keyword SORDER\n");
        noglobalshift = 1;
    }

    if (noglobalshift == 0) {
        if (fits_read_key
            (afptr, TINT, "SCHIPVAR", schipvariation, comment, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword SCHIPVAR\n");
            exit(2);
        }
        if (fits_read_key
            (afptr, TINT, "SCHORDER", schiporder, comment, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword SCHORDER\n");
            exit(2);
        }
    }

    /*
     * read table of star numbers in each chip 
     */
    /*
     * get the table of star info 
     */
    hdunum = 2;
    if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
        status = 0;
        fflush(stdout);
        fprintf(stderr, " error moving to table of star information \n");
        fits_report_error(stderr, status);      /* print error message */
        status = 0;
	return(0);
    } else {
        if (hdutype != BINARY_TBL) {
            fflush(stdout);
            fprintf(stderr,
                    " Error file format not as expected - Table of star info is not a binary table \n");
            exit(2);
        }

        if (fits_get_num_cols(afptr, &ncol, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            return (status);
        }
        if (ncol != 2) {
            fflush(stdout);
            fprintf(stderr, " error in table format, no. columns = %d \n",
                    ncol);
            exit(2);
        }

        if (fits_get_num_rows(afptr, &nrows, &status)) {
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }
        if (nrows != nxchip[0] * nychip[0]) {
            fflush(stdout);
            fprintf(stderr,
                    " error, number of table entries does not match array size \n");
            exit(2);
        }

        frow = 1;
        felem = 1;
        nelem = nxchip[0] * nychip[0];
        dnull = 0.;
        fnull = 0.;
        anynull = 0;
        if (fits_read_col
            (afptr, TINT, 2, frow, felem, nelem, &dnull, numstar, &anynull,
             &status)) {
            fflush(stdout);
            fprintf(stderr, " error reading dec column in table \n");
            fflush(stdout);
            fits_report_error(stderr, status);  /* print error message */
            exit(2);
        }
    }

    /*
     * repeat for x,y shifts coefficients 
     */
    /*
     * read table of star numbers in each chip 
     */
    /*
     * get the table of star info 
     */

    if (noglobalshift == 0) {
        hdunum = 3;
        if (fits_movabs_hdu(afptr, hdunum, &hdutype, &status)) {
            status = 0;
            fflush(stdout);
            fprintf(stderr,
                    " error moving to table of PSF shifts information \n");
            fits_report_error(stderr, status);  /* print error message */
            fprintf(stderr, " run program globalshifts first \n");
            exit(2);
        } else {
            if (hdutype != BINARY_TBL) {
                fflush(stdout);
                fprintf(stderr,
                        " Error file format not as expected - Table of star info is not a binary table \n");
                exit(2);
            }
            if (fits_get_num_cols(afptr, &ncol, &status)) {
                fflush(stdout);
                fits_report_error(stderr, status);      /* print error message */
                return (status);
            }
            if (ncol != 2) {
                fflush(stdout);
                fprintf(stderr, " error in table format, no. columns = %d \n",
                        ncol);
                exit(2);
            }

            if (fits_get_num_rows(afptr, &nrows, &status)) {
                fflush(stdout);
                fits_report_error(stderr, status);      /* print error message */
                exit(2);
            }

            *sncoeffs = nrows;
        }
    }

    fits_close_file(afptr, &status);

    if (status) {
        fflush(stdout);
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }



    return (0);

}
