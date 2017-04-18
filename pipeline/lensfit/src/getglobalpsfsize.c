#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int
getglobalpsfsize(char *psfname, int *correct_distortion, int *pwidth, int *pheight, int *order,
                 int *chipvariation, int *chiporder, int *ncoeffs,
                 int *nxchip, int *nychip, int *xchipsampling,
                 int *ychipsampling, int *big_gap, double *hxsize,
                 double *hysize, float *poserror)
{
    fitsfile *afptr;
    int status = 0;
    int anaxis;
    long anaxes[3] = { 1, 1, 1 };
    char comment[64];

    /*
     * open input images 
     */
    fits_open_file(&afptr, psfname, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    /*
     * read dimensions 
     */
    fits_get_img_dim(afptr, &anaxis, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis < 2 || anaxis > 3) {
        printf("Error: expected 2 or 3 dimensions in PSF fits file\n");
        exit(2);
    }

    fits_get_img_size(afptr, anaxis, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    if (anaxis == 2) {
        anaxes[2] = 1;
        order[0] = 0;
    } else {
        if (fits_read_key(afptr, TINT, "ORDER", order, comment, &status)) {
            fits_report_error(stderr, status);  /* print error message */
            printf(" keyword ORDER\n");
            return (status);
        }
    }

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
        (afptr, TINT, "CHIPVAR", chipvariation, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHIPVAR\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "CHORDER", chiporder, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword CHORDER\n");
        return (status);
    }
    if (fits_read_key(afptr, TINT, "NXCHIP", nxchip, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword NXCHIP\n");
        *nxchip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "NYCHIP", nychip, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword NYCHIP\n");
        *nychip = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "XSAMPL", xchipsampling, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword XSAMPL\n");
        *xchipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "YSAMPL", ychipsampling, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword YSAMPL\n");
        *ychipsampling = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TINT, "BIG_GAP", big_gap, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword BIG_GAP\n");
        *big_gap = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HXSIZE", hxsize, comment, &status)) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" keyword HXSIZE\n");
        *hxsize = 0;
        status = 0;
    }
    if (fits_read_key(afptr, TDOUBLE, "HYSIZE", hysize, comment, &status)) {
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

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        return (status);
    }

    pwidth[0] = (int) anaxes[0];
    pheight[0] = (int) anaxes[1];
    ncoeffs[0] = (int) anaxes[2];

    return (0);

}
