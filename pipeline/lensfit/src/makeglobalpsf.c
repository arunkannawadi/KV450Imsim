/*--------------------------------------------------------------------------

makeglobalpsf v7

Code to create a surface fit to the PSF, based on simpler stacking
algorithm, makepsf.  Here's the preamble for makepsf:

[Code to create PSFs suitable for lensfit, by stacking stars in an
image.  Requires a list of input star positions plus the matching
image.  Stars with max brightness more than 50 percent of the CCD
saturation level are not used.  Magnitude ranges and peak S/N may also be
specified.  A list of galaxy positions is also used to eliminate
neighbouring galaxies.  Stars are sub-pixel registered onto the
stacked PSF in a two-step interative process.  Individual stars are
cross-correlated with the stacked PSF and any that are a poor match
are thrown out.  The PSF is assumed to be position dependent, and is
defined in a series of boxes dividing up the input image.  Stars are
accumulated into the PSF for each box if they fall within some
distance of that box.  The sampling box size and smoothing are both
specified on the command line.  Output is a three-dimensional FITS
image, the first two dimensions are x and y position on the image,
the next dimension is used for each 2D PSF stacked image.
There is also a table of the x and y positions of the stars that have
been used for each PSF.]

The last part is replaced by a polynomial fit, which in this code is
a fit across all chips in an exposure, with optional chip-dependent
variation

History:

original makepsfsurface first created Lance Miller. 20/01/2007
modified to crosscorrelate initially with a delta function to optimise centering
LM 27/02/2007

modified to write out single FITS file containing all PSFs, plus smoothing
algorithm and additional options LM, 17/1/2008

command-line input modified to force S/N ratio to be compulsory, also
change to output information and time estimate included LM, 30/1/2008

addition of background value into extractdata requires change to makepsf
also, needed to get correct saturation value for bright stars LM, 30/1/2008

modification to getimagesize to return CCD PA (not used in makepsf
but needed for function compatibility) LM, 4/3/2008

Change to fit polynomial surface to pixel values LM, 7/5/2008
Also read bad pixel mask LM 12/5/2008
Improve cross-correlation method  LM 14/7/2008
Made compatible with new extractdata routine LM 23/7/08
Moved to extractpostagestamp routine LM 23/9/08
Added in option of WCS coords LM 26/9/08
Allowed option of selecting additional stars to be used LM sep 08
Modified memory allocation for varylpthresholdf arrays LM Dec 08

Modified to allow swarp correction of astrometric distortion, for
consistency with lensfit v5, LM 16/01/2009

Now writes info on numbers of stars used in each psf box into the 
coefficients file, so that"bad" regions may be flagged out in use.
LM 21/1/2009

This version make global polynomial fit to the ensemble of chips
in an exposure, with optional chip-dependence of low-order 
coefficients
LM Jan 2009

Corrected inconsistencies between swarp and non-swarp versions
LM 4/3/09

Various minor changes for consistency with lensfit modifications, LM  Apr-Oct 09

Corrected bug in output RA,dec values (stars were being muddled-up) LM 4 Nov 2009

ported to SVN system LM 5th August 2010

aligned with May 7 verion, LM 7th August 2010

corrected the delimiter bug and aligned with v7, LM 17 Nov 2010

Revision $Rev: 167 $ 
last changed by $Author: miller $
on $LastChangedDate: 2013-04-27 15:42:06 +0200 (Sat, 27 Apr 2013) $

---------------------------------------------------------------------------*/

#define DATE_STRING "$LastChangedDate: 2013-04-27 15:42:06 +0200 (Sat, 27 Apr 2013) $"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

/* swarp wcs include files (but note custom version of fitscat.h to avoid conflict with cfitsio */
#include "define.h"
#include "types.h"
#include "globals.h"
#include "./fitscat.h"
#include "fitswcs.h"
#include "data.h"
#include "field.h"
#include "header.h"
#include "interpolate.h"
#include "prefs.h"
#include "projapprox.h"
#include "resample.h"
#include "weight.h"
#include "wcs/wcs.h"
#include "prefs.h"

#include "lensutil.h"

#include <complex.h>
#include <fftw3.h>

int extractpostagestamp(float *, float *, int *, double *, double *,
                        float *, float *, float, float, float,
                        float, float, int, int, float, int *, float);
void svdfit2dsqc(double *, double *, double *, double *, int, int, int,
                 double *, double **, double **, double *);
int tbswarpextract(wcsstruct *, float *, float *, int *, double *, double *,
                   int *, float, float, float, int, int, ikernelstruct *,
                   float, int *, float, double *, double *, double *,
                   double *, double *, double **, double *, double *,
                   double *, double *, double **, double **);
int tbswarpextract_nodistortion(wcsstruct *, float *, float *, float *, int *,
                                double *, double *, float, float,
                                float, int, int, float,
                                int *, double *, float, double *, double *,
                                double *, double *, int);
int writeglobalcoeffs(char *, int, double **, int, int, int, int, double *, int,
                      int, int, int, int *, int, int, int, double, double, float);
int writepsf_fits(char *, double **, int, int, int, int, int *, float **,
                  float **, int, int);
int readpsfcatsize(char *);     // read size of catalogue
/* choose appropriate routine for catalogue format */
//int readpsfcat (char*, float*, float*, float*);
//int readpsfcat2_wcs (char*, double*, double*, float*);
int readpsfcat_simple_wcs(char *, double *, double *, float *);
int getimagesize(char *, char *, int *, float *, float *, float *, float *,
                 int *, int *, float);
void getdata_badpix(char *, char *, int *, float *, float *, float *);
int weightfilter(float *, float *, float *, float *, int *, int);
void dilatemask(float *, int *);
int xcorr_measure(fftw_plan, fftw_complex *, fftw_complex *,
                  int, int, int, int,
                  fftw_complex *, double *, fftw_complex *, double *);
int shiftft(int, int, fftw_complex *, double *);
void globalsvdfit(double *, double *, double *, double *, int *, int, int,
                  int, int, int, int, int, double *, double **, double **,
                  double *);
void globalreconstruct(double, double, int, int, int, int, int, int, double *,
                       double *);
int varylpthresholdf2(float *, int *, float, int, int, int, int, int, float,
                      float, int);
void mediansub(float *, float *, int *, float, float *);
int wcs_to_raw(wcsstruct *, double *, double *);
void load_dfield(fieldstruct *, catstruct *, char *, int, int);
ikernelstruct *init_ikernel(interpenum *, int);
char *create_psf_name(const char *, char *, const char *);
int write_psfsurface(char *, char *, double **, int, int, int, int, int *,
                     float **, float **, int, int);
void writestars(int, int, int, double **, double **, int *, double *,
                double *, float *, float *, float *, const char *,
                const char *, char *);
void get_poly_order(int, int, int, int, int *, int *);
int check_psf_norm(double **, int, int, int, int, int);

static int
compare(const void *ii, const void *jj)
{
    float *i, *j;
    i = (float *) ii;
    j = (float *) jj;
    if (*j > *i)
        return -1;
    if (*j < *i)
        return 1;
    return 0;
}

/* verbose flag */
#define VERBOSE 1

/* set whether or not to subtract a constant median background 
   (0=no 1=yes) */
#define SUBTRACT_MEDIAN 0

/* set whether or not to apply the median erosion filter to the weight image
   (same routine as weightimagefilter, used on CFHTLS weight images) */
#define FILTER_WEIGHTS 1

/* specify that weight (bad pixel) images are required */
#define WEIGHTS_REQUIRED 0

/* set number of stars to be used for averaging within each region */
#define NAVERAGE 20

/* set minimum number of stars required for subregions to be considered valid 
   (only used for warning) (regions with numbers equal to this are considered acceptable) 
#define PSF_STAR_LIMIT 1
*/

/* define whether input galaxy catalogue positions are in WCS (WCS=1)
   or pixel (WCS=0) coordinates */
#define WCS 1

/* define whether or not to use swarp astrometric corrections
   requires WCS 1 */
#define USE_SWARP 1

/* set whether or not to apply the internal swarp distortion correction (1=yes) 
   requires WCS=1 and USE_SWARP=1 also */
#define CORRECT_DISTORTION 0

/* set whether to output 3D fits images of the stars used (shifted and normalised) */
#define OUTPUT_STARS 1          // 1=yes

/* define whether input files are gzipped fits or plain fits
   0=not zipped, 1=weight files are zipped, 2=both data and weight files zipped */
#define GZIPPED 1

/* define whether to output individual psf surface files */
#define WRITE_INDIVIDUAL_SURFACE 0

/* define minimum accepted exposure time (must be the same as in lensfit) */
#define minexposuretime 0.

/* Set the instrument here by defining the flag that is checked later
   in the code KIDS, CFHT, SUPRIME or SINGLE */
#define SUPRIME (1)
/*---------------------------------------------------------------------------*/

/* global variables */
float poserror;

/* main program */
int
main(int argc, char *argv[])
{
    int halfpwidth, halfpheight, dim[2], badccd, imageid;
    int padwidth, padheight, halfpadwidth, cchwidth;
    int chipvariation, chiporder, chipnumber;
    int dist1;
    int niter;
    int flagbad;
    int i, j, k, ii;
    int x, y, x1, x2, y1, y2;
    int xmin, xmax, ymin, ymax;
    int iterate, warning;
    int pheight, pwidth, imagesize;
    int len, nobj, nobjt, nobj2;
    int boxsize, boxsmooth = 0, maxboxno, xboxnum, yboxnum, fitsize,
        halfboxsmooth;
    int psfnum, psfsize;
    int *num, *numstar;
    int **boxno, *box, *mainboxno;
    int *bad;
    int *listnum, *badnum;      /* numbers of good and bad stars in each psf */
    int nstars, nstars1, totalbad, totalgood, tot_nstar;
    int ix, iy, pixel, xx, yy;
    int num_off, num_wrong_mag, num_saturated, num_close, num_snratio;
    int num_badpix;
    int d, close, good;
    int ngoodpix, imax_intensity;
    int *nn;
    int *region, *objpix, nobjects;
    int lower_area_limit, lower_merging_limit, arealimit;
    int *offscale;
    int naxis = 2;
    int narg = 0;
    int fheight = 20;  // half-height of median filter in weightfilter
    int nchanged;
    int kimage, refimage, refchip;
    int image, nchip, xchipsampling, ychipsampling, big_gap, nchip_expected,
        nxchip, nychip;
    int **chip, xchip, ychip, ichip;

    time_t t1, t2, t3, t4, t5, timetest0;
    float memory, memloop;
    int timereported;

    double *shift, *xshift, *yshift, sumw;
    double *psfauto, dauto, cross, **norm;
    double **count, **newpsf;
    double *xfit, *yfit, *zfit, *wfit, *weight;
    double **xfits, **yfits, **zfits, **wfits;
    double xval, yval, rsq, rsq0, sumb;
    double xsize, ysize, hxsize, hysize;
    double *radegs, *decdegs;   // ra,dec in input catalogue
    double *rdegs, *ddegs;      // ra,dec arrays after stars off-image have been eliminated  
    double scalefactor[1];
    double cosmicraylimit;

    float maxlevel;
    float noise;
    float *opix, *changed, *weightfilterarray;
    float *apix, *temp, *badpix, *badtemp;
    float *sortarray, *sn;
    float fintensity_limit, fmax_intensity, scale;
    float *objx, *objy, *mag, *rmag;
    float brightmag, faintmag, snratio;
    float satlev, gain, satlimit, arcperpix, angle;
    float ximageoffset, yimageoffset;
    float **listx, **listy;
    double starweight;
    double **u, **v, *w, *avals;
    double *dobjx, *dobjy;
    int order, crossterm, ncoeffs, maxncoeffs;
    int *goodstar, nstar;

    double **psf, **npsfpixel, **acoeffs;
    double *c, *shiftpix;
    double *bpix, **dpix, **dbadpix;
    double maxvalue, distantmaxvalue, *normvalue, psfnorm, noisevalue;
    double rawpos[2], wcspos[2], wcsneg[2], wcscentre[2], *kern[2];

    fftw_complex *Bpix, *C, **Dpix, *Cpad, *shiftpixft;
    fftw_plan qplan, iqplan, padinv;

    char *evariable, *headevariable, *badpixdir, *psfdir;
    char *prefsname;
    char *catname;
    char **image_file;
    char *imagename, *headerfile;
    char *psfname, *coeffname, *fittedpsfname, *imageweightname;
    char *firstname, *badname, *shiftsname, *usedstarsname;
    char **argkey, **argval;
    char *pstr, delimiter[5], headdelimiter[5], *satname, *weight_suffix;

    FILE *shiftsfile, *usedstarsfile, *filep;

#ifdef KIDS
    char delims[] = "_";
    char delims2[] = "O";
#endif
#ifdef CFHT
    char delims[] = "_C";
    char delims2[] = "_C";
#endif
#ifdef SUPRIME
    char delims[] = "_O";
    char delims2[] = "_O";
#endif
#ifdef SINGLE
    char delims[] = "_C";
    char delims2[] = "_C";
#endif
    char *item = NULL;

    /*
     * variables for output fits file of stars used 
     */
    int status = 0;
    /*
     * table to identify each object 
     */
    char dotdelimiter[] = { "." };

    /*
     * swarp wcs variables 
     */
    catstruct *rawcat;
    rawcat = (catstruct *) NULL;
    fieldstruct *rawfield;
    rawfield = (fieldstruct *) NULL;
    wcsstruct *wcs_raw;
    wcs_raw = (wcsstruct *) NULL;
    ikernelstruct *ikernel;

    t1 = time(NULL);

    printf("\n %s version 7.2 (KiDS) ",argv[0]);
    char datedelims[]="()";
    char *dstring = NULL;
    char date_string[200];
    bzero(date_string,200);
    strcpy(date_string, DATE_STRING);
    dstring  = strtok(date_string, datedelims);
    if (dstring != NULL) dstring = strtok(NULL, datedelims);
    if (dstring != NULL) printf(" %s ",dstring);
    printf("\n");

    if (argc != 5 && argc != 7) {
        printf
            (" %s <imagelist> <global order> <chip order> <snratio> [<bright mag> <faint mag>] \n",
             argv[0]);
        exit(EXIT_FAILURE);
    }

    if (USE_SWARP == 1) {
        printf(" swarp astrometric corrections will be applied \n");
        if (WCS == 0) {
            printf
                (" WARNING:  WCS flag not set, using pixel object coordinates \n");
            printf
                (" even though WCS information is needed for the swarp correction \n");
        }

    if (CORRECT_DISTORTION == 1) {
      printf(" distortion corrections will be applied to the PSF\n");
        if (USE_SWARP == 0 || WCS == 0) {
            fprintf(stderr,
                    " CORRECT_DISTORTION is set but USE_SWARP or WCS are not\n");
            exit(EXIT_FAILURE);
        }
    } else {
      printf(" distortion corrections will not be applied to the PSF\n");
    }

        /*
         * read the swarp preferences file 
         */
        narg = 0;
        naxis = 2;
        QMALLOC(argkey, char *, 1);
        QMALLOC(argval, char *, 1);
        prefsname = getenv("SWARP_CONFIG");
        if (prefsname != NULL) {
            strcpy(prefs.prefs_name, prefsname);
            readprefs(prefs.prefs_name, argkey, argval, narg);
            prefs.verbose_type = QUIET; // turn off verbose reporting
            useprefs();
        } else {
            fprintf(stderr,
                    "configuration file environment variable SWARP_CONFIG not set \n");
            exit(EXIT_FAILURE);
        }
        free(argkey);
        free(argval);
        /*
         * initialise the convolution kernel 
         */
        ikernel = init_ikernel(prefs.resamp_type, naxis);
        for(i = 0; i < 2; i++)
            kern[i] = calloc(INTERP_MAXKERNELWIDTH, sizeof(double));
    } else {
        printf(" no swarp astrometric corrections will be applied \n");
    }
    /*
     * initialise the scalefactor 
     */
    *scalefactor = 0.;

    /*
     * test whether magnitude limits for PSF stars have been specified and read values
     * if so, or else set default values.  NB PSF star magnitude values need to be
     * given in the input PSF catalogue in the first case! 
     */

    if (argc == 7) {
        brightmag = atof(argv[5]);
        faintmag = atof(argv[6]);
        printf(" selecting PSF stars with %8.2f < m < %8.2f \n", brightmag,
               faintmag);
    } else {
        brightmag = -1.e10;
        faintmag = 1.e10;
    }

    order = atoi(argv[2]);

    if (strcmp(argv[3], "none") == 0) {
        chipvariation = 0;
        chiporder = 0;
        printf(" no chip variation allowed \n");
    } else {
        chipvariation = 1;
        chiporder = atoi(argv[3]);
        printf(" chip variation order %d \n", chiporder);
    }

    snratio = atof(argv[4]);
    if (snratio <= 0.)
        snratio = 0.;
    printf(" selecting stars with peak S/N ratio > %8.1f \n", snratio);
    if (snratio < 20.)
        printf(" WARNING: low S/N ratio not advised, choose S/N > 20 \n");

    /*
    if (WCS == 1) {
        printf(" reading WCS coordinates from input catalogue\n");
    } else {
        printf(" reading xy coordinates from input catalogue\n");
    }
    */

    /*
     * set some dimensions 
     */

    pwidth = 32;                /*size of each individual postage stamp (data and psf) */
    pheight = pwidth;           /* postage stamps must be square */
    halfpwidth = pwidth / 2 + 1;
    halfpheight = pheight / 2;
    psfsize = pheight * halfpwidth;
    padheight = 50 * pheight;
    padwidth = 50 * pwidth;
    halfpadwidth = 1 + padwidth / 2;
    cchwidth = pwidth / 4;

    if (padwidth <= pwidth || padheight <= pheight) {
        printf
            (" padded cross-correlation array size must be bigger than input size \n");
        exit(EXIT_FAILURE);
    }

    /*
     * initialise memory counters 
     */
    memory = 0.;
    memloop = 0;

    xfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "xfit");
    yfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "yfit");
    zfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "zfit");
    wfit = ml_calloc(pwidth * pheight, sizeof(double), &memory, "wfit");
    objpix = ml_calloc(pwidth * pheight, sizeof(double), &memory, "objpix");

    /*
     * set tolerance - 
     * (i) if any other object lies within this radius
     * the object will be rejected inside extractpostagestamp/swarpextract.
     * (ii) if any bad pixels lie within this radius the star will be rejected
     * (iii) star shifts larger than this amount will cause stars to be
     * rejected
     */

    // tolerance radius in pixels
    poserror = 6;

    evariable = ml_calloc(PATH_MAX, sizeof(char), &memory, "evariable");
    headevariable = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "headevariable");
    badpixdir = ml_calloc(PATH_MAX, sizeof(char), &memory, "badpixdir");
    psfdir = ml_calloc(PATH_MAX, sizeof(char), &memory, "psfdir");
    prefsname = ml_calloc(PATH_MAX, sizeof(char), &memory, "prefsname");
    catname = ml_calloc(PATH_MAX, sizeof(char), &memory, "catname");
    imagename = ml_calloc(PATH_MAX, sizeof(char), &memory, "imagename");
    headerfile = ml_calloc(PATH_MAX, sizeof(char), &memory, "headerfile");
    psfname = ml_calloc(PATH_MAX, sizeof(char), &memory, "psfname");
    fittedpsfname = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "fittedpsfname");
    usedstarsname = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "usedstarsname");
    firstname = ml_calloc(PATH_MAX, sizeof(char), &memory, "firstname");
    imageweightname = ml_calloc(PATH_MAX, sizeof(char), &memory,
                                "imageweightname");
    badname = ml_calloc(PATH_MAX, sizeof(char), &memory, "badname");
    shiftsname = ml_calloc(PATH_MAX, sizeof(char), &memory, "shiftsname");
    weight_suffix = ml_calloc(PATH_MAX, sizeof(char), &memory,
                              "weigh_suffix");
    satname = ml_calloc(PATH_MAX, sizeof(char), &memory, "satname");

    /*
     * create output PSF filename 
     */

    psfdir = getenv("PSF_DIR");
    strcpy(firstname, argv[1]);
    coeffname = create_psf_name(psfdir, firstname, dotdelimiter);

        /*
         * define the string delimiter used to allow the code to
         * create the filename for the .head file and weight file, 
	 * in the case where
         * USE_SWARP is defined to be 1.  The code will take the last
         * instance of this delimiter in the input filename string,
         * strip off the delimiter and following characters, and add
         * on the head file suffix defined in the swarp configuration
         * file (normally ".head").  e.g.  715231p_18C.sub.fits will
         * have head filename 715231p_18.head if delimiter is set to
         * "C".  If delimiter is set to "." the expected head file
         * would be 715231p_18C.sub.head The code will exit if this
         * head file is not found in the same path as the data file.
         */

    if (USE_SWARP == 1)
      {
#ifdef KIDS
        strcpy(weight_suffix, "F.weight.fits");
	strcpy(delimiter,"F");
	strcpy(headdelimiter,"O");
#endif
#ifdef CFHT
        strcpy(weight_suffix, "C.weight.fits");
	strcpy(delimiter,"C");
	strcpy(headdelimiter,"C");
#endif
#ifdef SINGLE
        strcpy(weight_suffix, "C.weight.fits");
	strcpy(delimiter,"C");
	strcpy(headdelimiter,"C");
#endif
#ifdef SUPRIME
	strcpy(weight_suffix, "OFCS.weight.fits");
	strcpy(delimiter,"O");
	strcpy(headdelimiter,"O");
#endif
      }
    else
      {
	strcpy(weight_suffix, ".weight.fits");
	strcpy(delimiter,".");
      }

    /*
     * set up the chip mosaic geometry
     */
#ifdef KIDS
    xchipsampling = 2040 + 70;
    ychipsampling = 4090 + 70;
    big_gap = 0;
    nxchip = 8;
    nychip = 4;
#endif
#ifdef CFHT
    xchipsampling = 2048 + 70;
    ychipsampling = 4612 + 70;
    big_gap = 425 - 70;
    nxchip = 9;
    nychip = 4;
#endif
#ifdef SUPRIME
    xchipsampling = 2013 + 70;  // total x size covered by chip plus gap
    ychipsampling = 4085 + 70;  // ditto for y
    big_gap = 0;                // allowance for extra gaps between some rows
    nxchip = 5;                 // number of chips on x axis
    nychip = 2;                 // number of chips on y axis
#endif
#ifdef SINGLE
  /* single chip geometry */
  xchipsampling = 4096;  // set this to the x size of the single image
  ychipsampling = 4096;  // same for y size
  big_gap = 0; // no gap between chips (because there's only one)
  nxchip = 1;  // one chip on x-axis
  nychip = 1;  // one chip on y-axis
#endif

    nchip_expected = nxchip * nychip;   // total number of chips
    hxsize = (double) (nxchip * xchipsampling) / 2.;    // half the total pixels on x axis
    hysize = (double) (nychip * ychipsampling + 2 * big_gap) / 2.;      // half the total pixels on y axis

    image_file = calloc(nchip_expected, sizeof(char *));

    /*
     * read the input file list 
     */
    filep = fopen(argv[1], "r");
    if (filep == NULL) {
        fprintf(stderr, " Error opening file %s \n", argv[1]);
        exit(EXIT_FAILURE);
    }

    nchip = 0;
    while (!feof(filep) && !ferror(filep)) {
        if (fscanf(filep, "%s", firstname) != 0) {
            if (!feof(filep) && !ferror(filep)) {
                if (nchip >= nchip_expected) {
                    fprintf(stderr, " too many images in file \n");
                    exit(EXIT_FAILURE);
                }
                image_file[nchip] = calloc(PATH_MAX, sizeof(char));
                strcpy(image_file[nchip], firstname);
                nchip++;
            }
        }
    }
    fclose(filep);
    status = 0;
    bzero(firstname, PATH_MAX);
    printf(" filenames read for %d input images \n", nchip);

    if (nchip<=0)
      {
	fflush(stdout);
	fprintf(stderr," no input images read from input list \n");
	exit(2);
      }

    if (nchip<nchip_expected)
      {
	printf(" WARNING: %d chips read from list but %d were expected \n",nchip,nchip_expected);
	fflush(stdout);
      }

    /*
     * get the order of the polynomial fit.  negative order is taken to mean
     * that cross-terms should be included 
     */
    get_poly_order(order, chiporder, chipvariation, nchip, &ncoeffs,
                   &crossterm);

    maxncoeffs = ncoeffs;
    if (maxncoeffs < 3)
        maxncoeffs = 3;


    u = ml_calloc((1 + maxncoeffs), sizeof(double *), &memory, "u");
    for(i = 0; i <= maxncoeffs; i++)
        u[i] = ml_calloc((1 + maxncoeffs), sizeof(double), &memory, "u[i]");

    v = (double **) calloc((1 + maxncoeffs), sizeof(double *));
    memory += (float) (1 + maxncoeffs) * sizeof(double);
    for(i = 0; i <= maxncoeffs; i++) {
        memory += (float) (1 + maxncoeffs) * sizeof(double);
        v[i] = (double *) calloc((1 + maxncoeffs), sizeof(double));
    }

    avals = (double *) calloc(1 + maxncoeffs, sizeof(double));
    w = (double *) calloc(1 + maxncoeffs, sizeof(double));
    memory += 2. * (1 + maxncoeffs) * sizeof(double);

    acoeffs = (double **) calloc((pwidth * pheight), sizeof(double *));
    memory += (float) (pwidth * pheight) * sizeof(double *);
    if (acoeffs == NULL) {
        printf("Memory allocation error for sub-image pointers\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < (pwidth * pheight); i++) {
        acoeffs[i] = (double *) calloc((ncoeffs), sizeof(double));
        memory += (float) (ncoeffs) * sizeof(double);
        if (acoeffs[i] == NULL) {
            printf("Memory allocation error for sub-image \n");
            exit(EXIT_FAILURE);
        }
    }


    /*
     * alocate memory for FFTW arrays etc 
     */

    Bpix = ml_calloc(halfpwidth * pheight, sizeof(fftw_complex), &memory,
                     "Bpix");
    shiftpixft =
        ml_calloc(halfpwidth * pheight, sizeof(fftw_complex), &memory,
                  "shiftpixt");
    bpix = ml_calloc(pwidth * pheight, sizeof(fftw_complex), &memory, "bpix");
    shiftpix = ml_calloc(pwidth * pheight, sizeof(fftw_complex), &memory,
                         "shiftpix");
    C = ml_calloc(halfpwidth * pheight, sizeof(fftw_complex), &memory, "C");
    c = ml_calloc(padwidth * padheight, sizeof(double), &memory, "c");
    Cpad = ml_calloc(halfpadwidth * padheight, sizeof(fftw_complex), &memory,
                     "Cpad");
    /*
     * temp array is used as spare storage inside extractdata 
     * for quadrant swapping when the postage stamps are extracted 
     */
    temp = ml_calloc(pwidth * pheight, sizeof(float), &memory, "temp");
    badtemp = ml_calloc(pwidth * pheight, sizeof(float), &memory, "badtemp");
    shift = ml_calloc(2, sizeof(double), &memory, "shift");

    /*
     * create FFTW plans 
     */
    qplan = fftw_plan_dft_r2c_2d(pwidth, pheight, bpix, Bpix, FFTW_MEASURE);
    iqplan =
        fftw_plan_dft_c2r_2d(pwidth, pheight, shiftpixft, shiftpix,
                             FFTW_MEASURE);
    padinv = fftw_plan_dft_c2r_2d(padwidth, padheight, Cpad, c, FFTW_MEASURE);

    /*
     * get environment variables holding directory names 
     */

    evariable = getenv("DATA_DIR");
    headevariable = getenv("HEAD_DIR");
    badpixdir = getenv("BADPIX_DIR");

    /*
     * create the psfstars catalogue filename 
     */
    catname = getenv("CATALOGUE_STARS");
    if (catname == NULL) {
        fprintf(stderr,
                " catalogue name not set, use environment variable CATALOGUE_STARS \n");
        exit(EXIT_FAILURE);
    }

    if (access(catname, F_OK) != 0) {
        printf(" Can't read catalogue %s \n", catname);
        exit(EXIT_FAILURE);
    }

    /*
     * read catalogue of objects 
     */

    nobjt = readpsfcatsize(catname);

    objx = ml_calloc(nobjt, sizeof(float), &memory, "objx");
    objy = ml_calloc(nobjt, sizeof(float), &memory, "objy");
    mag = ml_calloc(nobjt, sizeof(float), &memory, "mag");

    if (WCS == 0) {
        fprintf(stderr,
                " x,y star coordinate input not supported in this version \n");
        exit(EXIT_FAILURE);
    } else {
        /*
         * case where world coordinates have been supplied 
         */
        printf(" reading world coordinates from input catalogue\n");

        dobjx = ml_calloc(nobjt, sizeof(double), &memory, "dobjx");
        dobjy = ml_calloc(nobjt, sizeof(double), &memory, "dobjy");
        radegs = ml_calloc(nobjt, sizeof(double), &memory, "radegs");
        decdegs = ml_calloc(nobjt, sizeof(double), &memory, "decdegs");
        rdegs = ml_calloc(nobjt, sizeof(double), &memory, "rdegs");
        ddegs = ml_calloc(nobjt, sizeof(double), &memory, "ddegs");
        offscale = ml_calloc(nobjt, sizeof(int), &memory, "offscale");
        rmag = ml_calloc(nobjt, sizeof(float), &memory, "rmag");
        sn = ml_calloc(nobjt, sizeof(float), &memory, "sn");

        /*
         * read in coordinates assuming in decimal WCS 
         */
        nobj2 = readpsfcat_simple_wcs(catname, radegs, decdegs, rmag);
        if (nobj2 != nobjt) {
            printf(" error reading catalogue \n");
            printf(" nobjt = %d, nobj2 = %d \n", nobjt, nobj2);
            exit(EXIT_FAILURE);
        }

        if (nobjt > 0) {
            /*
             * check coords are in sensible range 
             */
            for(i = 0; i < nobjt; i++) {
                if (radegs[i] < 0. || radegs[i] > 360. || decdegs[i] < -90.
                    || decdegs[i] > 90.) {
                    fprintf(stderr,
                            " these don't look like world coordinates \n");
                    fprintf(stderr, " %lf %lf \n", radegs[i], decdegs[i]);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    if (nobjt > 0) {
        printf(" %d objects read from PSF catalogue %s \n", nobjt, catname);
    } else {
      fflush(stdout);
      fprintf(stderr," no objects read from PSF catalogue %s \n", catname);
      exit(EXIT_FAILURE);
    }


    fitsize = nobjt;
    wfits = ml_calloc(pwidth * pheight, sizeof(double *), &memory, "wfits");
    for(i = 0; i < (pwidth * pheight); i++)
        wfits[i] = ml_calloc(fitsize, sizeof(double), &memory, "wfits[i]");

    xfits = ml_calloc(pwidth * pheight, sizeof(double *), &memory, "xfits");
    for(i = 0; i < (pwidth * pheight); i++)
        xfits[i] = ml_calloc(fitsize, sizeof(double), &memory, "xfits[i]");

    yfits = ml_calloc(pwidth * pheight, sizeof(double *), &memory, "yfits");
    for(i = 0; i < (pwidth * pheight); i++)
        yfits[i] = ml_calloc(fitsize, sizeof(double), &memory, "yfits[i]");

    zfits = ml_calloc(pwidth * pheight, sizeof(double *), &memory, "zfits");
    for(i = 0; i < (pwidth * pheight); i++)
        zfits[i] = ml_calloc(fitsize, sizeof(double), &memory, "zfits[i]");

    chip = ml_calloc(pwidth * pheight, sizeof(int *), &memory, "chip");
    for(i = 0; i < (pwidth * pheight); i++)
        chip[i] = ml_calloc(fitsize, sizeof(int), &memory, "chip[i]");

    nn = ml_calloc(pwidth * pheight, sizeof(int), &memory, "nn");

    /*
     * create the filename for the psf catalogue file 
     */

    strcpy(firstname, argv[1]);
    if (!(pstr = strrchr(firstname, *dotdelimiter)))
        pstr = firstname + strlen(firstname);
    sprintf(pstr, "%s", ".surfacepsf.fits");

    if (psfdir != NULL) {
        strcpy(fittedpsfname, psfdir);
        len = strlen(fittedpsfname);
        if (strncmp(&fittedpsfname[len - 1], "/", 1) != 0) {
            strncat(fittedpsfname, "/", 1);
        }
        if (access(fittedpsfname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr, " Can't write psf files to %s \n", fittedpsfname);
	  exit(EXIT_FAILURE);
        }
        len = strlen(firstname);
        strncat(fittedpsfname, firstname, len);
    } else {
        strcpy(fittedpsfname, firstname);
    }

    //      printf(" creating psf name %s \n",fittedpsfname);

    if (access(fittedpsfname, F_OK) == 0) {
        printf(" *** PSF file exists: removing \n");
        remove(fittedpsfname);
    }


    numstar = calloc(nchip, sizeof(int));
    if (numstar == NULL) {
      fflush(stdout);
      fprintf(stderr, " error allocating memory for numstar \n");
      exit(EXIT_FAILURE);
    }

    /*
     * got through each input image and accumulate the star postage stamps etc 
     */

    tot_nstar = 0;

    /*
     * define a reference image for the scale factor from the middle
     * of the field. The scalefactor for all the images in this
     * exposure will be set from this reference.  Note that other
     * exposures being processed elsewhere will have a different
     * scalefactor
     */
    refchip = nxchip / 2 + (nychip / 2) * nxchip;
    // set a default halfway through the list (shouldn't be needed as
    // long as the ref chip is in the list)
    refimage = nchip / 2;
    // if more than one chip, then try to find the reference chip 
    if (nchip>1)
      {
	for(image = 0; image < nchip; image++) {
	  /*
	   * strip out the chip number and look for the refimage chip number 
	   */
	  bzero(firstname, PATH_MAX);
	  strcpy(firstname, image_file[image]);
	  item = strtok(firstname, delims);
	  item = strtok(NULL, delims2);
	  i = atoi(item) - 1;
	  if (i == refchip) {
            refimage = image;
            break;
	  }
	}
      }

    for(kimage = 0; kimage < nchip; kimage++) {
        // start with the reference image
        image = kimage + refimage;
        if (image >= nchip) image -= nchip;

        /*
         * strip out the chip number and remember it 
         */
        bzero(firstname, PATH_MAX);
        strcpy(firstname, image_file[image]);

	// work out the chip ID 
	if (nchip>1)
	  {
	    // if multiple chips, attempt to read the chip number from the filename
	    item = strtok(firstname, delims);
	    item = strtok(NULL, delims2);
	    i = atoi(item) - 1;
	    chipnumber = i;
	  }
	else
	  {
	    // only one chip expected, so set chip ID to 0
	    chipnumber = i = 0;
	  }
        ychip = i / nxchip;
        xchip = i - ychip * nxchip;
        printf("\n image %s chip %d %d \n", image_file[image], xchip + 1,
               ychip + 1);
        if (i < 0 || i >= nchip_expected) {
	  fflush(stdout);
	  fprintf(stderr," error reading chip number for image %d \n",kimage+1);
	  fprintf(stderr," chip number %d out of %d x %d \n",i+1,nxchip,nychip);
	  exit(EXIT_FAILURE);
        }

        /*
         * strip off .fits extension on first argument if supplied 
         */
        len = strlen(image_file[image]);
        if (len <= 0) {
            fprintf(stderr, " error reading arguments from command line \n");
            exit(EXIT_FAILURE);
        }
        if (strncmp(image_file[image] + len - 5, ".fits", 5) == 0
            || strncmp(image_file[image] + len - 5, ".FITS", 5) == 0) {
            len = len - 5;
        }

        bzero(firstname, PATH_MAX);
        strncpy(firstname, image_file[image], len);
        bzero(shiftsname, PATH_MAX);
        strcpy(shiftsname, firstname);
        strcat(shiftsname, "_shifts.log");
        if ((shiftsfile = fopen(shiftsname, "w")) == NULL) {
            fprintf(stderr, " failed to make file %s \n", shiftsname);
            exit(EXIT_FAILURE);
        }

        bzero(usedstarsname, PATH_MAX);
        strcpy(usedstarsname, firstname);
        strcat(usedstarsname, "_stars.log");
        if ((usedstarsfile = fopen(usedstarsname, "w")) == NULL) {
            fprintf(stderr, " failed to make file %s \n", usedstarsname);
            exit(EXIT_FAILURE);
        }

        bzero(imageweightname, PATH_MAX);

        /*
         * get the image weight filename 
         */
        if (badpixdir != NULL) {
            // printf(" Path to good/bad pixel images %s \n",badpixdir);
            strcpy(imageweightname, badpixdir);
            len = strlen(imageweightname);
            if (strncmp(&imageweightname[len - 1], "/", 1) != 0) {
                strncat(imageweightname, "/", 1);
                len++;
            }
            len = strlen(image_file[image]);
            strncat(imageweightname, image_file[image], len);
        } else {
            strcpy(imageweightname, image_file[image]);
        }

        /*
         * strip off end of input image name so that weight suffix can be added 
         */
        if (!(pstr = strrchr(imageweightname, *delimiter)))
            pstr = imageweightname + strlen(imageweightname);
        sprintf(pstr, "%s", weight_suffix);

        if (GZIPPED > 0) {
            strcat(imageweightname, ".gz");
        }

        printf(" weight file %s \n", imageweightname);


        if (evariable != NULL) {
            strcpy(imagename, evariable);
            len = strlen(imagename);
            if (strncmp(&imagename[len - 1], "/", 1) != 0) {
                strncat(imagename, "/", 1);
            }
            len = strlen(firstname);
            strncat(imagename, firstname, len);
        } else {
            strcpy(imagename, firstname);
        }
        strncat(imagename, ".fits", 5);
	if (VERBOSE == 1)
	  printf(" reading image data from %s \n", imagename);

        /*
         * by default, look for the head files in the same place as
         * the data files but if the head environment variable is set,
         * use that location instead
         */

        bzero(headerfile, PATH_MAX);

        if (USE_SWARP == 1) {
            if (headevariable == NULL && evariable != NULL) {
                headevariable = evariable;
            }

            if (headevariable != NULL) {
                strcpy(headerfile, headevariable);
                len = strlen(headerfile);
                if (strncmp(&headerfile[len - 1], "/", 1) != 0) {
                    strncat(headerfile, "/", 1);
                }
                strcat(headerfile, image_file[image]);
            } else {
                strcpy(headerfile, image_file[image]);
            }
            strcpy(satname, prefs.sat_keyword);
        } else {
            strcpy(satname, "SATLEVEL");
        }

        /*
         * use fitsio routines to get the image info and check the
         * "bad ccd" keyword (we could do this with the swarp fits
         * routines instead but leave this as it is just now)
         */
        if (getimagesize
            (imagename, satname, dim, &gain, &satlev, &arcperpix, &angle,
             &badccd, &imageid, (float) minexposuretime)) {
            fflush(stdout);
            fprintf(stderr, " error checking input image %s \n", imagename);
            exit(EXIT_FAILURE);
        }

	if (imageid <= 0) imageid = 1+xchip+ychip*nxchip;

        if (imageid != (1 + xchip + ychip * nxchip)) {
            fflush(stdout);
            fprintf(stderr, "error in identifying chip number, %s. \n",
                    imagename);
            fflush(stderr);
            exit(EXIT_FAILURE);
        }

        if (badccd == 1) {
            printf(" %s flagged as bad, skipping this image \n", imagename);
            continue;
        }

        if (USE_SWARP == 1) {
            /*
             * create headerfile name needed for the swarp code 
             */
            if (!(pstr = strrchr(headerfile, *headdelimiter)))
                pstr = headerfile + strlen(headerfile);
            sprintf(pstr, "%s", prefs.head_suffix);
            /*
             * test it can be opened 
             */
            filep = fopen(headerfile, "r");
            if (filep == NULL) {
                fprintf(stderr, "error opening scamp file %s \n", headerfile);
                printf(" error opening scamp file %s \n", headerfile);
                printf
                    (" without these files no further astrometric correction can be applied \n");
                printf
                    (" and yet USE_SWARP correction has been specified in the code \n");
                printf
                    (" either define USE_SWARP 0 or supply a set of scamp head files, one per image \n");
                fflush(stdout);
                exit(EXIT_FAILURE);
            }
            fclose(filep);
        }


        imagesize = dim[0] * dim[1];

	if (VERBOSE == 1)
	  {
	    printf(" image dimensions = %d x %d \n", dim[0], dim[1]);

	    //      printf (" CCD gain = %6.1f e/ADU \n",gain);

	    printf(" CCD saturation level = %8.0f ADU \n", satlev);
	  }

        /*
         * set maximum pixel value as 50 percent of CCD saturation level
         * in case of CCD non-linearity 
         */
        satlimit = 0.5 * satlev;

	/* set the maximum fraction of the total flux which is allowed to be in the 
	   maximum pixel - designed to cut out any cosmic rays that have made it through the
	   various filtering processes */
	cosmicraylimit = 0.5;  // 0.5 is essentially no cut: for CFHTLenS the value was 0.15

        /*
         * set maximum intensity used for object detection 
         */
        fmax_intensity = satlimit;

	if (VERBOSE == 1)
	  printf(" trying to allocate %d MB for images and arrays\n",
		 (int) (15. * imagesize * 4. / 1024. / 1024.));

        region = ml_calloc(imagesize, sizeof(int), &memloop,
                           "Error allocating memory for region mask %d\n", i);

        /*
         * allocate memory for input image 
         */
        apix = ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for image \n");
        badpix = ml_calloc(imagesize, sizeof(float), &memloop,
                           " error allocating memory for image \n");

	if (SUBTRACT_MEDIAN == 1)
	  sortarray = ml_calloc(imagesize, sizeof(float), &memloop,
				" error allocating memory for sortarray\n");

	if (FILTER_WEIGHTS == 1)
	  {
	    opix = (float*)ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for image \n");
	    changed = (float*)ml_calloc(imagesize, sizeof(float), &memloop,
                         "Error allocating memory for changed array \n");
	    weightfilterarray = (float*)ml_calloc(1+2*fheight, sizeof(float), &memloop,
                         "Error allocating memory for weightfilterarray \n");
	  }

        /*
         * use swarp code to get the full wcs info 
         */
        rawcat = read_cat(imagename);
        QCALLOC(rawfield, fieldstruct, 1);
        load_dfield(rawfield, rawcat, headerfile, 0, 0);        // modified version of load_field
        wcs_raw = copy_wcs(rawfield->wcs);      // get the wcs info

        /*
         * read image data and measure noise level 
         */
        getdata_badpix(imagename, imageweightname, dim, apix, badpix, &noise);

	if (FILTER_WEIGHTS == 1)
	  {
	    nchanged = weightfilter(badpix, opix, changed, weightfilterarray, dim, fheight);
	    if (VERBOSE==1)
	      printf(" %d bad pixels changed to good \n",nchanged);
	    dilatemask(opix,dim);
	  }
	else
	  {
	    opix = badpix;
	    nchanged = 0;
	  }

	// dilate the bad pixel mask by one pixel in all directions to take care of 
	// left-over unmasked parts
	dilatemask(badpix,dim);

	if (VERBOSE == 1)
	  {
	    printf(" estimated noise in image = %f ADU \n", noise);
	  }

        if (noise <= 0.) {
	  fflush(stdout);
	  fprintf(stderr, " Erroneous noise value \n");
	  exit(EXIT_FAILURE);
        }

        /*
         * subtract assumed constant median background 
         */

        maxlevel = satlev / 2.;

        if (SUBTRACT_MEDIAN == 1) {
            mediansub(apix, badpix, dim, maxlevel, sortarray);
        }

        /*
         * identify objects within this image 
         */

        /*
         * set area limits 
         */
        lower_area_limit = 1;   // try to eliminate single pixels eg cosmic rays
        lower_merging_limit = 5;        // anything bigger than this will be kept unmerged
        arealimit = 100000;
        /*
         * set lower intensity threshold 
         */
        fintensity_limit = 3. * noise;

        /*
         * maximum object detection pixel value has already been set to half of
         * the sat level.  Now allocate arrays large enough to hold the 
         * digitised pixel values.  This must match with assumed array sizes
         * inside varylpthresholdf 
         */
        if (noise > 0.) {
            scale = 2. / noise;
            imax_intensity = (int) (fmax_intensity * scale);
        } else {
            printf(" noise value < 0 \n");
            exit(EXIT_FAILURE);
        }

        nobjects = varylpthresholdf2(apix, region, noise, dim[0], dim[1],
                                     lower_area_limit, lower_merging_limit,
                                     arealimit, fintensity_limit,
                                     fmax_intensity, imax_intensity);
	if (VERBOSE == 1)
	  printf(" %d objects identified within image \n", nobjects);

        boxsize = 200;

        if (boxsize < 10) {
            printf(" PSF box size too small, must be more than 10 pixels \n");
            exit(EXIT_FAILURE);
        } else {
	  if (VERBOSE == 1)
            printf(" PSF sampling interval %d pixels, %f arcmin \n", boxsize,
                   (float) boxsize * arcperpix / 60.);
        }

        xsize = (double) dim[0];
        ysize = (double) dim[1];

        xboxnum = 1 + (int) ((float) (dim[0] - 1) / (float) boxsize);
        yboxnum = 1 + (int) ((float) (dim[1] - 1) / (float) boxsize);

	if (VERBOSE == 1)
	  {
	    printf(" number of PSF boxes in x direction = %f rounded up to %d \n",
		   ((float) dim[0] / (float) boxsize), xboxnum);
	    printf(" number of PSF boxes in y direction = %f rounded up to %d \n",
		   ((float) dim[1] / (float) boxsize), yboxnum);
	  }

        maxboxno = xboxnum * yboxnum;

        /*
         * convert WCS to pixel positions using swarp WCS library 
         */
        for(ii = 0; ii < nobjt; ii++) {
            dobjx[ii] = -999.;
            dobjy[ii] = -999.;
            offscale[ii] = 1;
            wcspos[0] = radegs[ii];
            wcspos[1] = decdegs[ii];
            if (wcs_to_raw(wcs_raw, wcspos, rawpos) == RETURN_OK) {
                if (rawpos[0] > 0. && rawpos[0] < (double) rawfield->width) {
                    if (rawpos[1] > 0.
                        && rawpos[1] < (double) rawfield->height) {
                        dobjx[ii] = rawpos[0];
                        dobjy[ii] = rawpos[1];
                        offscale[ii] = 0;
                    }
                }
            }
            //  printf(" pos check %lf %lf %lf %lf \n",radegs[ii],decdegs[ii],dobjx[ii],dobjy[ii]);
        }

        /*
         * eliminate objects off the image 
         */
        /*
         * nobjt is the total number of stars in the catalogue
         * nobj is the number of stars that lie on this particular chip 
         */
        nobj = 0;
        for(i = 0; i < nobjt; i++) {
            if (offscale[i] == 0) {
                rdegs[nobj] = radegs[i];
                ddegs[nobj] = decdegs[i];
                objx[nobj] = (float) dobjx[i];
                objy[nobj] = (float) dobjy[i];
                mag[nobj] = rmag[i];
                // printf(" %lf %lf %f %f \n",radegs[i],decdegs[i],objx[nobj],objy[nobj]);
                nobj++;
            } else {
                // printf(" %lf %lf off-image \n",radegs[i],decdegs[i]);
            }
        }

        box = ml_calloc(nobj, sizeof(int), &memloop, "box");
        xshift = ml_calloc(nobj, sizeof(double), &memloop, "xshift");
        yshift = ml_calloc(nobj, sizeof(double), &memloop, "yshift");
        weight = ml_calloc(nobj, sizeof(double), &memloop, "weight");
        normvalue = ml_calloc(nobj, sizeof(double), &memloop, "normvalue");
        goodstar = ml_calloc(nobj, sizeof(int), &memloop, "goodstar");

        /*
         * allocate maximum possible memory for psf identifier 
         */

        boxno = ml_calloc(maxboxno, sizeof(int *), &memloop,
                          "Memory allocation error for psf identifier\n");
        for(i = 0; i < maxboxno; i++)
            boxno[i] = ml_calloc(nobj, sizeof(int), &memloop,
                                 "Memory allocation error for psf identifier[i]\n");

        mainboxno = ml_calloc(nobj, sizeof(int), &memloop,
                              "Memory allocation error for mainboxno \n");

        /*
         * allocate memory for postage stamp images 
         */
        dpix = ml_calloc(nobj, sizeof(double *), &memloop,
                         "Memory allocation error for sub-image pointers\n");
        for(i = 0; i < nobj; i++)
            dpix[i] = ml_calloc(pwidth * pheight, sizeof(double),
                                &memloop,
                                "Memory allocation error for sub-image \n");
        dbadpix = ml_calloc(nobj, sizeof(double *), &memloop,
                            "Memory allocation error for sub-image pointers\n");
        for(i = 0; i < nobj; i++)
            dbadpix[i] = ml_calloc(pwidth * pheight, sizeof(double),
                                   &memloop,
                                   "Memory allocation error for sub-image \n");

        /*
         * set limits on portion of image that a PSF must occupy (excludes
         * a narrow band around the edges) 
         */
        xmin = pwidth / 2;
        ymin = pheight / 2;
        xmax = dim[0] - pwidth / 2;
        ymax = dim[1] - pheight / 2;

        /*
         * read, normalise, fft and store all the psfs 
         */
        psf = ml_calloc(maxboxno, sizeof(double *), &memloop,
                        "Memory allocation error for psf pointers\n");
        npsfpixel = ml_calloc(maxboxno, sizeof(double *), &memloop,
                              "Memory allocation error for psf pointers\n");
        for(i = 0; i < maxboxno; i++) {
            psf[i] = ml_calloc(pwidth * pheight, sizeof(double), &memloop,
                               "Memory allocation error for psf[i] \n");
            npsfpixel[i] = ml_calloc(pwidth * pheight, sizeof(double),
                                     &memloop,
                                     "Memory allocation error for psf[i] \n");
        }

        Dpix = ml_calloc(maxboxno, sizeof(fftw_complex *), &memloop,
                         "Memory allocation for Dpix pointers.\n");
        for(i = 0; i < maxboxno; i++)
            Dpix[i] = ml_calloc(halfpwidth * pheight, sizeof(fftw_complex),
                                &memloop,
                                "Memory allocation error for Dpix \n");

        bad = ml_calloc(maxboxno, sizeof(int), &memloop, "memory for bad\n");
        num = ml_calloc(nobj, sizeof(int), &memloop, "memory for num\n");
        listnum = ml_calloc(maxboxno, sizeof(int), &memloop,
                            "memory for listnum\n");
        badnum = ml_calloc(maxboxno, sizeof(int), &memloop,
                           "memory for badnum\n");
        listx = ml_calloc(maxboxno, sizeof(float *), &memloop,
                          "memory for listx\n");
        listy = ml_calloc(maxboxno, sizeof(float *), &memloop,
                          "memory for listy\n");
        for(i = 0; i < maxboxno; i++) {
            listx[i] = ml_calloc(nobj, sizeof(float), &memloop,
                                 "memory for listx[%d]\n", i);
            listy[i] = ml_calloc(nobj, sizeof(float), &memloop,
                                 "memory for listy[%d]\n", i);
        }

        newpsf = ml_calloc(maxboxno, sizeof(double *), &memloop,
                           "memory for newpsf\n");
        count = ml_calloc(maxboxno, sizeof(double *), &memloop,
                          "memory for count\n");
        for(i = 0; i < maxboxno; i++) {
            newpsf[i] = ml_calloc(pwidth * pheight, sizeof(double),
                                  &memloop, "memory for newpsf[%d]\n", i);
            count[i] = ml_calloc(pwidth * pheight, sizeof(double),
                                 &memloop, "memory for count[%d]\n", i);
        }
        psfauto = ml_calloc(maxboxno, sizeof(double), &memloop,
                            "memory for psfauto\n");
        norm = ml_calloc(maxboxno, sizeof(double *), &memloop,
                         "memory for norm pointers\n");
        for(i = 0; i < maxboxno; i++)
            norm[i] = ml_calloc(nobj, sizeof(double), &memloop,
                                "memory for norm[%d]", i);

	if (VERBOSE == 1)
	  printf(" total memory allocated %d MB \n",
		 (int) ((memory + memloop) / 1024. / 1024.));

        totalbad = 0;
        totalgood = 0;

        /*
         * loop through each star from the catalogue, find which psf
         * it belongs to: then extract the data
         */

        t2 = time(NULL);
        //  printf (" Initial setup time was %g secs \n", difftime(t2, t1) );

        /*
         * perform two iterations over the stacking/shifting
         * operation, the first time shifting stars so their peaks
         * align, the second time so the cross-correlation with the
         * resulting stacked psf is a maximum (helps eliminate noise
         * bias in the psf stacking)
         */

        timereported = 0;

        /*
         * guess at number of good stars as being 90 percent of input number 
         */
        nstars = (int) ((float) nobj * 0.9);
        if (nstars < 1)
            nstars = 1;

        niter = 2;

        for(iterate = 1; iterate <= niter; iterate++) {
            printf("\n");
            t4 = time(NULL);

            /*
             * choose the smoothing length so that there are about
             * NAVERAGE stars in a box on the edge of the image 
             */
            boxsmooth =
                (int) sqrt(2. * (float) xsize * (float) ysize * NAVERAGE /
                           (float) nstars);

            if (boxsmooth < boxsize) {
                boxsmooth = boxsize;
            }

	    if (VERBOSE == 1)
	      printf(" PSF smoothing length %d pixels, %f arcmin \n", boxsmooth,
		     (float) boxsmooth * arcperpix / 60.);

            halfboxsmooth = (boxsmooth - boxsize) / 2;

            nstars = 0;
            num_snratio = 0;
            num_saturated = 0;
            num_off = 0;
            num_wrong_mag = 0;
            num_close = 0;
            num_badpix = 0;

            for(i = 0; i < maxboxno; i++) {
                listnum[i] = 0;
                badnum[i] = 0;
                bad[i] = 0;
                for(j = 0; j < (pwidth * pheight); j++) {
                    newpsf[i][j] = 0.;
                    count[i][j] = 0.;
                }
            }

            /*
             * start looping through list of PSF stars 
             */

            for(i = 0; i < nobj; i++) {
                /*
                 * initialise arrays 
                 */
                num[i] = 0;
                boxno[0][i] = -1;
                sn[i] = weight[i] = 0.;

                close = 0;

                /*
                 * test if close to another star in PSF catalogue
                 */
                for(d = 0; d < nobj; d++) {
                    if (d != i) {
                        if (fabs(objx[i] - objx[d]) < pwidth / 2
                            && fabs(objy[i] - objy[d]) < pheight / 2) {
                            close = 1;  // has a close neighbour
                        }
                    }
                }

                if (close == 1)
                    num_close++;

                if (close == 0) {

                    x = (int) (objx[i] + 0.5) - 1;
                    y = (int) (objy[i] + 0.5) - 1;

                    if (x >= xmin && x < xmax && y >= ymin && y < ymax) {
                        if (mag[i] >= brightmag && mag[i] <= faintmag) {
                            ximageoffset = 0.;
                            yimageoffset = 0.;
                            maxvalue = 0.;
                            maxlevel = satlev / 2.;
                            if (USE_SWARP == 0)
                                ngoodpix = extractpostagestamp
                                    (apix, badpix, dim, dpix[i], dbadpix[i],
                                     temp, badtemp, objx[i], objy[i],
                                     poserror, ximageoffset, yimageoffset,
                                     pwidth, pheight, noise, region,
                                     maxlevel);
                            else {
                                /*
                                 * if this is the first time
                                 * swarpextract has been called, the
                                 * scalefactor will be zero, so report
                                 * the current image.  Thereafter the
                                 * scalefactor will be set to the
                                 * central value in this image for all
                                 * subsequent calls to swarpextract
                                 */
                                if (scalefactor[0] <= 0.)
                                    printf
                                        (" using image %s for the scalefactor reference value\n",
                                         image_file[image]);
                                // use the swarp astrometry, but decide whether to also correct the shear distortion or not
                                if (CORRECT_DISTORTION == 1) {
                                    ngoodpix = tbswarpextract
                                        (wcs_raw, apix, badpix, dim, dpix[i],
                                         dbadpix[i], objpix, objx[i], objy[i],
                                         poserror, pwidth, pheight, ikernel,
                                         noise, region, maxlevel, scalefactor,
                                         rawpos, wcspos, wcsneg, wcscentre,
                                         kern, xfit, yfit, zfit, wfit, u, v);
                                    // printf(" thread %d image %d object %d ngoodpix %d \n",nt,i,ii,ngoodpix); fflush(stdout);
                                } else {
				  flagbad = 1; // specify that background objects are to be flagged
                                    ngoodpix = tbswarpextract_nodistortion
				      (wcs_raw, apix, badpix, opix, dim, dpix[i],
                                         dbadpix[i], objx[i], objy[i],
                                         poserror, pwidth, pheight,
                                         noise, region, scalefactor, maxlevel,
                                         rawpos, wcspos, wcsneg, wcscentre, flagbad);
                                    // printf(" thread %d image %d object %d ngoodpix %d \n",nt,i,ii,ngoodpix); fflush(stdout);
                                }
                            }

                            // printf(" %f %f %lf %d \n",objx[i],objy[i],dpix[i][0],ngoodpix);
                            if (ngoodpix > 0) {
                                psfnorm = noisevalue = maxvalue = distantmaxvalue = 0.;
                                for(iy = 0; iy < pheight; iy++) {
                                    yy = iy > pheight / 2 ? iy - pheight : iy;
                                    for(ix = 0; ix < pwidth; ix++) {
                                        xx = ix >
                                            pwidth / 2 ? ix - pwidth : ix;
                                        if (xx * xx + yy * yy <= poserror*poserror) 
					  {
                                            pixel = ix + iy * pwidth;
                                            psfnorm += dpix[i][pixel];
                                            noisevalue += noise * noise;
                                            if (dpix[i][pixel] > maxvalue)
                                                maxvalue = dpix[i][pixel];
					  }
					else
					  {
					    // check maximum pixel outside star area
                                            pixel = ix + iy * pwidth;
                                            if (dpix[i][pixel] > distantmaxvalue)
                                                distantmaxvalue = dpix[i][pixel];
					  }
                                    }
                                }
                                if (noisevalue <= 0.) {
                                    fflush(stdout);
                                    fprintf(stderr,
                                            " negative noise value \n");
                                    exit(EXIT_FAILURE);
                                }
                                noisevalue = sqrt(noisevalue);
                                /*
                                 * printf(" maxvalues %lf %lf %f %f %d %f %f \n",
                                 * rdegs[i],ddegs[i],objx[i],objy[i],chipnumber+1,maxvalue/psfnorm,psfnorm/noisevalue);
                                 */
                                if (psfnorm > (double) (snratio * noisevalue)) {
				  // check star isn't saturated nor a likely cosmic ray hit
				  //  and also that its maximum is significantly larger than any value outside poserror
                                    if (maxvalue <= satlimit
                                        && maxvalue < cosmicraylimit * psfnorm
					&& maxvalue > 2.*distantmaxvalue) {
                                        sn[i] = psfnorm / noisevalue;
                                        x1 = (x - halfboxsmooth) / boxsize;
                                        if (x1 < 0)
                                            x1 = 0;
                                        x2 = (x + halfboxsmooth) / boxsize;
                                        if (x2 >= xboxnum)
                                            x2 = xboxnum - 1;
                                        y1 = (y - halfboxsmooth) / boxsize;
                                        if (y1 < 0)
                                            y1 = 0;
                                        y2 = (y + halfboxsmooth) / boxsize;
                                        if (y2 >= yboxnum)
                                            y2 = yboxnum - 1;
                                        for(iy = y1; iy <= y2; iy++) {
                                            for(ix = x1; ix <= x2; ix++) {
                                                psfnum = iy * xboxnum + ix;
                                                boxno[num[i]][i] = psfnum;
                                                if (psfnum < 0
                                                    || psfnum >= maxboxno) {
                                                    printf
                                                        (" star %d num %d psfnum %d \n",
                                                         i, num[i], psfnum);
                                                    exit(EXIT_FAILURE);
                                                }
                                                num[i]++;
                                            }
                                        }
                                        ix = x / boxsize;
                                        iy = y / boxsize;
                                        psfnum = iy * xboxnum + ix;
                                        if (psfnum < 0 || psfnum >= maxboxno) {
                                            printf(" star %d psfnum %d \n", i,
                                                   psfnum);
                                            exit(EXIT_FAILURE);
                                        }
                                        mainboxno[i] = psfnum;
                                        nstars++;
                                    } else {
                                        num_saturated++;
                                    }
                                } else {
                                    num_snratio++;
                                }
                            } else {
                                num_badpix++;
                            }
                        } else {
                            num_wrong_mag++;
                        }
                    } else {
                        num_off++;
                    }
                }
            }

	    if (VERBOSE == 1)
	      printf(" iteration %d \n", iterate);

            if (iterate == 1) {
                printf(" excluded stars: \n");
                printf(" %d stars outside image area \n", num_off);
                printf(" %d stars outside mag range \n", num_wrong_mag);
                printf(" %d stars with no good data \n", num_badpix);
                if (snratio > 0.)
                    printf(" %d stars with peak S/N ratio < %8.1f \n",
                           num_snratio, snratio);
                printf(" %d stars with peak > %8.0f ADU \n", num_saturated,
                       satlimit);
                printf(" %d stars in a close pair\n\n", num_close);
                printf(" %d stars matching selection criteria \n", nstars);
            }

            for(psfnum = 0; psfnum < maxboxno; psfnum++) {
                if (iterate == 1) {
                    /*
                     * forward transform initial estimate of PSF,
                     * gaussian at origin sigma=1.5 pixel
                     */
                    sumb = 0.;
                    rsq0 = 1.5;
                    for(iy = 0; iy < pheight; iy++) {
                        for(ix = 0; ix < pwidth; ix++) {
                            pixel = ix + iy * pwidth;
                            xx = ix;
                            if (xx > pwidth / 2)
                                xx -= pwidth;
                            yy = iy;
                            if (yy > pheight / 2)
                                yy -= pheight;
                            rsq = (double) (xx * xx + yy * yy);
                            bpix[pixel] = exp(-rsq / (2. * rsq0 * rsq0));
                            sumb += bpix[pixel];
                    }}
                    for(iy = 0; iy < pheight; iy++) {
                        for(ix = 0; ix < pwidth; ix++) {
                            pixel = ix + iy * pwidth;
                            bpix[pixel] = bpix[pixel] / sumb;
                        }
                    }
                } else {
                    /*
                     * forward transform estimate of PSF from previous iteration 
                     */
                    for(iy = 0; iy < pheight; iy++) {
                        for(ix = 0; ix < pwidth; ix++) {
                            pixel = ix + iy * pwidth;
                            bpix[pixel] = psf[psfnum][pixel];
                        }
                    }
                }

                fftw_execute(qplan);

                /*
                 * store the psf transform in Dpix 
                 */
                for(iy = 0; iy < pheight; iy++) {
                    for(ix = 0; ix < halfpwidth; ix++) {
                        pixel = ix + iy * halfpwidth;
                        Dpix[psfnum][pixel] = Bpix[pixel];
                    }
                }
            }

            for(psfnum = 0; psfnum < maxboxno; psfnum++) {
                for(iy = 0; iy < pheight; iy++) {
                    for(ix = 0; ix < pwidth; ix++) {
                        pixel = ix + iy * pwidth;
                        psf[psfnum][pixel] = 0.;
                        npsfpixel[psfnum][pixel] = 0.;
                    }
                }
            }

	    if (VERBOSE == 1)
	      printf(" cross-correlating...\n");

            nstars1 = 0;

            timetest0 = time(NULL);

            for(i = 0; i < nobj; i++) {
                /*
                 * if this star has been found to be useful for making any psf... 
                 */
                if (num[i] > 0) {
                    nstars1++;

                    /*
                     * prepare to FFT the current star 
                     */
                    for(iy = 0; iy < pheight; iy++) {
                        for(ix = 0; ix < pwidth; ix++) {
                            pixel = ix + iy * pwidth;
                            bpix[pixel] = dpix[i][pixel];
                        }
                    }

                    /*
                     * fft it 
                     */
                    fftw_execute(qplan);

                    /*
                     * choose the psf box that this star lies in, if the local PSF
                     * is not zero then xcorrelate with the
                     * local PSF and store the result 
                     */

                    psfnum = mainboxno[i];

                    if (psfnum >= 0 && creal(Dpix[psfnum][0]) > 0.) {

                        xcorr_measure(padinv, Bpix, Dpix[psfnum],
                                      pheight, pwidth, padheight, padwidth,
                                      C, c, Cpad, shift);

                        for(iy = 0; iy < pheight; iy++) {
                            for(ix = 0; ix < halfpwidth; ix++) {
                                pixel = ix + iy * halfpwidth;
                                shiftpixft[pixel] = Bpix[pixel];
                            }
                        }

                        /*
                         * if the shift is a sensible amount, shift
                         * the data by the cross-correlation amount
                         */
                        if (fabs(shift[0]) < poserror
                            && fabs(shift[1]) < poserror) {
                            shiftft(pheight, pwidth, shiftpixft, shift);
                            fftw_execute(iqplan);
                            for(iy = 0; iy < pheight; iy++) {
                                for(ix = 0; ix < pwidth; ix++) {
                                    pixel = ix + iy * pwidth;
                                    dpix[i][pixel] =
                                        shiftpix[pixel] / (pwidth * pheight);
                                }
                            }
                            // normalise star 
                            dauto = 0.;
                            good = 1;
                            for(iy = 0; iy < pheight; iy++) {
                                yy = iy > pheight / 2 ? iy - pheight : iy;
                                for(ix = 0; ix < pwidth; ix++) {
                                    xx = ix > pwidth / 2 ? ix - pwidth : ix;
                                    dist1 = xx * xx + yy * yy;
                                    if (dist1 <= poserror*poserror) {
                                        pixel = ix + iy * pwidth;
                                        if (dbadpix[i][pixel] > 0.)
                                            dauto += dpix[i][pixel];
                                        else
                                            good = 0;
                                    }
                                }
                            }
                            if (good == 1 && dauto > 0.) {
                                for(iy = 0; iy < pheight; iy++) {
                                    for(ix = 0; ix < pwidth; ix++) {
                                        pixel = ix + iy * pwidth;
                                        dpix[i][pixel] =
                                            dpix[i][pixel] / dauto;
                                    }
                                }
                            } else {
                                for(iy = 0; iy < pheight; iy++) {
                                    for(ix = 0; ix < pwidth; ix++) {
                                        pixel = ix + iy * pwidth;
                                        dbadpix[i][pixel] = 0.;
                                        dpix[i][pixel] = 0.;
                                    }
                                }
                            }
                        } else  // otherwise disqualify this entire star 
                        {
                            for(iy = 0; iy < pheight; iy++) {
                                for(ix = 0; ix < pwidth; ix++) {
                                    pixel = ix + iy * pwidth;
                                    dbadpix[i][pixel] = 0.;
                                    dpix[i][pixel] = 0.;
                                }
                            }
                        }

                        xshift[i] = shift[0];
                        yshift[i] = shift[1];

                        /*
                         * printf(" iter %d star %d shifts %f %f \n",
                         * iterate,i,xshift[i],yshift[i]); 
                         */

                    }

                    /*
                     * go through every psf box this star is
                     * associated with and add it onto the psf with
                     * the shift calculated above.  Use a modified S/N
                     * weight to avoid single bright stars dominating
                     */

                    for(k = 0; k < num[i]; k++) {
                        psfnum = boxno[k][i];
                        if (psfnum >= 0) {
                            for(iy = 0; iy < pheight; iy++) {
                                for(ix = 0; ix < pwidth; ix++) {
                                    pixel = ix + iy * pwidth;
                                    if (dbadpix[i][pixel] > 0.) {
                                        psf[psfnum][pixel] +=
                                            dpix[i][pixel] / (50. * 50. /
                                                              (sn[i] *
                                                               sn[i]) + 1.);
                                        npsfpixel[psfnum][pixel]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }


            /*
             * for each psf, measure the star and psf autocorrelation
             * function in central region, and star/psf
             * cross-correlation, hence get correlation coefficient to
             * test whether each star is a good match to the main psf
             */

            for(psfnum = 0; psfnum < maxboxno; psfnum++) {
                psfauto[psfnum] = 0.;
                for(iy = 0; iy < pheight; iy++) {
                    for(ix = 0; ix < pwidth; ix++) {
                        pixel = ix + iy * pwidth;
                        count[psfnum][pixel] = 0.;
                        if (npsfpixel[psfnum][pixel] > 0.) {
                            psf[psfnum][pixel] =
                                psf[psfnum][pixel] / npsfpixel[psfnum][pixel];
                            if ((ix < cchwidth || ix > (pwidth - cchwidth))
                                && (iy < cchwidth
                                    || iy > (pheight - cchwidth))) {
                                psfauto[psfnum] +=
                                    psf[psfnum][pixel] * psf[psfnum][pixel];
                            }
                        } else {
                            psf[psfnum][pixel] = 0.;
                        }
                    }
                }
            }

	    printf(" cross-correlation results\n");

            for(i = 0; i < nobj; i++) {
                for(k = 0; k < num[i]; k++) {
                    psfnum = boxno[k][i];
                    dauto = 0.;
                    cross = 0.;
                    norm[k][i] = 0.;
                    if (psfnum >= 0) {
                        /*
                         * calculate cross-correlation between psf and
                         * star only in central region
                         */
                        for(iy = 0; iy < pheight; iy++) {
                            for(ix = 0; ix < pwidth; ix++) {
                                pixel = ix + iy * pwidth;
                                if (npsfpixel[psfnum][pixel] > 0.) {
                                    if ((ix < cchwidth
                                         || ix > (pwidth - cchwidth))
                                        && (iy < cchwidth
                                            || iy > (pheight - cchwidth))) {
                                        dauto +=
                                            dpix[i][pixel] * dpix[i][pixel];
                                    }
                                }
                            }
                        }
                        for(iy = 0; iy < pheight; iy++) {
                            for(ix = 0; ix < pwidth; ix++) {
                                pixel = ix + iy * pwidth;
                                if (npsfpixel[psfnum][pixel] > 0.) {
                                    if ((ix < cchwidth
                                         || ix > (pwidth - cchwidth))
                                        && (iy < cchwidth
                                            || iy > (pheight - cchwidth))) {
				      cross +=
					dpix[i][pixel]*psf[psfnum][pixel] / 
					sqrt(dauto*psfauto[psfnum]);
				      norm[k][i] +=
					dpix[i][pixel]*psf[psfnum][pixel] /
					psfauto[psfnum];
                                    }
                                }
                            }
                        }

			printf("%d %lf %lf %f %f %f \n",i+1,rdegs[i],ddegs[i],objx[i],objy[i],cross);

                        /*
                         * deselect stars with low cross-correlation values 
                         */
                        if (cross < 0.85) {
                            badnum[psfnum]++;
                            boxno[k][i] = -1;
                        }
                    }
                }
            }


            /*
             * make a new stack of psf stars with bad ones removed 
             */

            for(i = 0; i < nobj; i++) {
                starweight = 1. / (1. + 50. * 50. / (sn[i] * sn[i]));
                for(k = 0; k < num[i]; k++) {
                    psfnum = boxno[k][i];
                    if (psfnum >= 0) {

                        /*
                         * store record of which stars have been used
                         * for this psf
                         */

                        listx[psfnum][listnum[psfnum]] = objx[i];
                        listy[psfnum][listnum[psfnum]] = objy[i];
                        listnum[psfnum]++;

                        /*
                         * co-add shifted stars to make new psfs 
                         */

                        for(iy = 0; iy < pheight; iy++) {
                            for(ix = 0; ix < pwidth; ix++) {
                                pixel = ix + iy * pwidth;

                                if (dbadpix[i][pixel] > 0.) {
                                    newpsf[psfnum][pixel] +=
                                        dpix[i][pixel] * starweight;
                                    count[psfnum][pixel] += starweight;
                                    totalgood++;
                                }

                            }
                        }
                    }
                }
            }

            /*
             * now deal with stacked PSFs - normalise 
             */

            for(psfnum = 0; psfnum < maxboxno; psfnum++) {
                bad[psfnum] = 0;
                for(i = 0; i < (pwidth * pheight); i++) {
                    if (count[psfnum][i] > 0.) {
                        psf[psfnum][i] = newpsf[psfnum][i] / count[psfnum][i];
                    } else {
                        bad[psfnum]++;
                    }
                }
                if (bad[psfnum] > 0) {
                    /*
                     * printf ("insufficient data to form psf %d \n",psfnum+1);
                     */
                    listnum[psfnum] = 0;
                    for(i = 0; i < (pwidth * pheight); i++) {
                        psf[psfnum][i] = 0.;
                    }
                } else {
                    /*
                     * normalise using central values
                     */
                    normvalue[0] = 0.;
                    for(y = 0; y < pheight; y++) {
                        yy = y > pheight / 2 ? y - pheight : y;
                        for(x = 0; x < pwidth; x++) {
                            xx = x > pwidth / 2 ? x - pwidth : x;
                            pixel = x + y * pwidth;
                            dist1 = xx * xx + yy * yy;
                            if (dist1 <= poserror*poserror)
                                normvalue[0] += psf[psfnum][pixel];
                        }
                    }
                    for(iy = 0; iy < pheight; iy++) {
                        for(ix = 0; ix < pwidth; ix++) {
                            pixel = ix + iy * pwidth;
                            psf[psfnum][pixel] =
                                psf[psfnum][pixel] / normvalue[0];
                        }
                    }
                }
            }

            /*
             * printf (" %f percent of pixels excluded as contaminated by other objects \n",
             * (100.*(float)totalbad/(float)(totalbad+totalgood)));
             */

            t5 = time(NULL);
            /*
             * printf (" time per iteration was %g secs \n", difftime(t5, t4) );
             */
            /*
             * write the PSF boxes for this image
             */
	    if (WRITE_INDIVIDUAL_SURFACE==1)
	      {
		if (write_psfsurface(psfdir, firstname, psf, xboxnum, yboxnum,
				     pwidth, pheight, listnum, listx, listy,
				     boxsize, boxsmooth) != 0) {
		  fprintf(stderr, " error writing PSF surface file.\n");
		  exit(EXIT_FAILURE);
		}
	      }

            /*
             * next iteration 
             */
        }

        /*
         * print out statistics and warnings 
         */

        t3 = time(NULL);
	if (VERBOSE==1)
	  printf(" psf creation time was %g secs \n", difftime(t3, t1));

        fprintf(usedstarsfile,
                " statistics of numbers of stars used in each sub-field PSF:\n");
        fprintf(usedstarsfile, "    PSF  n(x)  n(y)  number  number\n");
        fprintf(usedstarsfile, "                      used  rejected\n");
        for(i = 0; i < maxboxno; i++) {
            y = i / xboxnum;
            x = i - xboxnum * y;
            x++;
            y++;
            fprintf(usedstarsfile, " %6d %5d %5d %6d %6d\n", i + 1, x, y,
                    listnum[i], badnum[i]);
        }

        warning = 0;
        for(i = 0; i < maxboxno; i++) {
            if (listnum[i] > 0 && badnum[i] > listnum[i] && badnum[i] > 3) {
                warning = 1;
            }
        }

        if (warning == 1) {
            printf
                (" *** WARNING large fraction of stars rejected for at least one PSF \n");
            printf
                (" *** I strongly suggest review of star selection criteria \n");
        }

        /*
         * warning=0;
         * for (i=0; i<maxboxno; i++)
         * {
         * if (listnum[i] < PSF_STAR_LIMIT)
         * {
         * printf(" *** WARNING low numbers of stars selected for PSF subregion number %d\n",i+1);
         * warning=1;
         * }
         * }
         * 
         * if (warning==1)
         * {
         * printf (" *** these regions will be flagged as having no PSF \n");
         * printf (" *** I strongly suggest review of star selection criteria and checking of PSF \n");
         * }
         */

        /*
         * now accumulate fitting data for the good, shifted stars 
         */
        nstar = 0;
        for(i = 0; i < nobj; i++) {
            starweight = 1. / (1. + 50. * 50. / (sn[i] * sn[i]));
            normvalue[i] = weight[i] = 0.;
            x = (int) (objx[i] + 0.5) - 1;
            y = (int) (objy[i] + 0.5) - 1;
            ix = x / boxsize;
            iy = y / boxsize;
            psfnum = iy * xboxnum + ix;
            good = 0;
            box[i] = -1;
            for(k = 0; k < num[i]; k++) {
                if (psfnum == boxno[k][i]) {
                    good = 1;
                    box[i] = psfnum;
                }
            }
            if (good == 1) {
                /*
                 * this star has been accepted as good in the box it occupies 
                 */
                /*
                 * normalise it  - measure normalisation only within poserror 
                 */
                for(y = 0; y < pheight; y++) {
                    yy = y > pheight / 2 ? y - pheight : y;
                    for(x = 0; x < pwidth; x++) {
                        xx = x > pwidth / 2 ? x - pwidth : x;
                        dist1 = xx * xx + yy * yy;
                        if (dist1 <= poserror*poserror) {
                            pixel = x + y * pwidth;
                            if (dbadpix[i][pixel] <= 0.) {
                                good = 0;
                                normvalue[i] = 0.;
                                break;
                            }
                            normvalue[i] += dpix[i][pixel];
                        }
                    }
                }

                if (normvalue[i] <= 0.)
                    good = 0;

                if (good == 1 && normvalue[i] > 0.) {
                    /*
                     * divide by sum of data values excluding bad pixels - note that
                     * if a bad pixel occurs in a region of significant signal the
                     * normalisation will be biased, so this effectively assumes that
                     * bad pixels occur only in the sky background and not on the star
                     * itself - therefore assumes stars have already been filtered out
                     * if they are significantly affected by bad pixels 
                     */

                    for(j = 0; j < pwidth * pheight; j++) {
                        dpix[i][j] = dpix[i][j] / normvalue[i];
                    }
                    goodstar[nstar] = i;
                    nstar++;

                    /*
                     * define global x,y position in range -1 < val < 1 
                     */
                    xval = (objx[i] + xchip * xchipsampling) / hxsize - 1.;
                    yval = objy[i] + ychip * ychipsampling;
                    if (ychip > 0)
                        yval += big_gap;
                    if (ychip > 2)
                        yval += big_gap;
                    yval = yval / hysize - 1.;
                    // printf(" star position %lf %lf %d %d \n",xval,yval,xchip,ychip);

                    /*
                     * store the global information 
                     */
                    if (box[i] >= 0) {
                        for(j = 0; j < pwidth * pheight; j++) {
                            if (dbadpix[i][j] > 0.) {
                                /*
                                 * remember which chip we're on 
                                 */
                                chip[j][nn[j]] = xchip + nxchip * ychip;
                                wfits[j][nn[j]] = starweight;
                                /*
                                 * fill fit arrays 
                                 */
                                xfits[j][nn[j]] = xval;
                                yfits[j][nn[j]] = yval;
                                zfits[j][nn[j]] = dpix[i][j];
                                /*
                                 * count total number of useful values for this pixel 
                                 */
                                nn[j]++;
                                /*
                                 * store the star weight for output information only
                                 * (note only one weight stored for all pixels, as long as one
                                 * pixel at least is valid) 
                                 */
                                weight[i] = starweight;
                            }
                        }
                    }
                }
            }
        }

        if (OUTPUT_STARS == 1) {
            /*
             * write the good stars out to a fits file 
	     * NB careful to use rdegs,ddegs array not radegs,decdegs
             */
            writestars(pwidth, pheight, nstar, dbadpix, dpix, goodstar,
                       rdegs, ddegs, objx, objy, sn, psfdir, firstname,
                       imagename);
        }

        for(i = 0; i < nobj; i++) {
            if (WCS == 1) {
                fprintf(shiftsfile,
                        " %10.5lf %10.5lf %8.2lf %8.2lf %6.2lf %6.2lf %8.4lf \n",
                        rdegs[i], ddegs[i], objx[i], objy[i], xshift[i],
                        yshift[i], weight[i]);
            } else {
                fprintf(shiftsfile, " %8.2lf %8.2lf %6.2lf %6.2lf %8.4lf \n",
                        objx[i], objy[i], xshift[i], yshift[i], weight[i]);
            }
        }

        if (imageid > 0) {
            numstar[imageid - 1] = nstar;
        } else {
            printf(" Image chip ID <= 0 \n ");
        }
        printf(" image %d, number of good stars to fit surface = %d \n",
               imageid, nstar);
	fflush(stdout);

        tot_nstar += nstar;

        /*
         * free memory ready for next image 
         */

        free(region);
        free(apix);
        free(badpix);
	if (FILTER_WEIGHTS == 1)
	  {
	    free(opix);
	    free(changed);
	    free(weightfilterarray);
	  }
	if (SUBTRACT_MEDIAN == 1)
	  free(sortarray);
        free(box);
        free(xshift);
        free(yshift);
        free(weight);
        free(normvalue);
        free(goodstar);
        for(i = 0; i < maxboxno; i++) {
            free(boxno[i]);
            free(psf[i]);
            free(npsfpixel[i]);
            free(Dpix[i]);
            free(listx[i]);
            free(listy[i]);
            free(newpsf[i]);
            free(count[i]);
            free(norm[i]);
        }
        free(Dpix);
        free(npsfpixel);
        free(psf);
        free(boxno);
        free(mainboxno);
        free(newpsf);
        free(count);
        free(psfauto);
        free(norm);
        for(i = 0; i < nobj; i++)
            free(dpix[i]);
        free(dpix);
        for(i = 0; i < nobj; i++)
            free(dbadpix[i]);
        free(dbadpix);
        free(bad);
        free(num);
        free(listnum);
        free(badnum);
        free(listx);
        free(listy);

        free_cat(&rawcat, 1);
        free(rawfield);

        //printf(" memory freed ready for next image\n");

        memloop = 0.;

        fclose(shiftsfile);
        fclose(usedstarsfile);

        /*
         * next image 
         */
    }

    if (tot_nstar <= 0) {
        printf(" no useful stars found in mosaic, no PSF created \n");
        exit(0);
    }

    /*
     * now fit the global polynomial surface 
     */

    for(y = 0; y < pheight; y++) {
        for(x = 0; x < pwidth; x++) {
            pixel = x + y * pwidth;
            // printf (" fitting pixel %d %d \n",x,y);
            /*
             * normalise the weights 
             */
            sumw = 0.;
            for(i = 0; i < nn[pixel]; i++)
                sumw += wfits[pixel][i];
            for(i = 0; i < nn[pixel]; i++)
                wfits[pixel][i] = wfits[pixel][i] / sumw;

            globalsvdfit(xfits[pixel], yfits[pixel], zfits[pixel],
                         wfits[pixel], chip[pixel], nn[pixel], nchip, order,
                         crossterm, chipvariation, chiporder, ncoeffs, avals,
                         u, v, w);
            /*
             * put fit coefficients into global pixel array 
             */
            for(j = 0; j < ncoeffs; j++) {
                acoeffs[pixel][j] = avals[j + 1];       /* shift by 1 from NR routine */
            }
        }
    }

    printf(" surface coefficients calculated \n");

    /*
     * write out a 3D image of the fitted pixel values across the chip 
     */

    yboxnum = 7 * nychip;       // 7*4
    xboxnum = 3 * nxchip;       // 3*9

    psf = (double **) calloc(xboxnum * yboxnum, sizeof(double *));
    listx = (float **) calloc(xboxnum * yboxnum, sizeof(float *));
    listy = (float **) calloc(xboxnum * yboxnum, sizeof(float *));
    listnum = (int *) calloc(xboxnum * yboxnum, sizeof(int));

    psfnum = 0;
    int iix, iiy;
    for(iy = 0; iy < nychip; iy++) {
        for(iiy = 0; iiy < 7; iiy++) {
            for(ix = 0; ix < nxchip; ix++) {
                for(iix = 0; iix < 3; iix++) {
                    ichip = ix + iy * nxchip;
                    xval =
                        ix * xchipsampling + iix * dim[0] / 3 + dim[0] / 6;;
                    xval = xval / hxsize - 1.;
                    yval =
                        iy * ychipsampling + iiy * dim[1] / 7 + dim[1] / 14;
                    if (iy > 0)
                        yval += big_gap;
                    if (iy > 2)
                        yval += big_gap;
                    yval = yval / hysize - 1.;
                    psf[psfnum] =
                        (double *) calloc(pwidth * pheight, sizeof(double));
                    listx[psfnum] = (float *) calloc(1, sizeof(float));
                    listy[psfnum] = (float *) calloc(1, sizeof(float));
                    for(pixel = 0; pixel < (pwidth * pheight); pixel++) {
                        globalreconstruct(xval, yval, order, crossterm,
                                          chipvariation, chiporder, ichip,
                                          nchip, acoeffs[pixel],
                                          &psf[psfnum][pixel]);
                    }
                    psfnum++;
                }
            }
        }
    }


    /*
     * check the normalisation 
     */
    check_psf_norm(psf, psfnum, xboxnum, yboxnum, pwidth, pheight);

    /*
     * add a default star in to satisfy the next routine 
     */
    for(iy = 0; iy < yboxnum; iy++) {
        for(ix = 0; ix < xboxnum; ix++) {
            psfnum = ix + iy * xboxnum;
            if (listnum[psfnum] <= 0) {
                listnum[psfnum] = 0;
                listx[psfnum][listnum[psfnum]] = 0.;
                listy[psfnum][listnum[psfnum]] = 0.;
                listnum[psfnum]++;
            }
        }
    }

    if (writepsf_fits
        (fittedpsfname, psf, xboxnum, yboxnum, pwidth, pheight, listnum,
         listx, listy, boxsize, boxsmooth) != 0) {
        printf(" error writing psf \n");
        exit(EXIT_FAILURE);
    }


/*  write out polynomial coefficients to FITS file */

    if (crossterm == 1) {
        order = -order;
    }
    int correct_distortion = (int)CORRECT_DISTORTION;
    if (writeglobalcoeffs
        (coeffname, correct_distortion, acoeffs, ncoeffs, order, chipvariation, chiporder,
         scalefactor, pwidth, pheight, nxchip, nychip, numstar, xchipsampling,
         ychipsampling, big_gap, hxsize, hysize, poserror) != 0) {
        printf(" error writing psf \n");
        exit(EXIT_FAILURE);
    }

    return 0;
}

int
check_psf_norm(double **psf, int psfnum, int xboxnum, int yboxnum, int pwidth,
               int pheight)
{
    int ix, iy, pixel;
    int x, y, xx, yy, dist1;
    int psfwarning = 0;
    double psfnorm;
    double minpsfnorm = 100.;
    double maxpsfnorm = 0.;
    double minpsfnormlimit = 0.95;
    double maxpsfnormlimit = 1.05;

    for(iy = 0; iy < yboxnum; iy++) {
        for(ix = 0; ix < xboxnum; ix++) {
            psfnum = ix + iy * xboxnum;
            psfnorm = 0.;
            psfwarning = 0;
            for(y = 0; y < pheight; y++) {
                yy = y > pheight / 2 ? y - pheight : y;
                for(x = 0; x < pwidth; x++) {
                    xx = x > pwidth / 2 ? x - pwidth : x;
                    dist1 = xx * xx + yy * yy;
                    if (dist1 <= poserror*poserror) {
                        pixel = x + y * pwidth;
                        psfnorm += psf[psfnum][pixel];
                    }
                }
            }
            if (psfnorm < minpsfnormlimit) {
                printf(" WARNING PSF normalisation in bin %d %d = %lf \n", ix,
                       iy, psfnorm);
                psfwarning = 1;
            }
            if (psfnorm > maxpsfnormlimit) {
                printf(" WARNING PSF normalisation in bin %d %d = %lf \n", ix,
                       iy, psfnorm);
                psfwarning = 1;
            }
            if (psfwarning == 0) {
                if (psfnorm < minpsfnorm)
                    minpsfnorm = psfnorm;
                if (psfnorm > maxpsfnorm)
                    maxpsfnorm = psfnorm;
            }
        }
    }
    printf(" PSF normalisation min and max of decent locations: %lf %lf \n",
           minpsfnorm, maxpsfnorm);

    return psfwarning;
}


int
write_psfsurface(char *psfdir, char *firstname, double **f, int xboxnum,
                 int yboxnum, int pwidth, int pheight, int *listnum,
                 float **listx, float **listy, int boxsize, int boxsmooth)
{
    int len;
    char outname[PATH_MAX];


    if (psfdir != NULL) {
        strcpy(outname, psfdir);
        len = strlen(outname);
        if (strncmp(&outname[len - 1], "/", 1) != 0) {
            strncat(outname, "/", 1);
        }
        if (access(outname, F_OK) != 0) {
            printf(" Can't write psf files to %s \n", outname);
            exit(EXIT_FAILURE);
        }
        len = strlen(firstname);
        strncat(outname, firstname, len);
    } else {
        strcpy(outname, firstname);
    }
    strcat(outname, ".psfsurface.fits");

    if (access(outname, F_OK) == 0) {
        printf(" *** surface psf file file %s exists: removing \n", outname);
        remove(outname);
    }
    return writepsf_fits(outname, f, xboxnum, yboxnum, pwidth, pheight,
                         listnum, listx, listy, boxsize, boxsmooth);
}


void
writestars(int pwidth, int pheight, int nstar, double **dbadpix,
           double **dpix, int *goodstar, double *radegs, double *decdegs,
           float *objx, float *objy, float *sn, const char *psfdir,
           const char *firstname, char *imagename)
{
    int psize, testsize;
    int i, j, x, y;
    int xx, yy, pixel;
    int len;
    long lrow, anaxis = 3;
    long anaxes[3];
    long fpixel[3] = { 1, 1, 1 };
    double *otestout;
    char *stampname, *pstr;
    fitsfile *tfptr;
    char pdelimiter[] = { "/" };
    int testfields = 5;         /* number of table columns in table */
    char testextname[] = "star info";   /* extension name */
    char *testtype[] = { "RA", "DEC", "X", "Y", "SN" }; /* names of columns */
    char *testform[] = { "D", "D", "E", "E", "E" };     /* data types of columns */
    char *testunit[] = { "degs", "degs", "pixels", "pixels", "ratio" }; /* units */
    int status = 0;

    psize = pwidth * pheight;
    testsize = nstar * psize;
    otestout = (double *) calloc(testsize, sizeof(double));
    stampname = (char *) calloc(PATH_MAX, sizeof(char));
    if (!stampname || !otestout) {
        fprintf(stderr, "Not enough memory to write usedstars\n");
        return;
    }
    for(i = 0; i < nstar; i++) {
        for(y = 0; y < pheight; y++) {
            for(x = 0; x < pwidth; x++) {
                j = x + y * pwidth;
                xx = x - pwidth / 2;
                yy = y - pheight / 2;
                if (xx < 0)
                    xx += pwidth;
                if (yy < 0)
                    yy += pheight;
                pixel = (xx + yy * pwidth) + i * psize;
                if (dbadpix[goodstar[i]][j] > 0.)
                    otestout[pixel] = dpix[goodstar[i]][j];
                else
                    otestout[pixel] = -0.1;
            }
        }
    }

    anaxes[0] = pwidth;
    anaxes[1] = pheight;
    anaxes[2] = nstar;


    /*
     * create the filename for the fits image 
     */
    if (psfdir != NULL) {
        strcpy(stampname, psfdir);
        len = strlen(stampname);
        if (strncmp(&stampname[len - 1], "/", 1) != 0) {
            strncat(stampname, "/", 1);
        }
        if (access(stampname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr," Can't write psf files to %s \n", stampname);
            exit(EXIT_FAILURE);
        }
        len = strlen(firstname);
        strncat(stampname, firstname, len);
    } else {
        strcpy(stampname, firstname);
    }
    strcat(stampname, "_usedstars.fits");

    if (access(stampname, F_OK) == 0) {
        printf(" *** postage stamps file exists: removing \n");
        remove(stampname);
    }

    status = 0;
    /*
     * create the file and write the postage stamps array 
     */
    fits_create_file(&tfptr, stampname, &status);
    fits_create_img(tfptr, -32, anaxis, anaxes, &status);
    if (status) {
        printf(" error writing postagestamp test pixel data %s \n",
               stampname);
        fits_report_error(stderr, status);      /* print error message */
    }
    if (fits_write_pix(tfptr, TDOUBLE, fpixel, testsize, otestout, &status)) {
        printf(" error writing postagestamp test pixel data %s \n",
               stampname);
        fits_report_error(stderr, status);      /* print error message */
    }

    if (!(pstr = strrchr(imagename, *pdelimiter)))
        pstr = imagename;
    if (fits_update_key
        (tfptr, TSTRING, "IMGNAME", pstr, "image name", &status)) {
        printf(" error writing keyword into %s \n", stampname);
        fits_report_error(stderr, status);      /* print error message */
    }
    if (fits_update_key(tfptr, TINT, "PWIDTH", &pwidth, "width", &status)) {
        printf(" error writing keyword into %s \n", stampname);
        fits_report_error(stderr, status);      /* print error message */
    }
    if (fits_update_key(tfptr, TINT, "PHEIGHT", &pwidth, "height", &status)) {
        printf(" error writing keyword into %s \n", stampname);
        fits_report_error(stderr, status);      /* print error message */
    }
    if (fits_update_key(tfptr, TINT, "NSTAR", &pwidth, "number", &status)) {
        printf(" error writing keyword into %s \n", stampname);
        fits_report_error(stderr, status);      /* print error message */
    }
    fits_create_tbl(tfptr, BINARY_TBL, nstar, testfields, testtype,
                    testform, testunit, testextname, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        printf(" error writing table into %s \n", stampname);
    }
    for(i = 0; i < nstar; i++) {
        lrow = i + 1;
        fits_write_col(tfptr, TDOUBLE, 1, lrow, 1, 1,
                       &radegs[goodstar[i]], &status);
        fits_write_col(tfptr, TDOUBLE, 2, lrow, 1, 1,
                       &decdegs[goodstar[i]], &status);
        fits_write_col(tfptr, TFLOAT, 3, lrow, 1, 1,
                       &objx[goodstar[i]], &status);
        fits_write_col(tfptr, TFLOAT, 4, lrow, 1, 1,
                       &objy[goodstar[i]], &status);
        fits_write_col(tfptr, TFLOAT, 5, lrow, 1, 1, &sn[goodstar[i]],
                       &status);
    }
    fits_close_file(tfptr, &status);
    if (status) {
        printf(" error writing postagestamp test pixel data %s \n",
               stampname);
        fits_report_error(stderr, status);      /* print error message */
    }
    free(otestout);
    free(stampname);
    status = 0;
}


char *
create_psf_name(const char *psfdir, char *firstname, const char *dotdelimiter)
{
    char *pstr;
    int len;
    char *coeffname;

    coeffname = (char *) calloc(PATH_MAX, sizeof(char));

    if (!(pstr = strrchr(firstname, *dotdelimiter)))
        pstr = firstname + strlen(firstname);
    sprintf(pstr, "%s", ".psfcoeffs.fits");

    if (psfdir != NULL) {
        strcpy(coeffname, psfdir);
        len = strlen(coeffname);
        if (strncmp(&coeffname[len - 1], "/", 1) != 0) {
            strncat(coeffname, "/", 1);
        }
        if (access(coeffname, F_OK) != 0) {
	  fflush(stdout);
	  fprintf(stderr," Can't write psf files to %s \n", coeffname);
            exit(EXIT_FAILURE);
        }
        strcat(coeffname, firstname);
    } else {
        strcpy(coeffname, firstname);
    }
    printf(" creating psf name %s \n", coeffname);

    if (access(coeffname, F_OK) == 0) {
        printf(" *** PSF file exists: removing \n");
        remove(coeffname);
    }
    return coeffname;
}



void
getdata_badpix(char *imagename, char *badpixelname, int *dim, float *apix,
               float *badpix, float *mednoiseval)
{
    fitsfile *afptr, *bfptr;    /* FITS file pointers */
    int status = 0;             /* CFITSIO status value must be initialized to zero */
    int anaxis, bnaxis;
    long anaxes[2] = { 1, 1 }, fpixel[2] = {
    1, 1}, bnaxes[2] = {
    1, 1};
    int size, bsize;
    int i, n, ix, iy, ymin, ymax, xmin, xmax;
    int x, y, medbin;
    float sortval[10201];
    float median, quartile, noiseval[32];

/* open input image */
    fits_open_file(&afptr, imagename, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(EXIT_FAILURE);
    }

/* read dimensions */
    fits_get_img_dim(afptr, &anaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(EXIT_FAILURE);
    }

    if (anaxis != 2) {
        printf
	  ("Error: %d image axes found: images with other than 2 dimensions are not supported\n",anaxis);
        exit(EXIT_FAILURE);
    }

    if (dim[0] != anaxes[0]) {
        fprintf(stderr, " error reading image dimensions \n");
        exit(EXIT_FAILURE);
    }
    if (dim[1] != anaxes[1]) {
        fprintf(stderr, " error reading image dimensions \n");
        exit(EXIT_FAILURE);
    }

    size = dim[0] * dim[1];
    bsize = dim[0] * dim[1];

/* read input data into image array */

    if (fits_read_pix(afptr, TFLOAT, fpixel, size, NULL, apix, NULL, &status)) {
        printf(" error reading pixel data \n");
        exit(EXIT_FAILURE);
    }

/* close main image file */

    fits_close_file(afptr, &status);

    if (status) {
        fits_report_error(stderr, status);      /* print error message */
        exit(EXIT_FAILURE);
    }


/* repeat for the bad pixel mask (which may not exist) */

/* open input images */
    fits_open_file(&bfptr, badpixelname, READONLY, &status);

    if (status == 0) {

/* read dimensions */
        fits_get_img_dim(bfptr, &bnaxis, &status);
        fits_get_img_size(bfptr, 2, bnaxes, &status);

        if (status) {
            fits_report_error(stderr, status);  /* print error message */
            exit(EXIT_FAILURE);
        }

        if (bnaxis != 2) {
	  printf("Error: %d weight image axes found: images with other than 2 dimensions are not supported\n",bnaxis);
            exit(EXIT_FAILURE);
        }

        if (dim[0] != bnaxes[0]) {
            fprintf(stderr, " error reading bad pixel image dimensions \n");
            exit(EXIT_FAILURE);
        }
        if (dim[1] != bnaxes[1]) {
            fprintf(stderr, " error reading bad pixel image dimensions \n");
            exit(EXIT_FAILURE);
        }

/* read input data into image array */

        if (fits_read_pix(bfptr, TFLOAT, fpixel, bsize, NULL, badpix,
                          NULL, &status)) {
            printf(" error reading bad pixel mask \n");
            fits_report_error(stderr, status);  /* print error message */
            exit(EXIT_FAILURE);
        }

        fits_close_file(bfptr, &status);

        if (status) {
            fits_report_error(stderr, status);  /* print error message */
            exit(EXIT_FAILURE);
        }

    } else {
        if (WEIGHTS_REQUIRED > 0) {
            fprintf(stderr, " bad pixel/weight mask not found \n");
            printf(" bad pixel/weight mask not found %s \n", badpixelname);
            exit(EXIT_FAILURE);
        } else {
            printf(" bad pixel/weight mask not in use, continuing \n");
            for(i = 0; i < bsize; i++)
                badpix[i] = 1.;
            status = 0;
        }
    }

/* sample pixels, fill up an array of values of pixels drawn from
the central region of the image, sort it
and hence find the median and 25-percentile.  If the image has been sky subtracted,
the difference divided by 0.67 will be a good estimate of the rms noise.
If the image has not been subtracted but has a uniform background, this will
still work.  If there are significant background variations, these will contribute
to the noise measurement, making it larger than just the random noise. 
repeat at a number of locations across the image and take the median of these noise
estimates.
*/

    medbin = 0;

    for(iy = 0; iy < 8; iy++) {
        ymin = dim[1] / 16 + (iy * dim[1]) / 8 - 50;
        ymax = ymin + 101;
        if (ymin < 0)
            ymin = 0;
        if (ymax >= dim[1])
            ymax = dim[1] - 1;

        for(ix = 0; ix < 4; ix++) {
            xmin = dim[0] / 8 + (ix * dim[0]) / 4 - 50;
            xmax = xmin + 101;
            if (xmin < 0)
                xmin = 0;
            if (xmax >= dim[0])
                xmax = dim[0] - 1;

            n = 0;
            for(y = ymin; y < ymax; y++) {
                for(x = xmin; x < xmax; x++) {
                    i = y * dim[0] + x;
                    if (badpix[i] > 0.) {
                        sortval[n] = apix[i];
                        n++;
                    }
                }
            }

            if (n > 0) {
                qsort(sortval, n, sizeof(float), compare);
                median = sortval[(n / 2)];
                quartile = sortval[(n / 4)];
                if (median > quartile) {
                    noiseval[medbin] = (median - quartile) / 0.67;
                    medbin++;
                }
            }
        }
    }

    if (medbin > 0) {
        qsort(noiseval, medbin, sizeof(float), compare);
        *mednoiseval = noiseval[(medbin / 2)];
    } else {
        *mednoiseval = 0.;
    }

    status = 0;
    return;

}


void dilatemask(float *badpix, int *dim)
{
  int x, y, i, j;
  
  /* dilate the bad pixel mask by one pixel in all directions - takes care
     of inaccuracies in the bad pixel mask */
  for (y=0; y<dim[1]; y++)
    {
      for (x=1; x<dim[0]; x++)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      badpix[i-1]=0.;
	    }
	}
      for (x=dim[0]-2; x>=0; x--)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      badpix[i+1]=0.;
	    }
	}
    }
  for (x=0; y<dim[0]; x++)
    {
      for (y=1; y<dim[1]; y++)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      j = x + (y-1)*dim[0];
	      badpix[j]=0.;
	    }
	}
      for (y=dim[1]-2; y>=0; y--)
	{
	  i = x + y*dim[0];
	  if (badpix[i]<=0.)
	    {
	      j = x + (y+1)*dim[0];
	      badpix[j]=0.;
	    }
	}
    }
}
