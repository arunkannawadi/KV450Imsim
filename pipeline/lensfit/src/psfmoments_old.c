#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

void
psfmoments(double *psf, int pheight, int pwidth, double *psfe,
           double *centroid, double *moments)
{
    double Q1, Q2, Q11, Q12, Q22, Q_norm, denom;
    double dist, distsq, xxcen, yycen, xx, yy;
    double pval, sigma, sigmasq, rsq, w;
    int i, x, y, iter, niter;

    // initialise variables and set default return values
    Q1 = 0.;
    Q2 = 0.;
    Q11 = 0.;
    Q12 = 0.;
    Q22 = 0.;
    Q_norm = 0.;
    psfe[0] = 0.;
    psfe[1] = 0.;
    moments[0] = moments[1] = moments[2] = 0.;
    centroid[0] = centroid[1] = 0.;

    // sigma of Gaussian weight function for astrometry and moments
    sigma = 2.0; 
    sigmasq = sigma*sigma;
    niter = 3; // number of moments iterations for centroid
    // distance limit for moments summations
    dist = (double)pwidth/2; 
    distsq = dist*dist;

    // to match with the original lensfit tophat moments, set
    // niter=1, dist=8, sigma very large

    // iteratively find the weighted model centroid

    for (iter=0; iter<niter; iter++)
      {
	Q1 = Q2 = Q_norm = 0.; // initialise summations for this iteration
	for (y=0; y<pheight; y++)
	  {
	    yy = y>pheight/2 ? y-pheight : y; // allow for swapped quadrants
	    yycen = yy - centroid[1]; // centre the weight function on the current best guess centroid
	    for (x=0; x<pwidth; x++)
	      {
		i = y*pwidth + x;
		xx = x>pwidth/2 ? x-pwidth : x; // allow for swapped quadrants
		dist = xx*xx + yy*yy;
		if (dist < distsq) // force outer boundary to be circular centred on zero
		  {
		    xxcen = xx - centroid[0]; // centre the weight function on the current best guess centroid
		    rsq = xxcen*xxcen + yycen*yycen; // radius squared defined in oversampled units relative to best-guess centroid
		    w = exp(-rsq/2./sigmasq); // weight function
		    pval = psf[i];
		    Q1 += xx*pval*w;
		    Q2 += yy*pval*w;
		    Q_norm += pval*w;
		  }
	      }
	  }
	if (Q_norm>0.)
	  {
	    centroid[0] = Q1/Q_norm;
	    centroid[1] = Q2/Q_norm;
	  }
      }

    if (Q_norm > 0.) {
        /*
         * measure second moments 
         */
	Q_norm = 0.; // initialise summation
	for (y=0; y<pheight; y++)
	  {
	    yy = y>pheight/2 ? y-pheight : y; // allow for swapped quadrants
	    yycen = yy - centroid[1]; // centre the weight function on the centroid
	    for (x=0; x<pwidth; x++)
	      {
		i = y*pwidth + x;
		xx = x>pwidth/2 ? x-pwidth : x; // allow for swapped quadrants
		dist = xx*xx + yy*yy;
		if (dist < distsq) // force outer boundary to be circular centred on zero
		  {
		    xxcen = xx - centroid[0]; // centre the weight function on the centroid
		    rsq = xxcen*xxcen + yycen*yycen; // radius squared relative to best-guess centroid
		    w = exp(-rsq/2./sigmasq); // weight function
		    pval = psf[i];
                    Q11 += pval * xx*xx * w;
                    Q12 += pval * xx*yy * w;
                    Q22 += pval * yy*yy * w;
                    Q_norm += pval * w;
		  }
	      }
	  }

        denom = Q11 * Q22 - Q12 * Q12;
        if (denom > 0.) {
            denom = Q11 + Q22 + 2. * sqrt(denom);
        } else {
            denom = Q11 + Q22;
        }
        if (denom != 0.) {
            psfe[0] = (Q11 - Q22) / denom;
            psfe[1] = (2. * Q12) / denom;
        }

        if (psfe[0] > 1.) psfe[0] = 1.;
        if (psfe[1] > 1.) psfe[1] = 1.;
        if (psfe[0] < -1.) psfe[0] = -1.;
        if (psfe[1] < -1.) psfe[1] = -1.;

        if (Q_norm > 0.) {
            moments[0] = Q11 / Q_norm;
            moments[1] = Q22 / Q_norm;
            moments[2] = Q12 / Q_norm;
        }

    }

}
