#include "../header/copyright.h"
/*=============================================================================
 * FILE: disk.c
 *
 * PURPOSE: Contains all functions related to (protoplanetary) disk models.
 *   We adpot a rescaled version of the minimum-mass solar nebular model, with
 *   all disk parameters having a power law dependence on disk radius. The
 *   X-ray and cosmic ray ionization rate is taken from literature.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  - init_disk()         - read disk parameters
 *  - Ionization_disk()   - calculate the ionization rate
 *  - SurfDen_disk()      - disk surface density
 *  - Temp_disk()         - disk temperature
 *  - Height_disk()       - disk scale height
 *  - Rho_disk()          - mass density in the disk
 *  - SurfDenZ_disk()     - disk column density from the top
 *
 * REFERENCES:
 *   Bai, X.-N. & Goodman, J., 2009, ApJ, 701, 737
 *
 * History:
 *   Written by  Xuening Bai      Nov. 2010
==============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/chemproto.h"
#include "../header/prototypes.h"

#ifdef CHEMISTRY

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   There is no private function.
 *============================================================================*/

/*============================================================================*/
/*---------------------------- Public Functions ------------------------------*/

/*----------------------------------------------------------------------------*/
/* Read the disk parameters
 */
void init_disk(Nebula *Disk)
{
  //Disk->Mstar   = par_getd_def("disk","Mstar",  1.0);
  Disk->Mstar   = par_getd("disk","Mstar");

  //Disk->Sigma0  = par_getd_def("disk","Sigma",  1700.0);
  Disk->Sigma0  = par_getd("disk","Sigma");
  //Disk->pS      = par_getd_def("disk","pS",     1.5);
  Disk->pS      = par_getd("disk","pS");

  //Disk->Tdisk   = par_getd_def("disk","Tdisk",  280.0);
  Disk->Tdisk   = par_getd("disk","Tdisk");
  //Disk->pT      = par_getd_def("disk","pT",     0.5);
  Disk->pT      = par_getd("disk","pT");

  //Disk->CR_rate = par_getd_def("disk","CR_rate",1.0e-17);
  Disk->CR_rate = par_getd("disk","CR_rate");
  //Disk->Lx      = par_getd_def("disk","Lx",     1.0e30);
  Disk->Lx      = par_getd("disk","Lx");
  //Disk->Tx      = par_getd_def("disk","Tx",     5.0);
  Disk->Tx      = par_getd("disk","Tx");

  //Disk->RD_rate = par_getd_def("disk","RD_rate",7.0e-19);
  Disk->RD_rate = par_getd("disk","RD_rate");

  return;
}

/*----------------------------------------------------------------------------*/
/* Obtain the ionization rate using a solar nebula model
 * Reference: Fromang et al. (2002), Igea & Glassgold (1999)
 *            fitting formula from Bai & Goodman (2009)
 */
Real Ionization_disk(Nebula *Disk, Real radius, Real z)
{
  Real sigma, sigmaz, nh1, nh2;
  Real rate;
  Real r1, col1, pow1, r2, col2, pow2, coef3, coef5, lnTx;

  /* fitting coefficients for Tx=3keV and Tx=5keV */
  Real r13,col13,pow13,r23,col23,pow23;
  Real r15,col15,pow15,r25,col25,pow25;

  r13 = 3.0e-12;    col13 = 1.5e21;    pow13 = 0.40;
  r23 = 0.5e-15;    col23 = 7.0e23;    pow23 = 0.65;
  r15 = 2.0e-12;    col15 = 3.0e21;    pow15 = 0.50;
  r25 = 1.5e-15;    col25 = 1.0e24;    pow25 = 0.70;

/*---------- Cosmic ray ionization rate ----------------
 */

  /* disk column density */
  sigma  = SurfDen_disk (Disk, radius);
  sigmaz = SurfDenZ_disk(Disk, radius, z);

  rate = Disk->CR_rate * (exp(-sigmaz/96.0) + exp(-(sigma-sigmaz)/96.0));

/*---------- X-ray ionization rate ---------------------
 */

  /* hydrogen column density (mean atomic weight = 1.425) */
  nh1 =        sigmaz  * 6.02e23 / 1.425;
  nh2 = (sigma-sigmaz) * 6.02e23 / 1.425;

  /* The following fitting formula applies for 1 < TX < 8 (keV) */

  if ((Disk->Tx < 1.0) || (Disk->Tx > 8.0))
    ath_error("The X-ray source temperature must be between 1 and 8 keV!\n");

  lnTx = log(Disk->Tx);

  coef3 = (log(5.0)-lnTx)/(log(5.0)-log(3.0));
  coef5 = (lnTx-log(3.0))/(log(5.0)-log(3.0));

  r1   = exp(coef3*log(r13)   + coef5*log(r15));
  col1 = exp(coef3*log(col13) + coef5*log(col15));
  pow1 = exp(coef3*log(pow13) + coef5*log(pow15));
  r2   = exp(coef3*log(r23)   + coef5*log(r25));
  col2 = exp(coef3*log(col23) + coef5*log(col25));
  pow2 = exp(coef3*log(pow23) + coef5*log(pow25));

  rate += 2.0*(Disk->Lx/1.0e29)*pow(radius,-2.2)
          * (r1*(exp(-pow(nh1/col1,pow1)) + exp(-pow(nh2/col1,pow1)))
           + r2*(exp(-pow(nh1/col2,pow2)) + exp(-pow(nh2/col2,pow2))));

/*---------- Radioactive decay ionization rate ---------------------
 */

  rate += Disk->RD_rate;

  return rate;
}


/*------------------------------------------------------------------------------
 * Disk surface density as a function of radius (AU)
 */
Real SurfDen_disk(Nebula *Disk, Real radius)
{/* g/cm^2 */

  return Disk->Sigma0 * pow(radius, -Disk->pS);
}

/*------------------------------------------------------------------------------
 * Disk temperature as a function of radius (AU)
 */
Real Temp_disk(Nebula *Disk, Real radius)
{/* Kelvin */

  return Disk->Tdisk * pow(radius, -Disk->pT);
}

/*------------------------------------------------------------------------------
 * Isothermal sound speed as a function of radius (AU)
 */
Real Cs_disk(Nebula *Disk, Real radius)
{/* cm/s */

  Real T;

  T  = Temp_disk(Disk, radius);

  return 9.088e3*sqrt(T/MUN);
}

/*------------------------------------------------------------------------------
 * Keperian frequency as a function of radius (AU)
 */
Real Omega_disk(Nebula *Disk, Real radius)
{/* s^(-1) */

  return 1.99e-7*Disk->Mstar*pow(radius, -1.5);
}

/*------------------------------------------------------------------------------
 * Disk scale height as a function of radius (AU)
 */
Real Height_disk(Nebula *Disk, Real radius)
{/* AU */

  return 0.034 * sqrt(Disk->Tdisk/280.0/Disk->Mstar)
               * pow (radius, 0.5*(3.0-Disk->pT));
}

/*------------------------------------------------------------------------------
 * Disk mass density as a function of radius (AU) and height (H)
 */
Real Rho_disk(Nebula *Disk, Real radius, Real z)
{/* g/cm^3 */

  Real h, sigma;

  sigma = SurfDen_disk(Disk, radius);
  h     = Height_disk (Disk, radius);

  return 6.68e-14*sigma/sqrt(2.0*PI)/h*exp(-0.5*SQR(z));
}

/*------------------------------------------------------------------------------
 * Reference diffusion coefficient as a function of radius (AU)
 * = c_s^2/Omega
 */
Real eta0_disk(Nebula *Disk, Real radius)
{/* cm^2/s */

  Real Cs, Omega;

  Cs    = Cs_disk (Disk, radius);
  Omega = Omega_disk(Disk, radius);

  return SQR(Cs)/Omega;
}

/*------------------------------------------------------------------------------
 * Surface mass density from top to height z (H) as a function of radius (AU)
 */
Real SurfDenZ_disk(Nebula *Disk, Real radius, Real z)
{/* g/cm^2 */

  Real h, sigma, z0;

  z0    = fabs(z);
  sigma = SurfDen_disk(Disk, radius);

  return 0.5 * sigma * (1.0 - Erf(z0/1.41421356));
}



/* ionization as a function of sigma and rho */
Real Ionization_disk1(Nebula *Disk, Real radius, Real sigmaz)
{
  Real sigma,nh1, nh2;
  Real rate;
  Real r1, col1, pow1, r2, col2, pow2, coef3, coef5, lnTx;

  /* fitting coefficients for Tx=3keV and Tx=5keV */
  Real r13,col13,pow13,r23,col23,pow23;
  Real r15,col15,pow15,r25,col25,pow25;

  r13 = 3.0e-12;    col13 = 1.5e21;    pow13 = 0.40;
  r23 = 0.5e-15;    col23 = 7.0e23;    pow23 = 0.65;
  r15 = 2.0e-12;    col15 = 3.0e21;    pow15 = 0.50;
  r25 = 1.5e-15;    col25 = 1.0e24;    pow25 = 0.70;

/*---------- Cosmic ray ionization rate ----------------
 */

  /* disk column density */
  sigma  = SurfDen_disk (Disk, radius);

  rate = Disk->CR_rate * (exp(-sigmaz/96.0) + exp(-(sigma-sigmaz)/96.0));

/*---------- X-ray ionization rate ---------------------
 */

  /* hydrogen column density (mean atomic weight = 1.425) */
  nh1 =        sigmaz  * 6.02e23 / 1.425;
  nh2 = (sigma-sigmaz) * 6.02e23 / 1.425;

  /* The following fitting formula applies for 1 < TX < 8 (keV) */

  if ((Disk->Tx < 1.0) || (Disk->Tx > 8.0))
    ath_error("The X-ray source temperature must be between 1 and 8 keV!\n");

  lnTx = log(Disk->Tx);

  coef3 = (log(5.0)-lnTx)/(log(5.0)-log(3.0));
  coef5 = (lnTx-log(3.0))/(log(5.0)-log(3.0));

  r1   = exp(coef3*log(r13)   + coef5*log(r15));
  col1 = exp(coef3*log(col13) + coef5*log(col15));
  pow1 = exp(coef3*log(pow13) + coef5*log(pow15));
  r2   = exp(coef3*log(r23)   + coef5*log(r25));
  col2 = exp(coef3*log(col23) + coef5*log(col25));
  pow2 = exp(coef3*log(pow23) + coef5*log(pow25));

  rate += 2.0*(Disk->Lx/1.0e29)*pow(radius,-2.2)
          * (r1*(exp(-pow(nh1/col1,pow1)) + exp(-pow(nh2/col1,pow1)))
           + r2*(exp(-pow(nh1/col2,pow2)) + exp(-pow(nh2/col2,pow2))));

/*---------- Radioactive decay ionization rate ---------------------
 */

  rate += Disk->RD_rate;

  return rate;
}


#endif /* CHEMISTRY */
