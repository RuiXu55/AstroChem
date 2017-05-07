#include "../header/copyright.h"
/*=============================================================================
 * FILE: difusivity.c
 *
 * PURPOSE: Contains functions to calculate the magnetic diffusivities and
 *   related quantities.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   Cal_NIMHD()       - calculate the magnetic diffusivities
 *   Cal_recomb()      - calculate the recombination time
 *
 * REFERENCES:
 *   Wardle, M., 2007, ApSS, 311, 35
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
/* Calculate the strength of non-ideal MHD terms
 *  - B   is the strength of the presumed magnetic field (Gauss)
 *  - rho is the gas density (g/cm^3)
 *  - T   is temperature (K)
 * We assume the mean molecular weight of the neutrals is 2.34 m_H.
 */
void Cal_NIMHD(ChemEvln *Evln)
{
  int i, j, k;
  Chemistry *Chem = Evln->Chem;

  Real pre, n15, mratio, mu, rate, beta;
  Real sig_O, sig_H, sig_P, sig_perp;

  SpeciesInfo *Spe;

  /* preliminaries */
  n15 = Evln->rho/(MUN * 1.672e-9); /*  n / 10^15 cm^(-3)    */
  pre = 14.4 / Evln->B;             /*  ec / B               */

  /* calculate the conductivities */

  sig_O = 0.0;
  sig_H = 0.0;
  sig_P = 0.0;

  for (i=0; i<Chem->Ntot; i++)
  {
    Spe = &(Chem->Species[i]);

    if (Spe->charge != 0)
    {
      if (i == 0)
      { /* electron */
        rate = 8.3e-9 * MAX(1.0,sqrt(Evln->T/100.0));
      }
      else if (i < Chem->GrInd)
      { /* ions */
        mu = MUN * Spe->mass/(MUN + Spe->mass);

        rate = 2.0e-9 * sqrt(1.0/mu);
      }
      else
      { /* grains */
        rate = MAX(1.3e-9*fabs(Spe->charge),
                   4.0e-3*SQR(Spe->gsize)*sqrt(Evln->T/100.0));
//                 1.6e-7*SQR(Spe->gsize)*sqrt(Evln->T/100.0));
      }

      mratio = (MUN+Spe->mass)/Spe->mass;
      beta = (9.59e-12 / (rate * n15)) * (Spe->charge * Evln->B / MUN) * mratio;

      sig_O += pre * Evln->NumDen[i] * Spe->charge * beta;
      sig_H += pre * Evln->NumDen[i] * Spe->charge / (1.0 + SQR(beta));
      sig_P += pre * Evln->NumDen[i] * Spe->charge * beta / (1.0 + SQR(beta));
    }
  }

//fprintf(stderr,"sigO=%e,sigH=%e,sigP=%e\n",sig_O,sig_H,sig_P);

  /* calculate the diffusivities */

  sig_perp = sqrt(SQR(sig_H) + SQR(sig_P));

  Evln->eta_O = 7.15e19 / sig_O;
  Evln->eta_H = 7.15e19 * sig_H / SQR(sig_perp);
  Evln->eta_A = 7.15e19 * sig_P / SQR(sig_perp) - Evln->eta_O;

  return;
}

/*----------------------------------------------------------------------------*/
/* Calculate recombination time based on the rate of change of magnetic
 * diffusivities. The recombination time is calculated by:
 *
 *   t_recomb = eta / (d_eta/d_t)|_t0
 *
 * Input arguments:
 *  - Bmin, Bmax = the min and max strength of the magnetic field (Gauss)
 *  - nB         = number of magnetic field bins in the calculation
 *  - Dt         = differential time (must be much smaller than t_{O,H,A})
 *
 * Output arguments:
 *  - t_{O,H,A}  = recombination time based on Ohmic, Hall & ambipolar
 *                 diffusivities.
 * Warning: t_O is simply a number, while t_H and t_A are arrays of size nB,
 *          and must be pre-allocated!
 */
void Cal_recomb(ChemEvln *Evln, Real Bmin, Real Bmax, int  nB,   Real Dt,
                                           Real *t_O, Real *t_H, Real *t_A)
{
  int i;
  Real dlnB, dttry;

  Chemistry *Chem = Evln->Chem;
  ChemEvln myEvln;

  /* duplicate the original evolution model */

  dup_chemevln(Evln, &myEvln);

  /* deplete the ionization rate */

  IonizationCoeff(&myEvln, 0.0, 0.0, 1);

  /* evolve the chemistry network for Dt */

  dttry = 0.1 * Dt;
  evolve( Dt, dttry, ChemErr);

  /* calculate the recombination time */

  dlnB = log(Bmax/Bmin)/nB;

  for (i=0; i<nB; i++)
  {
    Evln->B  = Bmin * exp((i+0.5)*dlnB);
    myEvln.B = Evln->B;

    Cal_NIMHD(Evln);
    Cal_NIMHD(&myEvln);

    *t_O   = Dt * fabs(Evln->eta_O/(myEvln.eta_O-Evln->eta_O));
    t_H[i] = Dt * fabs(Evln->eta_H/(myEvln.eta_H-Evln->eta_H));
    t_A[i] = Dt * fabs(Evln->eta_A/(myEvln.eta_A-Evln->eta_A));
  }

  final_chemevln(&myEvln);

  return;
}

#endif /* CHEMISTRY */

