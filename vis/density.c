#include "../header/copyright.h"
/*=============================================================================
 * FILE: density.c
 *
 * PURPOSE: Contains functions to calculate the number densities of the chemical
 *   species and related quantities.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_numberden()  - calculate the initial number densities
 *   reset_numberden() - scale the number density to new densities
 *   denscale()        - calculate the number density variation scale
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
/* Calculate the initial number density
 * Rho is the mass density, in unit of g/cm^3
 */
void init_numberden(ChemEvln *Evln, Real rho, int verbose)
{

  Real sum = 0.0;
  Chemistry   *Chem = Evln->Chem;
  ElementInfo *Ele  = Chem->Elements;
  SpeciesInfo *Spe;
  int i, k;

  Evln->rho = rho;
  Evln->t   = 0.0;

/* Calculate abundance - number density ratio */

  sum = 0.0;
  for (i=0; i<Chem->N_Ele+Chem->NGrain; i++) {
    sum += Ele[i].mass * Ele[i].abundance;
  }

  Evln->Abn_Den = (sum*1.672e-24)/rho; /* == 1/n_H */

/* Calculate initial number density */
  ath_pout(verbose,"\n");
  ath_pout(verbose,"Initial number density:\n");
  ath_pout(verbose,"\n");

  for (i=0; i<Chem->Ntot; i++)
    Evln->NumDen[i] = 0.0;

  for (i = 0; i < Chem->N_Ele+Chem->NGrain; i++)
  {/* Initialize with the first single-element species of each element */

    k = Ele[i].single[0];
    Spe = &(Chem->Species[k]);

    Evln->NumDen[k] = Ele[i].abundance / (Evln->Abn_Den * Spe->composition[i]);
    ath_pout(verbose, "[%s]: = %e\n", Spe->name, Evln->NumDen[k]);
  }

/* Calculate the density variation scale */

  denscale(Evln, verbose);

  return;
}

/*----------------------------------------------------------------------------*/
/* Reset the number density with the new rho
 */
void reset_numberden(ChemEvln *Evln, Real rho_new, int verbose)
{
  int i;
  Real ratio;
  Chemistry   *Chem = Evln->Chem;

  ratio = rho_new/Evln->rho;

  Evln->rho *= ratio;

  Evln->Abn_Den /= ratio+TINY_NUMBER;

  for (i=0; i<Chem->Ntot; i++)
    Evln->NumDen[i] *= ratio;

  denscale(Evln, verbose);

  return;
}

/*----------------------------------------------------------------------------*/
/* user provided routine to calculate the maximum density variation scale,
 * which is set by the maximum possible density of the species
 * vb: verbose
 */
void denscale(ChemEvln *Evln, int vb)
{
  int i, j;
  Real den, denmin, denmaxt = 0.0;
  Chemistry *Chem = Evln->Chem;

  ath_pout(vb,"\n");
  ath_pout(vb,"Number density variation scales for each species:\n");
  ath_pout(vb,"\n");

  for (i=1; i<Chem->Ntot; i++) /* exclude electron */
  {
    denmin = 1.0e23;

    for (j=0; j<Chem->N_Ele + Chem->NGrain; j++)
    {/* maximum number density possile for this species */

      if (Chem->Species[i].composition[j] > 0)
      {
        den = Chem->Elements[j].abundance
              /(Evln->Abn_Den*Chem->Species[i].composition[j]);

        if (den < denmin)
          denmin = den;
      }
    }

    /* For normal species, set it to 1% of the maximum possible density */

    Evln->DenScale[i] = 0.01*denmin/(Real)(Chem->Ntot);

    if (denmin > denmaxt)  denmaxt = denmin;

    ath_pout(vb,"[%7s]: density scale is %e\n",
                       Chem->Species[i].name, Evln->DenScale[i]);
  }

  /* For e-, set it to relatively large value */
  Evln->DenScale[0] = 1.0e-7 * denmaxt;
  ath_pout(vb,"[     e-]: density scale is %e\n",Evln->DenScale[0]);
  ath_pout(vb,"\n");

  return;
}

#endif /* CHEMISTRY */
