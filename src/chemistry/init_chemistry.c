#include "../header/copyright.h"
/*=============================================================================
 * FILE: init_chemistry.c
 * PURPOSE: It wraps up all the functions to initiate a chemistry model. 
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_chemistry()  - initialize a chemistry model
 *   init_chemevln()   - initialize a chemical evolution calculation
 *   final_chemistry() - finalize the chemistry model
 *   final_chemevln()  - finalize a chemical evolution calculation
 *
 * History:
 *   Written by  Xuening Bai      Nov. 2010
==============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/chemproto.h"
#include "../header/prototypes.h"

#ifdef CHEMISTRY

/*----------------------------------------------------------------------------*/
/* Initiate the chemistry structure
 */
void init_chemistry(Chemistry *Chem)
{

  /* Read all elements, species and dust properties
   * Arrays of Elements, Species, GrSize, GrFrac, NumDen, DenScale are initiated
   */
  init_species(Chem);

  /* Real and construct all reaction equations
   * Arrays of Reactions and K are initiated
   */
  init_reactions(Chem);

  /* Construct chemical evolution equations 
   * The Equations array is initiated
   */
  init_equations(Chem);

  return;
}

/*----------------------------------------------------------------------------*/
/* Initiate the chemical evolution arrays
 */
void init_chemevln(Chemistry *Chem, ChemEvln *Evln)
{
  if (Chem != NULL) {

    Evln->t = 0.0;

    Evln->Chem = Chem;

    Evln->NumDen   = (Real*)calloc_1d_array(Chem->Ntot, sizeof(Real));
    if (Chem->NGrain > 0) {
      Evln->GrAvail  = (Real*)calloc_1d_array(Chem->Ntot,sizeof(Real));
    }
    Evln->DenScale = (Real*)calloc_1d_array(Chem->Ntot, sizeof(Real));
    Evln->K        = (Real*)calloc_1d_array(Chem->NReaction, sizeof(Real));
    Evln->rate_adj = (Real*)calloc_1d_array(Chem->NReaction, sizeof(Real));
  }
  else {
    ath_error("[init_chemevln]: The Chemistry model is NULL!\n");
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* Initiate the chemical evolution arrays by duplication
 */
void dup_chemevln(ChemEvln *Evln, ChemEvln *Evln_new)
{
  int i;
  Chemistry *Chem = Evln->Chem;

  if (Chem != NULL) {

    Evln_new->t    = Evln->t;

    Evln_new->Chem = Evln->Chem;

    Evln_new->NumDen   = (Real*)calloc_1d_array(Chem->Ntot, sizeof(Real));
    if (Chem->NGrain > 0) {
      Evln->GrAvail  = (Real*)calloc_1d_array(Chem->Ntot, sizeof(Real));
    }
    Evln_new->DenScale = (Real*)calloc_1d_array(Chem->Ntot, sizeof(Real));
    Evln_new->K        = (Real*)calloc_1d_array(Chem->NReaction, sizeof(Real));
    Evln_new->rate_adj = (Real*)calloc_1d_array(Chem->NReaction, sizeof(Real));

    for (i=0; i<Chem->Ntot; i++)
    {
      Evln_new->NumDen[i]   = Evln->NumDen[i];
      Evln_new->DenScale[i] = Evln->DenScale[i];
    }

    for (i=0; i<Chem->NReaction; i++)
    {
      Evln_new->K[i] = Evln->K[i];
      Evln_new->rate_adj[i] = Evln->rate_adj[i];
    }

    Evln_new->rho     = Evln->rho;
    Evln_new->T       = Evln->T;
    Evln_new->B       = Evln->B;
    Evln_new->Abn_Den = Evln->Abn_Den;
  }
  else {
    ath_error("[dup_chemevln]: The Chemistry model is NULL!\n");
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* Finalize the chemistry structure
 */
void final_chemistry(Chemistry *Chem)
{
  int i;

  for (i=0; i<Chem->Ntot; i++)
    free(Chem->Species[i].composition);

  for (i=0; i<Chem->NReaction; i++)
    free(Chem->Reactions[i].coeff);

  for (i=0; i<Chem->Ntot; i++)
    free(Chem->Equations[i].EqTerm);

  if (Chem->NGrain > 0) {
    free(Chem->GrSize);
    free(Chem->GrFrac);
  }
  free(Chem->Elements);
  free(Chem->Species);
  free(Chem->Reactions);
  free(Chem->Equations);

  return;
}

/*----------------------------------------------------------------------------*/
/* Finalize the chemical evolution arrays
 */
void final_chemevln(ChemEvln *Evln)
{
  Chemistry *Chem = Evln->Chem;
  if (Chem->NGrain > 0) {
    free(Evln->GrAvail);
  }

  Evln->Chem = NULL;

  free(Evln->K);
  free(Evln->rate_adj);
  free(Evln->NumDen);
  free(Evln->DenScale);

  return;
}

#endif /* CHEMISTRY */
