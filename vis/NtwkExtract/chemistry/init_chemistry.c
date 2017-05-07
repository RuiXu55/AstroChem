#include "../copyright.h"
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
#include "chemistry.h"
#include "../prototypes.h"
#include "prototypes.h"

#ifdef CHEMISTRY

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

  free(Chem->Elements);
  free(Chem->Species);
  free(Chem->Reactions);

  return;
}

#endif /* CHEMISTRY */
