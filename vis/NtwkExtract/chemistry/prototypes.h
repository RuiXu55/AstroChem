#ifndef CHEMISTRY_PROTOTYPES_H
#define CHEMISTRY_PROTOTYPES_H 
//#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the src/chemistry
 *   directory.
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../defs.h"

/*----------------------------------------------------------------------------*/
#ifdef CHEMISTRY

/*----------------------------------------------------------------------------*/
/* init_chemistry.c */
void final_chemistry(Chemistry *Chem);

/*----------------------------------------------------------------------------*/
/* init_reactions.c */
void init_reactions(Chemistry *Chem);

void InsertReactionInit(Chemistry *Chem);
int  CheckReaction(Chemistry *Chem, ReactionInfo Reaction);
void PrintReaction(Chemistry *Chem, int i, Real K);

/*----------------------------------------------------------------------------*/
/* init_species.c */
void init_species(Chemistry *Chem);
int  FindSpecies(Chemistry *Chem, char name[NL_SPE]);

/*----------------------------------------------------------------------------*/
/* reprocess_list.c */
void reprocess_list(char *list);

#endif /* CHEMISTRY */

#endif /* CHEMISTRY_PROTOTYPES_H */
