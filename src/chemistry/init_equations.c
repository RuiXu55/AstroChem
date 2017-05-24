#include "../header/copyright.h"
/*=============================================================================
 * FILE: init_equations.c
 * PURPOSE: Construct chemical evolution equations from the entire chemical
 *   reaction network. All equation information are stored in Chem->Equations.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_equations()
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
#include <string.h>
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/chemproto.h"
#include "../header/prototypes.h"

#ifdef CHEMISTRY

#define DNE 10

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   InsertTerm()    - insert a reaction term into the evolution equation
 *============================================================================*/
void InsertTerm(Chemistry *Chem, int a, int nreactant, int k, int sign);
void OutputEquation(Chemistry *Chem);

/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*----------------------------------------------------------------------------*/
/* The Third  step of initializing the chemical model
 * Construct chemical evolution equations from all reactions
 */
void init_equations(Chemistry *Chem)
{
  int i, k, l;
  int a, b, c, d, e, f;
  int nreactant;

/* Memory allocation */
  Chem->Equations = (EquationInfo*)calloc_1d_array(Chem->Ntot,
                                                    sizeof(EquationInfo));
  for (i=0; i<Chem->Ntot; i++)
  {
    Chem->Equations[i].NTerm = 0;
    Chem->Equations[i].EqTerm_size = DNE;
    Chem->Equations[i].EqTerm = (EquationTerm*)calloc_1d_array(DNE,
                                                     sizeof(EquationTerm));
  }

/* Construct equation */
  for (i=0; i<Chem->Ntot; i++)
    Chem->Equations[i].NTerm = 0;

  l = 0;
  for (i=0; i<Chem->NReaction; i++)
  {
    nreactant = 2;
    a = Chem->Reactions[i].reactant[0];
    b = Chem->Reactions[i].reactant[1];
    c = Chem->Reactions[i].product[0];
    d = Chem->Reactions[i].product[1];
    e = Chem->Reactions[i].product[2];
    f = Chem->Reactions[i].product[3];
    if (Chem->Reactions[i].use == 1)
    {
      if ((b >= 0) || (b == -10))
      {
        if (b == -10)
          nreactant = 1;  /* no second reactant */
        InsertTerm(Chem, a, nreactant, i, -1);
        InsertTerm(Chem, c, nreactant, i, 1);
        if (b >= 0) InsertTerm(Chem, b, nreactant, i, -1);
        if (d >= 0) InsertTerm(Chem, d, nreactant, i, 1);
        if (e >= 0) InsertTerm(Chem, e, nreactant, i, 1);
        if (f >= 0) InsertTerm(Chem, f, nreactant, i, 1);
      }
      /* neutral+grain and desorption reactions, contain all kinds of grains */
      else
      {
        InsertTerm(Chem, a, 1, i, -1);
        InsertTerm(Chem, c, 1, i, 1);
      }
      l += 1;
    }
  }

  return;
}

/*============================================================================*/
/*----------------------------- Private Functions ----------------------------*/

/*----------------------------------------------------------------------------*/
/* Inserting a reaction term in species a from reaction k
 * sign = 1 (a is a product) or -1 (a is a reactant)
 */
void InsertTerm(Chemistry *Chem, int a, int nreactant, int k, int sign)
{
  int i, m;
  EquationInfo *Eq;

  Eq = &(Chem->Equations[a]);

  m = Eq->NTerm;
  Eq->NTerm += 1;

  if (Eq->NTerm > Eq->EqTerm_size)
  {/* reallocate if there are more reaction terms */
    Eq->EqTerm = (EquationTerm*)realloc(Eq->EqTerm, 
                                 (Eq->EqTerm_size+DNE)*sizeof(EquationTerm));
    Eq->EqTerm_size += DNE;
  }

  /* Insert equation term */
  Eq->EqTerm[m].N = nreactant;
  for (i=0; i<nreactant; i++)
    Eq->EqTerm[m].lab[i] = Chem->Reactions[k].reactant[i];

  Eq->EqTerm[m].ind = k;
  Eq->EqTerm[m].dir = sign;
  Eq->EqTerm[m].type = Chem->Reactions[k].rtype;

  return;
}

/*----------------------------------------------------------------------------*/
/* Output all equations
 */
void OutputEquation(Chemistry *Chem)
{
  FILE *AllE;
  EquationTerm *Tm;
  int i,k,l;
  AllE = fopen("AllEquations.txt","w");

  for (i=0; i<Chem->Ntot; i++)
  {
    fprintf(AllE,"d[%s]/dt = ",Chem->Species[i].name);
    for (k=0; k<Chem->Equations[i].NTerm; k++)
    {
      Tm = &(Chem->Equations[i].EqTerm[k]);

      if (Tm->dir > 0)
        fprintf(AllE," + ");
      else
        fprintf(AllE," - ");

      fprintf(AllE,"K_%d",Tm->ind+1);

      for (l=0; l<Tm->N; l++)
        fprintf(AllE,"[%s]",Chem->Species[Tm->lab[l]].name);
    }
    fprintf(AllE,"\n\n");
  }
  fclose(AllE);
}

#undef DNE

#endif /* CHEMISTRY */
