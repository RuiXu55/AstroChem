#include "../copyright.h"
/*=============================================================================
 * FILE: init_reactions.c
 *
 * PURPOSE: Initialize all chemical reactions. It reads ionization reactions
 *   and other chemical reactions from standard input file "ChemReactions.txt"
 *   in the same directory as the executable. Below is a sample input format:
 *
 *  # Number of Ionization Reactions
 *    1
 *  # List of Ionization Reactions
 *  # R1    R2      P1      P2      P3      P4      ratio
 *    H2    0       H2+     e-      0       0       1.0
 *  # Number of Gas-phase Reactions
 *    3
 *  # List of Gas-phase Reactions
 *  # Type  R1    R2    P1    P2   P3   P4   alpha   beta  gamma  Tmin   Tmax
 *    RR    H2+   e-    H2    0    0    0    3.0E-6  -0.5   0     1      100000
 *    CE    H2+   Mg    Mg+   H2   0    0    3.0E-9  0      0     1      100000
 *    RR    Mg+   e-    Mg    0    0    0    3.0E-11 -0.5   0     1      100000
 *
 *  This code further construct all grain related reactions, termed into 
 *  different types. For details on reaction construction, see Bai & Goodman
 *  (2009).
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_reactions()
 *   InsertReactionInit()
 *   CheckReaction()
 *   FindInverse()
 *   PrintReaction()
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
#include "../defs.h"
#include "chemistry.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef CHEMISTRY

#define DNR 20 /* number of reactions added each time */

/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*----------------------------------------------------------------------------*/
/* Second step of initializing the chemical model
 * Read all chemical reactions, construct all grain related reactions
 */
void extract_reactions(Chemistry *Chem)
{
  FILE *fp;
  int i, k, n, ind, val;
  int speclab[10];
  char line[MAXLEN],name[NL_SPE];
  Real coef;

  fp = fopen("AllReactions.txt","r");

  if(fp == NULL){
    ath_perr(-1,"[extract_reactions]: Unable to open the AllReactions.txt!\n");
    return;
  }

  Chem->NReaction = 0;
  Chem->Reaction_size = DNR;
  Chem->Reactions = (ReactionInfo*)calloc_1d_array(DNR, sizeof(ReactionInfo));

/*----------- Read and construct other gas-phase chemical reactions ----------*/

  fscanf(fp, "%d\n", &k);	/* Number of chemical reactions; */
  if (k < 0)
    ath_error("[init_reactions]: Number of chemical reactions must be >=0!\n");

  for (i=0; i<k; i++)
  {
    n = Chem->NReaction;
    InsertReactionInit(Chem);
    val = 0;

    fscanf(fp, "%d", &ind);
    fscanf(fp, "%s", Chem->Reactions[n].rtype);

    fscanf(fp, "%s", name);
    Chem->Reactions[n].reactant[0] = FindSpecies(Chem, name); /* reactant 1 */
    val = MIN(val, Chem->Reactions[n].reactant[0]);

    fscanf(fp, "%s", name);
    Chem->Reactions[n].reactant[1] = FindSpecies(Chem, name); /* reactant 2 */
    //if (strcmp(name,"PHOTON") != 0)
    val = MIN(val, Chem->Reactions[n].reactant[1]);

    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[0] = FindSpecies(Chem, name);  /* product 1 */
    val = MIN(val, Chem->Reactions[n].product[0]);

    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[1] = FindSpecies(Chem, name);  /* product 2 */
    if (strcmp(name,"PHOTON") != 0)
      val = MIN(val, Chem->Reactions[n].product[1]);

    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[2] = FindSpecies(Chem, name);  /* product 3 */
    val = MIN(val, Chem->Reactions[n].product[2]);

    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[3] = FindSpecies(Chem, name);  /* product 4 */
    val = MIN(val, Chem->Reactions[n].product[3]);

    fscanf(fp,"%d",&ind);

    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].alpha = coef;      /* coefficient 1 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].beta = coef;       /* coefficient 2 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].gamma = coef;      /* coefficient 3 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].Tmin = coef;       /* coefficient 4: Tmin */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].Tmax = coef;       /* coefficient 5: Tmax */

    fgets(line,MAXLEN,fp); 

    if (val == -10) /* reaction contain unknown species */
    {
      Chem->NReaction -= 1;
      continue;
    }
    else if (CheckReaction(Chem, Chem->Reactions[n]) < 0)
    {
       PrintReaction(Chem, n, 0.0);
       ath_perr(-1,"[init_reactions]: Error in reaction %d!\n", n+1);
    }
  }

  fclose(fp);

/*----------------------- End of constructing reactions ----------------------*/

  /* Output All the reactions for checking purpose */
  n = Chem->NReaction;

  fp = fopen("ChemReactions.txt","w");
  fprintf(fp,"# UMIST 2012, no negative ions, no photoreactions\n");
  fprintf(fp,"# Number of Ionization Reactions\n");
  fprintf(fp,"4\n");
  fprintf(fp,"# List of Ionization Reactions\n");
  fprintf(fp,"# R1    R2      P1      P2      P3      P4      ratio\n");
  fprintf(fp,"H2    0       H2+     e-      0       0       0.97\n");
  fprintf(fp,"H2    0       H+      H       e-      0       0.03\n");
  fprintf(fp,"H     0       H+      e-      0       0       0.50\n");
  fprintf(fp,"He    0       He+     e-      0       0       0.84\n");
  fprintf(fp,"# Number of Gas-phase Reactions\n");

  fprintf(fp,"%d\n", n);
  fprintf(fp,"# List of Gas-phase Reactions\n");
  fprintf(fp,"# type	R1	R2	P1	P2	P3	P4	alpha		beta		gamma		Tmin		Tmax\n");

  for (i=0; i<n; i++)
  {
    fprintf(fp,"%s	",Chem->Reactions[i].rtype);

    k = Chem->Reactions[i].reactant[0];
    fprintf(fp, "%s	", Chem->Species[k].name);
    Chem->Species[k].nreact += 1;

    k = Chem->Reactions[i].reactant[1];
    if (k<0) fprintf(fp, "0	");
    else {
      fprintf(fp, "%s	", Chem->Species[k].name);
      Chem->Species[k].nreact += 1;
    }

    k = Chem->Reactions[i].product[0];
    if (k<0) fprintf(fp, "0	");
    else {
      fprintf(fp, "%s	", Chem->Species[k].name);
      Chem->Species[k].nprod += 1;
    }

    k = Chem->Reactions[i].product[1];
    if (k<0) fprintf(fp, "0	");
    else {
      fprintf(fp, "%s	", Chem->Species[k].name);
      Chem->Species[k].nprod += 1;
    }

    k = Chem->Reactions[i].product[2];
    if (k<0) fprintf(fp, "0	");
    else {
      fprintf(fp, "%s	", Chem->Species[k].name);
      Chem->Species[k].nprod == 1;
    }

    k = Chem->Reactions[i].product[3];
    if (k<0) fprintf(fp, "0	");
    else {
      fprintf(fp, "%s	", Chem->Species[k].name);
      Chem->Species[k].nprod += 1;
    }

    fprintf(fp,"%e	",Chem->Reactions[i].coeff[0].alpha);
    fprintf(fp,"%e	",Chem->Reactions[i].coeff[0].beta);
    fprintf(fp,"%e	",Chem->Reactions[i].coeff[0].gamma);
    fprintf(fp,"%e	",Chem->Reactions[i].coeff[0].Tmin);
    fprintf(fp,"%e\n",    Chem->Reactions[i].coeff[0].Tmax);
  }

  fprintf(fp,"# Number of Grain Surface Reactions\n");
  fprintf(fp,"0\n");
  fprintf(fp,"# List of Grain Surface Reactions\n");
  fprintf(fp,"# R1    R2      P1      P2      Ea\n");
  fprintf(fp," H     H       H2      0       0.0\n");
  fclose(fp);


  for (k=0; k<Chem->Ntot; k++) {
    ath_pout(0,"k=%3d, [%12s]: nreact = %3d, nproduct = %3d\n", k,
             Chem->Species[k].name,Chem->Species[k].nreact,Chem->Species[k].nprod);
  }

  return;
}


/*----------------------------------------------------------------------------*/
/* Initial step before inserting a reaction
 */
void InsertReactionInit(Chemistry *Chem)
{
  int n;
  ReactionInfo *React;

  Chem->NReaction += 1;
  n = Chem->NReaction-1;
  if (n >= Chem->Reaction_size)
  {
    Chem->Reactions = (ReactionInfo*)realloc(Chem->Reactions,
                      (Chem->Reaction_size+DNR)*sizeof(ReactionInfo));
    Chem->Reaction_size += DNR;
  }

  React = &(Chem->Reactions[n]);

  React->coeff = (Coefficient*)realloc(React->coeff, sizeof(Coefficient));
  React->reactant[0] = -1;
  React->reactant[1] = -1;
  React->product[0] = -1;
  React->product[1] = -1;
  React->product[2] = -1;
  React->product[3] = -1;
  React->coeff[0].alpha = 0.0;
  React->coeff[0].beta = 0.0;
  React->coeff[0].gamma = 0.0;
  React->coeff[0].Tmin = 0.0;
  React->coeff[0].Tmax = 0.0;

  return;
}


/*---------------------------------------------------------------------------*/
/* Checek whether a reaction satisfy conservation laws
 */
int CheckReaction(Chemistry *Chem, ReactionInfo Reaction)
{
  int i, r1, r2, p1, p2, p3, p4;
  int Comp;

  r1 = Reaction.reactant[0];    r2 = Reaction.reactant[1];
  p1 = Reaction.product[0];     p2 = Reaction.product[1];
  p3 = Reaction.product[2];     p4 = Reaction.product[3];

  for (i=0; i<Chem->N_Ele; i++)
  {
    Comp = 0;
    if (r1>=0) Comp += Chem->Species[r1].composition[i];
    if (r2>=0) Comp += Chem->Species[r2].composition[i];
    if (p1>=0) Comp -= Chem->Species[p1].composition[i];
    if (p2>=0) Comp -= Chem->Species[p2].composition[i];
    if (p3>=0) Comp -= Chem->Species[p3].composition[i];
    if (p4>=0) Comp -= Chem->Species[p4].composition[i];
    if (Comp != 0) return -1;   /* This reaction does not conserve element! */
  }

  if (r1>=0) Comp += Chem->Species[r1].charge;
  if (r2>=0) Comp += Chem->Species[r2].charge;
  if (p1>=0) Comp -= Chem->Species[p1].charge;
  if (p2>=0) Comp -= Chem->Species[p2].charge;
  if (p3>=0) Comp -= Chem->Species[p3].charge;
  if (p4>=0) Comp -= Chem->Species[p4].charge;

  if (Comp != 0) return -1;     /* This reaction does not conserve charge! */

  return 0;
}

/*---------------------------------------------------------------------------*/
/* Print Reaction i
 */
void PrintReaction(Chemistry *Chem, int i, Real K)
{
  int k;

  ath_pout(0, "Reaction %4d (K=%10e): ", i+1, K);

  k = Chem->Reactions[i].reactant[0];
  ath_pout(0, "%7s ", Chem->Species[k].name);

  k = Chem->Reactions[i].reactant[1];
  if (k<0) ath_pout(0, "          -> ");
  else ath_pout(0, "+ %7s -> ", Chem->Species[k].name);

  k = Chem->Reactions[i].product[0];
  ath_pout(0, "%7s ", Chem->Species[k].name);

  k = Chem->Reactions[i].product[1];
  if (k>=0) ath_pout(0, "+ %7s ", Chem->Species[k].name);

  k = Chem->Reactions[i].product[2];
  if (k>=0) ath_pout(0, "+ %7s ", Chem->Species[k].name);

  k = Chem->Reactions[i].product[3];
  if (k>=0) ath_pout(0, "+ %7s ", Chem->Species[k].name);

  ath_pout(0,"\n");
  return;
}


#undef DNR

#endif /* CHEMISTRY */

