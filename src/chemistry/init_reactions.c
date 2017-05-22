#include "../header/copyright.h"
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
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/chemproto.h"
#include "../header/prototypes.h"

#ifdef CHEMISTRY

#define DNR 20 /* number of reactions added each time */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   AddReactionTRange()   - add another temperature range
 *   CheckInverse()        - check if two reactions are inverse to each other
 *   ReactionCmp()         - compare current reaction with the previous one
 *============================================================================*/
void AddReactionTRange(Chemistry *Chem, int n);
int  CheckInverse(Chemistry *Chem, int i, int j);
int  ReactionCmp(Chemistry *Chem, int n);


/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*----------------------------------------------------------------------------*/
/* Second step of initializing the chemical model
 * Read all chemical reactions, construct all grain related reactions
 */
void init_reactions(Chemistry *Chem)
{
  FILE *fp;
  int h, i, j, k, l, m, n, sgn;
  int Ns_Gr = 2*Chem->GrCharge+1;
  int speclab[10];
  char line[MAXLEN],name[NL_SPE],fname[20],mantle[NL_SPE];
  char name1[NL_SPE],name2[NL_SPE],name3[NL_SPE],name4[NL_SPE];
  Real coef;

  sprintf(fname,"%s",par_gets("job","read_reaction"));
  fp = fopen(fname,"r");

  if(fp == NULL){
    ath_perr(-1,"[init_reactions]: Unable to open the ChemReactions.txt!\n");
    return;
  }

  Chem->NReaction = 0;
  Chem->Reaction_size = DNR;
  Chem->Reactions = (ReactionInfo*)calloc_1d_array(DNR, sizeof(ReactionInfo));

/*------------------ Read and construct ionization reactions -----------------*/

  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);
  fscanf(fp, "%d\n", &k);	/* Number of ionization reactions */
  if (k < 0)
    ath_error("[init_reactions]: Number of ionization reactions must >=0!\n");

  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);

  for (i=0; i<k; i++)
  {
    n = Chem->NReaction;
    InsertReactionInit(Chem);

    Chem->Reactions[n].rtype = 0;   /* Type is Ionization */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].reactant[0] = FindSpecies(Chem, name); /* reactant 1 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].reactant[1] = FindSpecies(Chem, name); /* reactant 2 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[0] = FindSpecies(Chem, name);  /* product 1 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[1] = FindSpecies(Chem, name);  /* product 2 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[2] = FindSpecies(Chem, name);  /* product 3 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[3] = FindSpecies(Chem, name);  /* product 4 */
    fscanf(fp, "%lf\n", &coef);
    Chem->Reactions[n].coeff[0].gamma = coef; /* coefficient 3: probability */
    Chem->Reactions[n].use = 1;

    if (CheckReaction(Chem, Chem->Reactions[n]) < 0)
    {
       PrintReaction(Chem,n,0.0);
       ath_perr(-1,"[init_reactions]: Error in reaction %d!\n", n+1);
    }

  }

/*----------- Read and construct other gas-phase chemical reactions ----------*/

  fgets(line,MAXLEN,fp);
  fscanf(fp, "%d\n", &k);	/* Number of chemical reactions; */
  if (k < 0)
    ath_error("[init_reactions]: Number of chemical reactions must be >=0!\n");

  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);

  for (i=0; i<k; i++)
  {
    n = Chem->NReaction;
    InsertReactionInit(Chem);

    Chem->Reactions[n].rtype = 1;   /* Type is gas-phase reaction */
    fscanf(fp, "%s", name);         /* sub-type */

    if (strcmp(name,"PH") == 0)     /* photo-reaction */
      Chem->Reactions[n].rtype = 10;

    fscanf(fp, "%s", name);
    Chem->Reactions[n].reactant[0] = FindSpecies(Chem, name); /* reactant 1 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].reactant[1] = FindSpecies(Chem, name); /* reactant 2 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[0] = FindSpecies(Chem, name);  /* product 1 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[1] = FindSpecies(Chem, name);  /* product 2 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[2] = FindSpecies(Chem, name);  /* product 3 */
    fscanf(fp, "%s", name);
    Chem->Reactions[n].product[3] = FindSpecies(Chem, name);  /* product 4 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].alpha = coef;      /* coefficient 1 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].beta = coef;       /* coefficient 2 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].gamma = coef;      /* coefficient 3 */
    fscanf(fp, "%lf", &coef);
    Chem->Reactions[n].coeff[0].Tmin = coef;       /* coefficient 4: Tmin */
    fscanf(fp, "%lf\n", &coef);
    Chem->Reactions[n].coeff[0].Tmax = coef;       /* coefficient 5: Tmax */
    Chem->Reactions[n].use = 1;

    if (CheckReaction(Chem, Chem->Reactions[n]) < 0)
    {
       PrintReaction(Chem, n, 0.0);
       ath_perr(-1,"[init_reactions]: Error in reaction %d!\n", n+1);
    }

    /* If this reaction is the same as last one, but temperature ranges are 
     * different, we should combine them */
    if (ReactionCmp(Chem, n) == 0)
    {
      n -= 1;
      m = Chem->Reactions[n].NumTRange;
      AddReactionTRange(Chem, n);

      /* Sort by temperature (Insert sorting method) */
      j = 0;
      while (j < m)
        if (  Chem->Reactions[n+1].coeff[0].Tmin 
            > Chem->Reactions[ n ].coeff[j].Tmin )   j++;
        else break;

      for (l=m-1; l>=j; l--)
        Chem->Reactions[n].coeff[l+1] = Chem->Reactions[n].coeff[l];

      Chem->Reactions[n].coeff[j] = Chem->Reactions[n+1].coeff[0];
      Chem->NReaction -= 1;
    }
  }

/*---------------------- Construct grain-phase reactions ---------------------*/

  /* e- + gr(n+): 2 * Chem->GrCharge reactions */

  for (k=0; k<Chem->NGrain; k++)
  {
    m = Chem->GrInd + k * Ns_Gr;               /* Ns_Gr = 2*GrCharge+1 */
    for (i=1; i<=2*Chem->GrCharge; i++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 2;            /* Type is ion-grain reaction */
      Chem->Reactions[n].reactant[0] = 0;            /* reactant 1: electron */
      Chem->Reactions[n].reactant[1] = m + i;        /* reactant 2: grain */
      Chem->Reactions[n].product[0] = m + i - 1;     /* product 1: grain */
                                                     /* product 2-4: none */
      /* alpha is set to electron mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[0].mass;
      Chem->Reactions[n].coeff[0].gamma = 1.0;       /* probability = 1 */
      Chem->Reactions[n].use = 1;
    }
  }

  /* ion+/- + gr(n+/-): reactions */

  for (j=Chem->N_Neu_f+1; j<=3*Chem->N_Neu_f; j++)
  {
    sgn = Chem->Species[j].charge;
    for (k=0; k<Chem->NGrain; k++)
    {
      m = Chem->GrInd + Chem->GrCharge + k * Ns_Gr;
      // l==1  charge ++ or --, l==2 +- or -+
      for (l=1; l<=2; l++)
        for (i=1; i<=Chem->GrCharge; i++)
        {
          n = Chem->NReaction;
          InsertReactionInit(Chem);
          Chem->Reactions[n].rtype = 2;      /* Type is ion-grain reaction */
          Chem->Reactions[n].reactant[0] = j;      /* reactant 1: ion */
                    
          if ((sgn==1) && (l==1))   // X+ + gr(n+) -> X[m] + gr(n+1)+, n=0,1
          {
           Chem->Reactions[n].reactant[1] = m + i - 1;  // reactant 2: grain 
           Chem->Reactions[n].product[0] =
                     j + k*(Chem->N_Neu_f+Chem->N_Neu+ 
										 Chem->N_Neu_s)+ Chem->ManInd;      // product 1: mantle species 
           Chem->Reactions[n].product[1] = m + i;       // product 2: grain 
           Chem->Reactions[n].use = 1;
          }

          if ((sgn==-1) && (l==1))  /* X- + gr(n-) -> X[m] + gr(n+1)-, n=0,1 */
          {
            Chem->Reactions[n].reactant[1] = m - i + 1;  /* reactant 2: grain */
            Chem->Reactions[n].product[0] =
								j + k*(Chem->N_Neu_f + Chem->N_Neu + 
								Chem->N_Neu_s)+Chem->ManInd;            /* product 1: mantle species */
            Chem->Reactions[n].product[1] = m - i;       /* product 2: grain */
            Chem->Reactions[n].use = 1;
          }

          if ((sgn==1) && (l==2))   /* X+ + gr(n-) -> X + gr(n-1)-, n=1,2 */
          {
            Chem->Reactions[n].reactant[1] = m - i;      /* reactant 2: grain */
            Chem->Reactions[n].product[0] =
                         j - Chem->N_Neu_f;   /* product 1: neutral counterpart */
            Chem->Reactions[n].product[1] = m - i + 1;   /* product 2: grain */
            Chem->Reactions[n].use = 1;
          }

          if ((sgn==-1) && (l==2))   /* X- + gr(n+) -> X + gr(n-1)+, n=1,2 */
          {
            Chem->Reactions[n].reactant[1] = m + i;      /* reactant 2: grain */
            Chem->Reactions[n].product[0] =
                         j - 2*Chem->N_Neu_f; /* product 1: neutral counterpart */
            Chem->Reactions[n].product[1] = m + i - 1;   /* product 2: grain */
            Chem->Reactions[n].use = 1;
          }

          /* alpha is set to ion mass! */
          Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
          /* beta is set to binding energy */
          Chem->Reactions[n].coeff[0].beta  = Chem->Species[j].Eb;
          Chem->Reactions[n].coeff[0].gamma = 1.0;       /* probability is 1 */

        }
      }
    }

  /* Ion + gr(n+): Chem->N_Neu * 3 * Chem->GrCharge reactions */

  for (j=Chem->NeuInd+Chem->N_Neu; j<Chem->NeuInd+2*Chem->N_Neu; j++)
    for (k=0; k<Chem->NGrain; k++)
    {
      m = Chem->GrInd + Chem->GrCharge + k * Ns_Gr;

      for (l=1; l<=2; l++)
        for (i=1; i<=Chem->GrCharge; i++)
        {
          n = Chem->NReaction;
          InsertReactionInit(Chem);
          Chem->Reactions[n].rtype = 2;      /* Type is ion-grain reaction */
          Chem->Reactions[n].reactant[0] = j;      /* reactant 1: ion */

          if (l == 1)	/* X+ + gr(n+) -> X[m] + gr(n+1)+, n=0,1 */
          {
           Chem->Reactions[n].reactant[1] = m + i - 1;  /* reactant 2: grain */
           Chem->Reactions[n].product[0] =
						 (j-Chem->NeuInd-Chem->N_Neu)+ 
						 Chem->ManInd + Chem->N_Neu_f;     // product 1: mantle species 
             //j + (k+1)*Chem->N_Neu;  /* product 1: mantle species */
           Chem->Reactions[n].product[1] = m + i;       /* product 2: grain */
           Chem->Reactions[n].use = 1;
          }

          if (l == 2)	/* X+ + gr(n-) -> X + gr(n-1)-, n=1,2 */
          {
            Chem->Reactions[n].reactant[1] = m - i;      /* reactant 2: grain */
            Chem->Reactions[n].product[0] =
                         j - Chem->N_Neu;   /* product 1: neutral counterpart */
            Chem->Reactions[n].product[1] = m - i + 1;   /* product 2: grain */
            Chem->Reactions[n].use = 1;
          }                                              /* product 3-4: none */
          /* alpha is set to ion mass! */
          Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
          /* beta is set to binding energy */
          Chem->Reactions[n].coeff[0].beta  = Chem->Species[j].Eb;
          Chem->Reactions[n].coeff[0].gamma = 1.0;       /* probability is 1 */

        }
    }

  /* Ion with no neutral counterpart + gr(n-): NumSpeCharge * Chem->GrCharge
   * reactions */

  for (j=Chem->SIonInd; j<Chem->SIonInd+Chem->N_Ion_s; j++)
  {/* X+ + gr(m-) -> gr((m-1)-) + Y + Z + ... */

    /* Search for dissociative recombination reactions for this ion */
    h = 0;	coef = 0.0;
    for (k=0; k<Chem->NReaction; k++)
    {
      if (Chem->Reactions[k].rtype == 1)
      if (  (Chem->Reactions[k].reactant[0] == j)
         && (Chem->Reactions[k].reactant[1] == 0))
      {
        speclab[h] = k;		h += 1;
        coef += Chem->Reactions[k].coeff[0].alpha;   /* Rate of reaction */
      }
    }

    for (k=0; k<Chem->NGrain; k++)
    {
      m = Chem->GrInd + Chem->GrCharge + k * Ns_Gr;
      for (i=1; i<=Chem->GrCharge; i++)
      for (l=0; l<h; l++)
      {
        n = Chem->NReaction;
        InsertReactionInit(Chem);
        Chem->Reactions[n].rtype = 2;    /* Type is ion-grain reaction */
        Chem->Reactions[n].reactant[0] = j;            /* reactant 1: ion */
        Chem->Reactions[n].reactant[1] = m - i;        /* reactant 2: grain- */
        Chem->Reactions[n].product[0] = m - i + 1;     /* product 1: grain-+1 */
        Chem->Reactions[n].product[1] =
               Chem->Reactions[speclab[l]].product[0]; /* product 2 */
        Chem->Reactions[n].product[2] =
               Chem->Reactions[speclab[l]].product[1]; /* product 3 */
        Chem->Reactions[n].product[3] =
               Chem->Reactions[speclab[l]].product[2]; /* product 4 */
        Chem->Reactions[n].use = 1;
        /* alpha is set to ion mass! */
        Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
        /* probability (branching ratio) of this reaction */
        Chem->Reactions[n].coeff[0].gamma =
               Chem->Reactions[speclab[l]].coeff[0].alpha/coef;
      }
    }
  }

  /* Neutral + gr(n+): Chem->N_Neu + Chem->N_Neu_s reactions */

  for (j=1; j<=Chem->N_Neu_f; j++)
    for (k=0; k<Chem->NGrain; k++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 3;    /* Type is neutral-grain reaction */
      Chem->Reactions[n].reactant[0] = j;      /* reactant 1: neutral */
      Chem->Reactions[n].reactant[1] = -k-1;   /* reactant 2: ALL GRAINS(k+1) */
      Chem->Reactions[n].product[0] =
                     j + k*(Chem->N_Neu_f+Chem->N_Neu+ 
										 Chem->N_Neu_s)+ Chem->ManInd;      // product 1: mantle species 
      Chem->Reactions[n].product[1] = -k-1;    /* product 2: ALL GRAINS(k+1) */
                                               /* product 3-4: none */
      /* alpha is neutral mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
      /* beta is binding energy */
      Chem->Reactions[n].coeff[0].beta = Chem->Species[j].Eb;
      Chem->Reactions[n].use = 1;
    }

  for (j=Chem->NeuInd; j<Chem->NeuInd+Chem->N_Neu; j++)
    for (k=0; k<Chem->NGrain; k++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 3;    /* Type is neutral-grain reaction */
      Chem->Reactions[n].reactant[0] = j;      /* reactant 1: neutral */
      Chem->Reactions[n].reactant[1] = -k-1;   /* reactant 2: ALL GRAINS(k+1) */
      Chem->Reactions[n].product[0] =
			  j-Chem->NeuInd+ Chem->ManInd + Chem->N_Neu_f;     // product 1: mantle species 
      Chem->Reactions[n].product[1] = -k-1;    /* product 2: ALL GRAINS(k+1) */
                                               /* product 3-4: none */
      /* alpha is neutral mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
      /* beta is binding energy */
      Chem->Reactions[n].coeff[0].beta = Chem->Species[j].Eb;
      Chem->Reactions[n].use = 1;
    }

  /* Adsorption */

  for (j=Chem->SNeuInd; j<Chem->SNeuInd+Chem->N_Neu_s; j++)
    for (k=0; k<Chem->NGrain; k++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 3;   /* Type is neutral-grain reaction */
      Chem->Reactions[n].reactant[0] = j;      /* reactant 1: neutral */
      Chem->Reactions[n].reactant[1] = -k-1;   /* reactant 2: ALL GRAINS(k+1) */
      Chem->Reactions[n].product[0] =
				j-Chem->SNeuInd + Chem->ManInd + 
				Chem->N_Neu_f + Chem->N_Neu;           /* product 1: mantle species */
      Chem->Reactions[n].product[1] = -k-1;    /* product 2: ALL GRAINS(k+1) */
                                               /* product 3-4: none */
      /* alpha is neutral mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
      /* beta is binding energy */
      Chem->Reactions[n].coeff[0].beta = Chem->Species[j].Eb;
      Chem->Reactions[n].use = 1;
    }

  /* Desorption: Chem->N_Neu+NumSpecNeutral reactions */
  for(j=Chem->ManInd+1;j<=Chem->ManInd+Chem->N_Neu_f;j++)
	{
    for (k=0; k<Chem->NGrain; k++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 4;           /* Type is desorption */
      Chem->Reactions[n].reactant[0] = 
        j + k*(Chem->N_Neu_f+Chem->N_Neu+Chem->N_Neu_s);    /* reactant 1: mantle species */
      Chem->Reactions[n].reactant[1] = -k-1;  /* reactant 2: ALL GRAINS(k+1) */
      Chem->Reactions[n].product[0] = j - Chem->ManInd; /* product 1: neutral */
      Chem->Reactions[n].product[1] = -k-1;   /* product 2: ALL GRAINS(k+1) */
                                              /* product 3-4: none */
      /* alpha is neutral mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
      /* beta is binding energy */
      Chem->Reactions[n].coeff[0].beta = Chem->Species[j].Eb;
      Chem->Reactions[n].use = 1;

    }
  }

	for (j=Chem->ManInd+Chem->N_Neu_f; j<Chem->ManInd+Chem->N_Neu_f+Chem->N_Neu;j++)
    for (k=0; k<Chem->NGrain; k++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 4;           /* Type is desorption */
      Chem->Reactions[n].reactant[0] =
				j+k*(Chem->N_Neu_f+Chem->N_Neu+Chem->N_Neu_s);  /* reactant 1: mantle species */
      Chem->Reactions[n].reactant[1] = -k-1;  /* reactant 2: ALL GRAINS(k+1) */
      Chem->Reactions[n].product[0] = 
				j - Chem->ManInd-Chem->N_Neu_f 
				+ 2*Chem->N_Neu_f+1; /* product 1: neutral */
      Chem->Reactions[n].product[1] = -k-1;   /* product 2: ALL GRAINS(k+1) */
                                              /* product 3-4: none */
      /* alpha is neutral mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
      /* beta is binding energy */
      Chem->Reactions[n].coeff[0].beta = Chem->Species[j].Eb;
      Chem->Reactions[n].use = 1;
    }

	for (j=Chem->ManInd+Chem->N_Neu_f+Chem->N_Neu; 
			j<Chem->ManInd+Chem->N_Neu_f+Chem->N_Neu+Chem->N_Neu_s;j++)
    for (k=0; k<Chem->NGrain; k++)
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 4;           /* Type is desorption */
      Chem->Reactions[n].reactant[0] =
				j+k*(Chem->N_Neu_f+Chem->N_Neu+Chem->N_Neu_s);  /* reactant 1: mantle species */
      Chem->Reactions[n].reactant[1] = -k-1;  /* reactant 2: ALL GRAINS(k+1) */
      Chem->Reactions[n].product[0] =
				j - Chem->ManInd-Chem->N_Neu_f + 
				2*Chem->N_Neu_f+Chem->N_Neu+1; /* product 1: neutral */
      Chem->Reactions[n].product[1] = -k-1;   /* product 2: ALL GRAINS(k+1) */
                                              /* product 3-4: none */
      /* alpha is neutral mass! */
      Chem->Reactions[n].coeff[0].alpha = Chem->Species[j].mass;
      /* beta is binding energy */
      Chem->Reactions[n].coeff[0].beta = Chem->Species[j].Eb;
      Chem->Reactions[n].use = 1;

    }

  /* Grain-Grain reaction: Chem->GrCharge * Chem->GrCharge reactions */
  m = Chem->GrInd + Chem->GrCharge;    /* Corresponds to first neutral charge */
  for (k=0; k<Chem->NGrain; k++)
  for (j=m+1; j<=m+Chem->GrCharge; j++) /* positive charge */
  for (l=0; l<Chem->NGrain; l++)
  for (i=m-Chem->GrCharge; i<=m-1; i++) /* negative charge */
  {
    if (k != l) /* collision of two different types of charges */
      for (h=0; h<=1; h++)
      {
        n = Chem->NReaction;
        InsertReactionInit(Chem);
        Chem->Reactions[n].rtype = 5;        /* Type is grain-grain reaction */
        Chem->Reactions[n].reactant[0] =
                             j + k * Ns_Gr;  /* reactant 1: positive charge */
        Chem->Reactions[n].reactant[1] =
                             i + l * Ns_Gr;  /* reactant 2: negative charge */

        if (h == 0) /* case 1 */
        {
          Chem->Reactions[n].product[0] =
                             m + k * Ns_Gr;  /* product 1: neutral charge */
          Chem->Reactions[n].product[1] =
                         i+j-m + l * Ns_Gr;  /* product 2: net charge */
          /* Probability of this reaction (brancing ratio) */
          Chem->Reactions[n].coeff[0].alpha =
                            Chem->GrSize[l]/(Chem->GrSize[l]+Chem->GrSize[k]);
        }

        if (h == 1) /* case 2 */
        {
          Chem->Reactions[n].product[0] =
                             m + l * Ns_Gr;  /* product 1: neutral charge */
          Chem->Reactions[n].product[1] =
                         i+j-m + k * Ns_Gr;  /* product 2: net charge */
          /* Probability of this reaction (branching ratio) */
          Chem->Reactions[n].coeff[0].alpha =
                            Chem->GrSize[k]/(Chem->GrSize[l]+Chem->GrSize[k]);
        }
                                             /* product 3-4: none */
        Chem->Reactions[n].use = 1;
      }
    else /* collision between two grains of the same type */
    {
      n = Chem->NReaction;
      InsertReactionInit(Chem);
      Chem->Reactions[n].rtype = 5;         /* Type is grain-grain reaction */
      Chem->Reactions[n].reactant[0] =
                            j + k * Ns_Gr;  /* reactant 1: positive charge */
      Chem->Reactions[n].reactant[1] =
                            i + k * Ns_Gr;  /* reactant 2: negative charge */
      Chem->Reactions[n].product[0] =
                            m + k * Ns_Gr;  /* product 1: neutral charge */
      Chem->Reactions[n].product[1] =
                        i+j-m + k * Ns_Gr;  /* product 2: net charge */
                                            /* product 3-4: none */
      Chem->Reactions[n].coeff[0].alpha = 1.0;   /* Probability=1 */
      Chem->Reactions[n].use = 1;
    }
  }

/*----------- Read and construct grain-surface reactions ----------*/

  if (Chem->NGrain > 0)
  {
    fgets(line,MAXLEN,fp);
    fscanf(fp, "%d\n", &k);       /* Number of grain-surface reactions; */
    if (k < 0)
      ath_error("[init_reactions]: Number of grain-surface reactions must be >=0!\n");

    fgets(line,MAXLEN,fp);
    fgets(line,MAXLEN,fp);

    for (i=0; i<k; i++)
    {
      fscanf(fp, "%s", name1);
      fscanf(fp, "%s", name2);
      fscanf(fp, "%s", name3);
      fscanf(fp, "%s", name4);
      fscanf(fp, "%lf", &coef);

      for (j=0;j<Chem->NGrain;j++)
      {
        /* First include surface reaction whose products are mantal species */
        sprintf(mantle,"[m%d]",j+1);
        n = Chem->NReaction;
        InsertReactionInit(Chem);

        Chem->Reactions[n].rtype = 6;   /* Type is grain surface reaction */
        /* reactants are mantle species */
        strcpy(name, name1); strcat(name,mantle);
        Chem->Reactions[n].reactant[0] = FindSpecies(Chem,name); 
        if (Chem->Reactions[n].reactant[0] == -10)
          ath_perr(-1,"[init_reactions]: Error in grain-surface reaction %d!\n", i+1);
        strcpy(name, name2); strcat(name,mantle);
        Chem->Reactions[n].reactant[1] = FindSpecies(Chem,name);
        if (Chem->Reactions[n].reactant[1] == -10)
          ath_perr(-1,"[init_reactions]: Error in grain-surface reaction %d!\n", i+1);
        /* products are mental species */
        strcpy(name, name3); strcat(name,mantle);
        Chem->Reactions[n].product[0]  = FindSpecies(Chem,name);
        strcpy(name, name4); strcat(name,mantle);
        Chem->Reactions[n].product[1]  = FindSpecies(Chem,name);
        Chem->Reactions[n].product[2]  = -10;
        Chem->Reactions[n].product[3]  = -10;
        /* relevant coeffieicnts */
        Chem->Reactions[n].coeff[0].alpha = coef; /* activation energy Ea (in Kelvin) */
        Chem->Reactions[n].coeff[0].beta  = (Real)(j); /* grain index */;
        Chem->Reactions[n].coeff[0].gamma = 0.90;  /* branching ratio is 90% */
        Chem->Reactions[n].use = 1;

        if (CheckReaction(Chem, Chem->Reactions[n]) < 0)
        {
          PrintReaction(Chem, n, 0.0);
          ath_perr(-1,"[init_reactions]: Error in grain-surface reaction %d!\n", i+1);
        }

        /* Next include the same reaction with product in the gas phase */
        n = Chem->NReaction;
        InsertReactionInit(Chem);

        Chem->Reactions[n].rtype = 6;   /* Type is grain surface reaction */
        /* reactants are mantle species */
        Chem->Reactions[n].reactant[0] = Chem->Reactions[n-1].reactant[0];
        Chem->Reactions[n].reactant[1] = Chem->Reactions[n-1].reactant[1];
        /* products are gas-phase species */
        Chem->Reactions[n].product[0]  = FindSpecies(Chem,name3);
        Chem->Reactions[n].product[1]  = FindSpecies(Chem,name4);
        Chem->Reactions[n].product[2]  = -10;
        Chem->Reactions[n].product[3]  = -10;
        /* relevant coeffieicnts */
        Chem->Reactions[n].coeff[0].alpha = coef; /* activation energy Ea (in Kelvin) */
        Chem->Reactions[n].coeff[0].beta  = (Real)(j); /* grain index */;
        Chem->Reactions[n].coeff[0].gamma = 0.10;  /* branching ratio is 10% */
        Chem->Reactions[n].use = 1;
      }
    } 
  }


  fclose(fp);

/*----------------------- End of constructing reactions ----------------------*/

  /* Find inverse of all reactions */
  m = FindInverse(Chem);

  /* Output All the reactions for checking purpose */
  n = Chem->NReaction;

  ath_pout(0,"Number of reactions: %d\n", n);
  ath_pout(0,"Reversible reaction pairs: %d\n", m);
  ath_pout(0,"\n");
  ath_pout(0,"List of reactions:\n");
  ath_pout(0,"\n");

  for (i=0; i<n; i++)
    PrintReaction(Chem,i,Chem->Reactions[i].coeff[0].alpha);


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
  React->reactant[0] = -10;
  React->reactant[1] = -10;
  React->product[0] = -10;
  React->product[1] = -10;
  React->product[2] = -10;
  React->product[3] = -10;
  React->coeff[0].alpha = 0.0;
  React->coeff[0].beta = 0.0;
  React->coeff[0].gamma = 0.0;
  React->coeff[0].Tmin = 0.0;
  React->coeff[0].Tmax = 0.0;
  React->NumTRange = 1;
  React->inverse = -1;

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

  for (i=0; i<Chem->N_Ele_tot; i++)
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

/*----------------------------------------------------------------------------*/
/* Find inverse of all reactions
 */
int FindInverse(Chemistry *Chem)
{
  int i, j, r, n, ni;

  n = Chem->NReaction;
  ni = 0;  /* number of reversible reaction pairs */

  for (i=0; i<n; i++)
    if (Chem->Reactions[i].inverse < 0)
    {
      j = i+1;
      r = 0;
      while ((j<n) && (r == 0))
      {
        r = CheckInverse(Chem, i, j);
        if (r == 1)
        {
          Chem->Reactions[i].inverse = j;
          Chem->Reactions[j].inverse = i;
        }
        j++;
      }
      if (r == 1) ni ++;
    }

  return ni;
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


/*============================================================================*/
/*----------------------------- Private Functions ----------------------------*/

/*----------------------------------------------------------------------------*/
/* Check whether reaction i and j are inverse reactions (1: Yes; 0: No)
 */
int CheckInverse(Chemistry *Chem, int i, int j)
{
  if (Chem->Reactions[i].product[2] >= 0) return 0;  /* More than 2 products */
  if (Chem->Reactions[j].product[2] >= 0) return 0;  /* More than 2 products */

  if (  (Chem->Reactions[i].reactant[0] == Chem->Reactions[j].product[0])
     && (Chem->Reactions[i].reactant[1] == Chem->Reactions[j].product[1]))
    return 1;

  if (  (Chem->Reactions[i].reactant[0] == Chem->Reactions[j].product[1])
     && (Chem->Reactions[i].reactant[1] == Chem->Reactions[j].product[0]))
    return 1;

  return 0;
}


/*---------------------------------------------------------------------------*/
/* Add another temperature range of the nth reaction
 */
void AddReactionTRange(Chemistry *Chem, int n)
{
  int m = Chem->Reactions[n].NumTRange;

  Chem->Reactions[n].NumTRange += 1;
  Chem->Reactions[n].coeff = (Coefficient*)realloc(Chem->Reactions[n].coeff,
                                                   (m+1)*sizeof(Coefficient));

  return;
}


/*----------------------------------------------------------------------------*/
/* Compare whether the nth reaction is the same as the previous one
 */
int ReactionCmp(Chemistry *Chem, int n)
{
  int i, p, q;

  if (n == 0)
    return -1;  /* Not the same */

  for (i=0; i<2; i++)
  {
    p = Chem->Reactions[n-1].reactant[i];
    q = Chem->Reactions[n].reactant[i];
    if (strcmp(Chem->Species[p].name, Chem->Species[q].name) != 0)
      return -1;
  }

  for (i=0; i<4; i++)
  {
    p = Chem->Reactions[n-1].product[i];
    q = Chem->Reactions[n].product[i];
    if (strcmp(Chem->Species[p].name, Chem->Species[q].name) != 0)
      return -1;
  }

  return 0;     /* They are the same */
}

#undef DNR

#endif /* CHEMISTRY */

