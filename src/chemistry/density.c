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

  FILE *fp;
  Real sum = 0.0;
  char line[MAXLEN],fname[20],name[NL_SPE];
  Chemistry   *Chem = Evln->Chem;
  ElementInfo *Ele  = Chem->Elements;
  SpeciesInfo *Spe;
  int i,j,k,p,initcond,N;
  Real sumgas, grtot, ratio,abun,ChargeDen;

  Evln->rho = rho;
  Evln->t   = 0.0;
  sumgas    = 0.0;
  ChargeDen = 0.0;
  grtot     = 0.0;
  /* set initial abundance from input file initial.txt
   * if initcond is larger than zero */
  initcond  = par_getd_def("problem","initcond", 0);
  /* Calculate initial number density */
  ath_pout(verbose,"\n");
  ath_pout(verbose,"Initial number density:\n");
  ath_pout(verbose,"\n");

  /* initialize species # density as zero */
  for (i=0; i<Chem->Ntot; i++)
     Evln->NumDen[i] = 0.0;

  // initcond>0 user defined initial condition
  if(initcond>0) 
  {
  /* Reset element abundance to zero */
  for (i=0; i<Chem->N_Ele_tot; i++)
     Chem->Elements[i].abundance = 0.0;

  sprintf(fname,"%s",par_gets("job","init_abundance"));
  fp = fopen(fname,"r");
  if(fp == NULL)
    ath_perr(-1,"[construct_species]:Unable to open the inital abundance file!\n");
  
  fgets(line,MAXLEN,fp);
  /* N is total number of initial species */
  fscanf(fp, "%d\n", &N);
  if (N <= 0)
    ath_error("[density]: # of species with initial abundance must be positive!\n");

  fgets(line,MAXLEN,fp);
  for (i=0; i<N; i++)
  {
    fscanf(fp, "%s", name);
    fscanf(fp, "%lf\n",&abun);
    for (p=1;p<Chem->Ntot;p++) {
      if(strcmp(Chem->Species[p].name, name)==0)
      {
        /* read species relative abundance */
        Evln->NumDen[p] = abun;
        Spe = &(Chem->Species[p]);
        sumgas += Spe->mass*abun;

        /* calculate element abundance */
        for (j=0; j<Chem->N_Ele_tot; j++)
          Chem->Elements[j].abundance += abun*Spe->composition[j];
        break;
      }
    }// end loop for p
  }/* end iteration for i over all input species */

  /* Calculate the relative charge density */
  for (i=0; i<Chem->Ntot; i++)
    if ((i != 0) && (Chem->Species[i].charge != 0))
      ChargeDen += Evln->NumDen[i] * Chem->Species[i].charge;
  /* Make up for the charge density */
  Evln->NumDen[0] = ChargeDen;

  /* Recalculate grain abundance for user defined
   * initial condition 
   * */
  for (i=0; i<Chem->NGrain; i++)
    grtot += Chem->GrFrac[i];
  ratio = sumgas / (1.0-grtot);
  for (i=0; i<Chem->NGrain; i++)
  {
    k = Chem->N_Ele+i;  /* label of grain element */
    Chem->Elements[k].abundance
      = Chem->GrFrac[i]/Chem->Elements[k].mass*ratio;
  }

  /* Rescale relative # density to real # density*/
  sum = 0.0;
  for (i=0; i<Chem->N_Ele+Chem->NGrain; i++) {
    sum += Ele[i].mass * Ele[i].abundance;
  }
  Evln->Abn_Den = (sum*1.672e-24)/rho; /* == 1/n_H */

  /* calculate species number density */
  for (p=0;p<Chem->Ntot;p++)
  {
    Evln->NumDen[p] /= Evln->Abn_Den;
    Spe = &(Chem->Species[p]);
    ath_pout(verbose, "[%s]: = %e\n", Spe->name, Evln->NumDen[p]);
  }

  /* initialize dust grain with single element */
  for (i = Chem->N_Ele; i < Chem->N_Ele+Chem->NGrain; i++)
  {/* Initialize with the first single-element species of each element */
    k = Ele[i].single[0];
    Spe = &(Chem->Species[k]);
    Evln->NumDen[k] = Ele[i].abundance / (Evln->Abn_Den * Spe->composition[i]);
    ath_pout(verbose, "[%s]: = %e\n", Spe->name, Evln->NumDen[k]);
   }

  /* output element abundance for comparison */
  for (j=0; j<Chem->N_Ele_tot; j++)
    ath_pout(0,"%3s abundance: %10e\n",
      Chem->Elements[j].name,Chem->Elements[j].abundance/Evln->Abn_Den);

  /* initialize species as their first single element */
  }else{
   /* Calculate abundance - number density ratio */
   sum = 0.0;
   for (i=0; i<Chem->N_Ele+Chem->NGrain; i++) {
     sum += Ele[i].mass * Ele[i].abundance;
   }
   Evln->Abn_Den = (sum*1.672e-24)/rho; /* == 1/n_H */

    for (i = 0; i < Chem->N_Ele+Chem->NGrain; i++)
    {/* Initialize with the first single-element species of each element */

    k = Ele[i].single[0];
    Spe = &(Chem->Species[k]);
      Evln->NumDen[k] = Ele[i].abundance / (Evln->Abn_Den * Spe->composition[i]);
    ath_pout(verbose, "[%s]: = %e\n", Spe->name, Evln->NumDen[k]);
   }
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
