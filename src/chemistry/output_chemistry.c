#include "../header/copyright.h"
/*=============================================================================
 * FILE: output_chemistry.c
 *
 * PURPOSE: Contains all the output functions.
 *
 *   To set up an output, use
 *      init_chemout(Chemistry *Chem, ChemOutput *ChemOut, int mode,
 *                   char *bname, Real bvalue, char *id)
 *      where
 *        mode  : a number that should be 1, 2 or 3, corresponding to
              1 - will output the number density of all species
              2 - will output the magnetic diffusivities
              3 - will output the recombination time
 *        *bname: the base name of the output file 
 *        bvalue: an arbitrary user specified number that is useful for
 *                identifying what is being output in this file (e.g., radius
 *                in the disk, etc.).
 *        *id   : the id number of the file (in case of mutiple outputs)
 *
 *   To execute the output, use one of the following:
 *     output_nspecies() -> mode=1
 *     output_etaB()     -> mode=2
 *     output_recomb()   -> mode=3
 *     output_etafit()   -> mode=4
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  - init_chemout()    - initialize the output structure
 *  - final_chemout()   - finalize the output structure
 *  - output_nspecies() - output the number density of selected species (mode=1)
 *  - output_etaB()     - output the magnetic diffusivities (mode=2)
 *  - output_etafit()   - output the fitting prameters to diffusivities (mode=4)
 *  - output_recomb()   - output the recombination time (mode=3)
 *  - ChemSet_allspecies() - select all species
 *  - ChemSet_allgrain()   - select all grain species
 *  - ChemSet_allgas()     - select all gas-phase species
 *  - ChemSet_selected()   - select species from reading input file
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

#define NCOL 7

/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*------------------------------------------------------------------------------
 * Initiate the output structure
 */
void init_chemout(Chemistry *Chem, ChemOutput *ChemOut, int mode, 
                                   char *bname, Real bvalue, char *id)
{
  int i;

  sprintf(ChemOut->outbase, "%s-%s",
                            par_gets("job","outbase"),par_gets("job","outid"));

  ChemOut->mode    = mode;
  ChemOut->lab     = 0;
  ChemOut->nsp     = 0;

  strcpy(ChemOut->bname, bname);
  ChemOut->bvalue = bvalue;

  if (mode == 1)
  { /* species number densities */
    sprintf(ChemOut->outid,"nsp-%s",id);

    ChemOut->ind = (int*)calloc_1d_array(Chem->Ntot, sizeof(int));

    for (i=0; i<Chem->Ntot; i++)
      ChemOut->ind[i] = i;
  }
  else if (mode == 2)
  { /* magnetic diffusivities */
    sprintf(ChemOut->outid,"eta-%s",id);

    ChemOut->ind = NULL;
  }
  else if (mode == 3)
  { /* recombination time */
    sprintf(ChemOut->outid,"rcb-%s",id);

    ChemOut->ind = NULL;
  }
  else if (mode == 4)
  { /* magnetic diffusion fitting */
    sprintf(ChemOut->outid,"fit-%s",id);

    ChemOut->ind = NULL;
  }
  else
  {
    ath_error("[init_chemout]: Output mode must be equal to 1, 2 or 3!\n");
  }

  return;
}

/*------------------------------------------------------------------------------
 * Finalize the output structure
 */
void final_chemout(ChemOutput *ChemOut)
{

  if (ChemOut->ind != NULL)
    free(ChemOut->ind);

  return;
}


/*------------------------------------------------------------------------------
 * Output the number density of selected species
 * (e.g., as a function of time or position)
 */
void output_nspecies(ChemEvln *Evln, ChemOutput *ChemOut,
                                           char *pname,   Real value)
{
  int i, j, rem;
  Chemistry *Chem = Evln->Chem;

  char fname[50];
  FILE *fp;

  if (ChemOut->mode != 1) {
    ath_error("[output_chem]: Outputing number densities requires mode = 1!\n");
  }

/* open the file */
  sprintf(fname,"%s.%s.dat", ChemOut->outbase, ChemOut->outid);

  if (ChemOut->lab == 0) {
    if ((fp = fopen(fname,"w")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }
  else {
    if ((fp = fopen(fname,"a+")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }

/* print the header */
  if (ChemOut->lab == 0) {
    fprintf(fp,"# %s: %10e\n", ChemOut->bname, ChemOut->bvalue);
    fprintf(fp,"# Number density as a function of %4s\n", pname);
    fprintf(fp,"# %d species in total\n",ChemOut->nsp);
    fprintf(fp,"# Species list has %d lines", (ChemOut->nsp-1)/NCOL +1);

    for (i=0; i<ChemOut->nsp; i++)
    {
      rem = i - (i/NCOL)*NCOL;

      if (rem == 0) {
        fprintf(fp,"\n#");
      }

      j = ChemOut->ind[i];
      if (rem == 0){
        fprintf(fp," %8s   ", Chem->Species[j].name);
      }else{
        fprintf(fp,"%10s   ", Chem->Species[j].name);
      }
    }

    fprintf(fp,"\n");
  }

  ChemOut->lab++;

/* print the data */
  fprintf(fp,"%10e ", value);
  fprintf(fp,"%10e",Evln->zeta_eff);

  for (i=0; i<ChemOut->nsp; i++) 
  { 
    rem = i - (i/NCOL)*NCOL;
    if (rem == 0) {
      fprintf(fp,"\n");
    }

    j = ChemOut->ind[i];
    fprintf(fp,"%10e ", Evln->NumDen[j]);    
  }

  fprintf(fp,"\n");

  fclose(fp);

  return;
}

/*------------------------------------------------------------------------------
 * Read the number density of selected species from output
 * (e.g., as a function of time or position)
 */
void read_nspecies(ChemEvln *Evln, ChemOutput *ChemOut, int nskip)
{
  int i, j, myind;
  Chemistry *Chem = Evln->Chem;
  double value,drho,drhomax,ChargeDen;
  char line[MAXLEN],fname[50];
  FILE *fp;

  if (ChemOut->mode != 1) {
    ath_error("[output_chem]: Outputing number densities requires mode = 1!\n");
  }

/* open the file */
  sprintf(fname,"%s.%s.dat", ChemOut->outbase, ChemOut->outid);

  if ((fp = fopen(fname,"r")) == NULL)
    ath_error("[output_chem]: Error opening file %s...\n", fname);

/* read header */
  for (i=0;i<(ChemOut->nsp-1)/NCOL+5;i++)
    fgets(line,MAXLEN,fp);

/* skip early data */
  for (i=0; i<(ChemOut->nsp+2)*nskip; i++)
    fscanf(fp,"%lf",&value);

/* read data */
  fscanf(fp,"%lf",&value);
  fscanf(fp,"%lf",&(Evln->zeta_eff));

  for (i=0; i<ChemOut->nsp; i++)
  {
    j = ChemOut->ind[i];
    fscanf(fp,"%lf", &(Evln->NumDen[j]));
  }

  fclose(fp);

/* make adjustment for charge conservation 
 * charge does not conserve due to floating point error in the output */
  drhomax = 0.0;
  myind = 0;
  ChargeDen=0.0;
  for (i=0; i<Chem->Ntot; i++)
  {
    if (Chem->Species[i].charge != 0)
    {
      drho = Evln->NumDen[i] * Chem->Species[i].charge;
      ChargeDen += drho;
      if (fabs(drho) > drhomax){
        drhomax = drho;
        myind = i;
      }
    }
  }

  Evln->NumDen[myind] -= ChargeDen/Chem->Species[myind].charge;

  return;
}

/*------------------------------------------------------------------------------
 * Output the magnetic diffusivities as a function of B at fixed time/position
 * (i.e., at chemical equilibrium)
 */
void output_etaB(ChemEvln *Evln, ChemOutput *ChemOut, char *pname, Real value,
                                     Real Omega, int nB)
{
  int i;
  Real Be,Bi,Bmin,Bmax,rho,rhomin,rhomax;
  Real dlnB, B, vA2, Cs2, eta0;
  Chemistry *Chem = Evln->Chem;

  char fname[50];
  FILE *fp;

/* Figure out min and max B */
  rhomin = par_getd("problem","rhomin");
  rhomax = par_getd("problem","rhomax");
  /* calculate the transition field strength  */
  Be = 1.2e-3 *rhomax/1.0e-11*MIN(1.0,sqrt(100/Evln->T));
  Bi = 0.82   *rhomax/1.0e-11;

  Bmin = 1.e-5 * Be;
  Bmax = 100000.0* Bi;
  //Bmin = par_getd("problem","Bmin");
  //Bmax = par_getd("problem","Bmax");
/* checkpoints */
  if (ChemOut->mode != 2) {
    ath_error("[output_chem]: Outputing etas requires mode = 2!\n");
  }

  if (nB <= 0) {
    ath_error("[output_chem]: nB must be positive!\n");
  }

/* open the file */
  sprintf(fname,"%s.%s.dat", ChemOut->outbase, ChemOut->outid);

  if (ChemOut->lab == 0) {
    if ((fp = fopen(fname,"w")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }
  else {
    if ((fp = fopen(fname,"a+")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }

/* print the header */
  if (ChemOut->lab == 0) {
    fprintf(fp,"# %s: %10e\n", ChemOut->bname, ChemOut->bvalue);
    fprintf(fp,"# nB = %d\n", nB);
    fprintf(fp,"# %5s         B-field  ", pname);
    fprintf(fp," eta_O        eta_H        eta_A\n");
  }

  ChemOut->lab++;

/* print the data */
  dlnB = log(Bmax/Bmin)/nB;

  for (i=0; i<nB; i++)
  {
    B = Bmin * exp((i+0.5)*dlnB);

    fprintf(fp,"%10e %10e ", value, B);

    /* calculate magnetic diffusivities */
    Evln->B = B;
    Cal_NIMHD(Evln);

    fprintf(fp,"%10e %10e %10e\n", Evln->eta_O, Evln->eta_H, Evln->eta_A);
  }

  fclose(fp);

  return;
}

/*------------------------------------------------------------------------------
 * Output the magnetic diffusivities as a function of B at fixed time/position
 * (i.e., at chemical equilibrium)
 */
void output_etafit(ChemEvln *Evln, ChemOutput *ChemOut,
                                     char *pname, Real value)
{
  int i;
  Real B, Bi, Be;    /* Bi: B for betai=1; Be: B for betae=1 */
  Real Q_H1,Q_H2,Q_A1,Q_A2; /* eta_H=Q_H*B, eta_A=Q_A*B^2 */
  Chemistry *Chem = Evln->Chem;

  char fname[50];
  FILE *fp;

/* checkpoints */
  if (ChemOut->mode != 4) {
    ath_error("[output_chem]: Outputing etaB-fit requires mode = 4!\n");
  }

/* open the file */
  sprintf(fname,"%s.%s.dat", ChemOut->outbase, ChemOut->outid);

  if (ChemOut->lab == 0) {
    if ((fp = fopen(fname,"w")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }
  else {
    if ((fp = fopen(fname,"a+")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }

/* print the header */
  if (ChemOut->lab == 0) {
    fprintf(fp,"# %s: %10e\n", ChemOut->bname, ChemOut->bvalue);
    fprintf(fp,"# rho: %10e\n",Evln->rho);
    fprintf(fp,"# %5s        B_i          ", pname);
//    fprintf(fp,"# %5s      ", pname);
    fprintf(fp,"eta_O        Q_H1         Q_A1         Q_H2         Q_A2\n");
  }

  ChemOut->lab++;

/* calculate the transition field strength  */
  Be = 1.2e-3 * Evln->rho/1.0e-11*MIN(1.0,sqrt(100/Evln->T));
  Bi = 0.82   * Evln->rho/1.0e-11;

  fprintf(fp,"%10e %10e ", value, Bi);
//  fprintf(fp,"%10e ", value);

/* calculate magnetic diffusivities in the Ohmic regime */  
  B = 0.03*Be;
  Evln->B = B;
  Cal_NIMHD(Evln);
    
  fprintf(fp,"%10e %10e %10e ", Evln->eta_O,Evln->eta_H/B,Evln->eta_A/SQR(B));

/* calculate magnetic diffusivities in the AD regime */
  B = 300.0*Bi;
  Evln->B = B;
  Cal_NIMHD(Evln);
     
  fprintf(fp,"%10e %10e\n", Evln->eta_H/B, Evln->eta_A/SQR(B));
  
  fclose(fp);
  
  return;
} 


/*------------------------------------------------------------------------------
 * Output the magnetic diffusivities at fixed B
 * (e.g., as a function of time or position)
 */
void output_recomb(ChemEvln *Evln, ChemOutput *ChemOut, char *pname, Real value,
                                   Real Omega, Real Bmin, Real Bmax, int nB)
{
  int i;
  Real Dt, dlnB, B, t_O, *t_H, *t_A;
  Chemistry *Chem = Evln->Chem;

  char fname[50];
  FILE *fp;

/* checkpoints */
  if (ChemOut->mode != 3) {
    ath_error("[output_chem]: Outputing t_recomb requires mode = 3!\n");
  }

  if ((Bmax < Bmin) || (Bmin <= 0)) {
    ath_error("[output_chem]: Bmax >= Bmin > 0 must be satisfied!\n");
  }

  if (nB <= 0) {
    ath_error("[output_chem]: nB must be positive!\n");
  }

/* initialization */
  Dt = 1.0e-7 * OneYear;

  t_H = (Real*)calloc_1d_array(nB, sizeof(Real));
  t_A = (Real*)calloc_1d_array(nB, sizeof(Real));

/* open the file */
  sprintf(fname,"%s.%s.dat", ChemOut->outbase, ChemOut->outid);

  if (ChemOut->lab == 0) {
    if ((fp = fopen(fname,"w")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }
  else {
    if ((fp = fopen(fname,"a+")) == NULL)
      ath_error("[output_chem]: Error opening file %s...\n", fname);
  }

/* print the header */
  if (ChemOut->lab == 0) {
    fprintf(fp,"# %s: %lf\n", ChemOut->bname, ChemOut->bvalue);
    fprintf(fp,"# nB = %d\n", nB);
    fprintf(fp,"# %5s         B-field      t_dyn       ", pname);
    fprintf(fp," t_O          t_H          t_A\n");
  }

  ChemOut->lab++;

/* calculate the recombination time */
  Cal_recomb(Evln, Bmin, Bmax, nB, Dt, &t_O, t_H, t_A);

/* print the data */
  dlnB = log(Bmax/Bmin)/nB;

  for (i=0; i<nB; i++)
  {
    B = Bmin * exp((i+0.5)*dlnB);

    fprintf(fp,"%12e %12e %12e ", value, B, 1.0/Omega);

    fprintf(fp,"%12e %12e %12e\n", t_O, t_H[i], t_A[i]);

  }

  free_1d_array(t_H);
  free_1d_array(t_A);

  fclose(fp);

  return;
}

/*------------------------------------------------------------------------------
 * Auxilary routine for output_nspecies:
 *   Get the indices array for all species
 */
void ChemSet_allspecies(Chemistry *Chem, ChemOutput *ChemOut)
{
  ChemOut->nsp = Chem->Ntot;

  return;
}

/*------------------------------------------------------------------------------
 * Auxilary routine for output_nspecies:
 *   Get the indices array for all grain species
 */
void ChemSet_allgrain(Chemistry *Chem, ChemOutput *ChemOut)
{
  int i;

  ChemOut->nsp = Chem->Ntot - Chem->GrInd;

  for (i=0; i<ChemOut->nsp; i++)
  {
    ChemOut->ind[i] = i + Chem->GrInd;
  }

  return;
}

/*------------------------------------------------------------------------------
 * Auxilary routine for output_nspecies:
 *   Get the indices array for all gas phase species
 */
void ChemSet_allgas(Chemistry *Chem, ChemOutput *ChemOut)
{
  ChemOut->nsp = Chem->GrInd-1;

  return;
}

/*------------------------------------------------------------------------------
 * Auxilary routine for output_nspecies:
 *   Get the indices array for selected species (read from the input file)
 */
void ChemSet_selected(Chemistry *Chem, ChemOutput *ChemOut)
{
  int i, n=0;

  for (i=0; i<Chem->Ntot; i++)
  {
    if (par_exist("out_species",Chem->Species[i].name))
    {
      ChemOut->ind[n] = i;
      n++;
    }
  }

  ChemOut->nsp = n;

  return;
}

#undef NCOL

#endif /* CHEMISTRY */

