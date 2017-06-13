/*=============================================================================
 * FILE: evolve.c
 *
 * PURPOSE: Contains functions to evolve the number densities of all species
 *   based on the chemistry model.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   evolve() - evolve the number densities for a given time period
 *   jacobi() - calculate the Jacobi matrix of the ODEs
 *   derivs() - calculate the rate of density change
 *   EleMakeup() - density makeup for charge/element conservation
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
#include <sys/types.h>
#include <time.h>
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/prototypes.h"
#include "../header/chemproto.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   ???Makeup()    - Makeup charge ane element densities for conservation
 *============================================================================*/

int EleMakeup_sub(ChemEvln *Evln, int q, Real dn);
int ChargeMakeup (ChemEvln *Evln, Real dne);
void cparray(Real *ary1, Real *ary2, int N);
int AdDesRate(ChemEvln *Evln);

/*============================================================================*/
/*---------------------------- Public Functions ------------------------------*/

/*----------------------------------------------------------------------------*/
/* Evolve the equation using numerical recipes routine stiff
 * Input parameters:
 *   te:       evolution time (s)
 *   dttry:    trial time step (s)
 *   err:      error level (%)
 * Return:
 *   status (0: good; -1: error)
 */
int evolve(ChemEvln *Evln, Real te, Real *dttry, Real err)
{
  int i, k, status, verbose, p, p1;
  Real *numden, *dn_o_dt, *numden1;
  Real dt, dtn, t=0.0, tp;
  Real myerr;
  clock_t c0, c1; /* Timing the code */
  Chemistry *Chem = Evln->Chem;
  Real Tn;

  dn_o_dt  = (Real*)calloc_1d_array(Chem->Nsp, sizeof(Real));
  numden   = (Real*)calloc_1d_array(Chem->Nsp, sizeof(Real));
  if (Chem->NGrain>0){
    status = AdDesRate(Evln);
    numden1  = (Real*)calloc_1d_array(2*Chem->NNeuT, sizeof(Real));
  }

  /* evolve the number densities */
  c0 = clock();
  ath_pout(0,"\n");
  ath_pout(0,"Chemical evolution started...\n");
  ath_pout(0,"At t=%e yr, Abn(e-)=%e, next dt=%e yr.\n",
          Evln->t/OneYear, Evln->NumDen[0]*Evln->Abn_Den, *dttry/OneYear);

  tp = Evln->t*1.5; /* to control the output of log information */
  int ind=0;
  Real dt1 =0;
  while (t < te)
  {
    *dttry = MIN(*dttry, te - t); /* timing control */

    myerr = err;

    derivs(Evln, Evln->NumDen, dn_o_dt); /* Reaction rates */
    cparray(Evln->NumDen, numden, Chem->Nsp);

    status = stifbs(Evln, numden, dn_o_dt,Chem->Nsp, &t, *dttry, myerr,
                                      Evln->DenScale, &dt, &dtn);
    
    if (status != 0) {
      // if it fails, try another solver...
      cparray(Evln->NumDen, numden, Chem->Nsp);
      status = stifkr(Evln, numden, dn_o_dt, Chem->Nsp, &t, *dttry, myerr,
                                      Evln->DenScale, &dt, &dtn);
    }

    if (status != 0)
    {
      ath_pout(0, "At t=%e yr, calculation fails...\n",Evln->t/OneYear);
      break;
    }

    cparray(numden, Evln->NumDen, Chem->Nsp);

    Evln->t += dt;

    ath_pout(0,"Evln.t=%2.4f\n",Evln->t/OneYear);
    /* update numden of mantle species */
    ind += 1;
    dt1 += dt;
    if (Chem->NGrain>0 && (ind>0 || dt1/OneYear>1e3) )
    {
      ath_pout(0,"Evln.t=%2.4f, Dt=%10e\n",Evln->t/OneYear,dt1/OneYear);
      for (k=0;k<Chem->NNeuT;k++){
	p   = ConvertInd(k);
	p1  = ConvertInd(k+Chem->NNeuT);
	Tn  = Evln->NumDen[p]+Evln->NumDen[p1];
	numden1[k]  = (Evln->NumDen[p]-Tn*Chem->AdDesRate[k].des
          /(Chem->AdDesRate[k].ad+Chem->AdDesRate[k].des))*
	  exp(-(Chem->AdDesRate[k].ad+Chem->AdDesRate[k].des)*dt1)+
	  Tn*Chem->AdDesRate[k].des/(Chem->AdDesRate[k].des+Chem->AdDesRate[k].ad);
	numden1[k+Chem->NNeuT] = Tn-numden1[k];
      }
      for (k=0;k<2*Chem->NNeuT;k++){
        p = ConvertInd(k);
        Evln->NumDen[p] = numden1[k];
      }
      ind=0;
      dt1 = 0.0;
    }

    if (Evln->t > tp) /* verbose control */
    {
      verbose = 0;
      tp = Evln->t * 1.5;
    } else
    {
      verbose = 1;
    }

    EleMakeup(Evln, verbose); /* Imposing conservation laws */

    *dttry = dtn;  /* recommended time step for the next cycle */

    ath_pout(verbose, "At t=%e yr, Abn(e-)=%e, next dt=%e yr.\n",
                   Evln->t/OneYear, Evln->NumDen[0]*Evln->Abn_Den, dtn/OneYear);

    c1 = clock();

    /* if evolution is too time consuming, or element makeup fails, quit */
    if ((c1-c0)/CLOCKS_PER_SEC > 2*3600.0)
    {
      status = 1;
      break; 
    }
  }

/* finalize and return the status */

  free_1d_array(dn_o_dt);
  free_1d_array(numden);

  if ((status >=0) || (Evln->t > 0.5*te))
    ath_pout(0,"Evolution completed at t=%e yr, with Abn(e-)=%e.\n",
                               Evln->t/OneYear, Evln->NumDen[0]*Evln->Abn_Den);
  else
    ath_pout(0,"Evolution terminated at t=%e yr, with Abn(e-)=%e.\n",
                               Evln->t/OneYear, Evln->NumDen[0]*Evln->Abn_Den);
  ath_pout(0,"\n");

  return status;
}

/*----------------------------------------------------------------------------*/
/* user provided routine for calculating the Jacobi matrix
 */
void jacobi(ChemEvln *Evln, Real *numden, Real **jacob)
{
  int i, j, k, l, n, p, a;
  Real Jt;
  Chemistry *Chem = Evln->Chem;
  EquationTerm *EqTerm;
  /* Initialization */
  for (i=0; i<Chem->Nsp; i++) {
  for (j=0; j<Chem->Nsp; j++) {
    jacob[i][j] = 0.0;
  }}

  /* Calculation */
  for (i=0; i<Chem->Nsp; i++)  /* loop over all species */
  {
    n = Chem->Equations[i].NTerm;

    for (j=0; j<n; j++)	/* loop over all reactions of this species */
    {

      EqTerm = &(Chem->Equations[i].EqTerm[j]);

      for (k=0; k<EqTerm->N; k++)  /* loop over all reactants */
      {
        if ((EqTerm->type) != 3 && (EqTerm->type !=4)){
          Jt = Evln->K[EqTerm->ind] * EqTerm->dir;

          for (l=0; l<EqTerm->N; l++)
          {
            p = EqTerm->lab[l];

            if (l != k) Jt *= numden[p];
          }

          p = EqTerm->lab[k];
          if (p>=Chem->Nsp)
            ath_perr(0,"Error! index too large in Jacobian!\n");
          jacob[i][p] += Jt;  /* Obtain the Jacobian from this term */
	}
      }//end for k
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* user provided routine for calculating time derivatives of number density
 */
void derivs(ChemEvln *Evln, Real *numden, Real *drv)
{
  int i, j, k, p;
  Real sum;
  Real rate;
  Chemistry *Chem = Evln->Chem;
  EquationTerm *EqTerm;

  for (k=0; k<Chem->Nsp; k++)
  {
    sum  = 0.0;

    for (i=0; i<Chem->Equations[k].NTerm; i++)
    {
      EqTerm = &(Chem->Equations[k].EqTerm[i]);

      rate = Evln->K[EqTerm->ind] * EqTerm->dir;
      if ((EqTerm->type) != 3 && (EqTerm->type !=4)){
        for (j=0; j<EqTerm->N; j++)
        {
          p = EqTerm->lab[j];
          rate *= numden[p];
        }
        sum += rate;
      }
    }

    drv[k] = sum;
  }

  return;
}

/*---------------------------------------------------------------------------*/
/* Make up the element number density to enforce conservation laws
 */
int EleMakeup(ChemEvln *Evln, int verbose)
{
  int i, j, k, l, status=0;
  Real den, denmax, disp, frac;

  Real *EleNumDen;        /* Number density of each element */
  Real ChargeDen = 0.0;   /* Total charge number density excluding electron */

  Chemistry *Chem = Evln->Chem;
  Real    *NumDen = Evln->NumDen;

/* Initialization */

  EleNumDen = (Real*)calloc_1d_array((Chem->N_Ele+Chem->NGrain), sizeof(Real));

  for (i=0; i<Chem->N_Ele + Chem->NGrain; i++) {
    EleNumDen[i] = 0.0;
  }

/* For those with negative density, set them to zero */

  for (i=0; i<Chem->Ntot; i++) {
    if (NumDen[i] < 0.0)
    {
      ath_pout(verbose, "Warning: At t=%e yr, [%s] = %e < 0!\n",
                            Evln->t/OneYear, Chem->Species[i].name,NumDen[i]);
      NumDen[i] = 0.0;
    }
  }

/* Calculate the elemental density */

  for (i=0; i<Chem->Ntot; i++)
  {
    for (j=0; j<Chem->N_Ele+Chem->NGrain; j++)
    {
     if (Chem->Species[i].composition[j] > 0)
        EleNumDen[j] += NumDen[i]*Chem->Species[i].composition[j];
    }
  }

/* Make up for the element densities */

  for (i=0; i<Chem->N_Ele; i++)
  {
    /* Calculate the discrepency */

    disp = (EleNumDen[i] - Chem->Elements[i].abundance/Evln->Abn_Den);

    ath_pout(verbose,"Discrepancy for %3s : %e over %e\n",
      Chem->Elements[i].name, disp, Chem->Elements[i].abundance/Evln->Abn_Den);

    /* if abundance is smaller than the true value, then increase
     * the number densities of its single-element species
     */
    if (disp < 0.0)
    {
      den = 0.0;

      /* Calculate the element number density from single-element species */
      for (j=0; j<Chem->Elements[i].numsig; j++)
      {
        l = Chem->Elements[i].single[j];
        den += NumDen[l]*Chem->Species[l].composition[i];
      }
      frac = disp/den;

      for (j=0; j<Chem->Elements[i].numsig; j++)
      {
        l = Chem->Elements[i].single[j];

        NumDen[l] *= (1.0 - frac);
      }
    }
    /* if abundance is larger than the true value, then reduce
     * the number densities of all neutral species containing this element
     */
    else
    {
      status = EleMakeup_sub(Evln, i, disp);
    }

    if (status < 0)
      return status;
  }

/* Make up for the grain densities */
  for (i=Chem->N_Ele; i<Chem->N_Ele+Chem->NGrain; i++)
  {
    /* Calculate the discrepency */
    disp = (EleNumDen[i] - Chem->Elements[i].abundance/Evln->Abn_Den);

    ath_pout(verbose,"Discrepancy for %3s : %e over %e\n",
      Chem->Elements[i].name, disp, Chem->Elements[i].abundance/Evln->Abn_Den);

    /* Density make up */
    frac = disp / EleNumDen[i];

    for (j=Chem->GrInd; j<Chem->ManInd; j++)
    {
      if (Chem->Species[j].composition[i] > 0)
        NumDen[j] *= (1.0-frac);
    }
  }

/* Calculate the charge density */

  for (i=0; i<Chem->Ntot; i++)
  {
    if ((i != 0) && (Chem->Species[i].charge != 0))
      ChargeDen += NumDen[i] * Chem->Species[i].charge;
  }

/* Make up for the charge density */

  if (ChargeDen >= 0.0)
  {
    NumDen[0] = ChargeDen;
  }
  else
  {
    NumDen[0] = 0.0;

    status = ChargeMakeup(Evln, -ChargeDen);
  }

  free_1d_array(EleNumDen);

  return status;
}


/*============================================================================*/
/*------------------------------ PRIVATE FUNCTIONS ---------------------------*/

/*---------------------------------------------------------------------------*/
/* Density make up for single element species
 * q:  name of the element
 * dn: density makeup
 */
int EleMakeup_sub(ChemEvln *Evln, int q, Real dn)
{
  int i, j, k=0;
  Real den, frac, dni;
  Chemistry *Chem = Evln->Chem;
  SpeciesInfo *Species = Chem->Species;

/* Find the neutral species containing this element
 * and have the largest number density
 */
  den = 0.0;

  for (i=0; i<Chem->Ntot; i++)
  {
    if ((Species[i].composition[q] > 0) && (Species[i].charge == 0))
    {
      den += Evln->NumDen[i] * Species[i].composition[q];
    }
  }

  if (den < dn)
  {
    ath_perr(0,"Error! Can not make up for [%s]!\n", Chem->Elements[q].name);
    return -1;
  }

/* Density makeup */

  frac = dn / den;

  for (i=0; i<Chem->Ntot; i++)
  {
    if ((Species[i].composition[q] > 0) && (Species[i].charge == 0))
    {
      dni = Evln->NumDen[i] * frac;

      Evln->NumDen[i] *= (1.0-frac);

      for (j=0; j<Chem->N_Ele+Chem->NGrain; j++)
      {
        if (Species[i].composition[j] > 0)
        {
          if (j != q)
          {
            k = Chem->Elements[j].single[0];
            Evln->NumDen[k] +=
                  dni*Species[i].composition[j]/Species[k].composition[j];
          }
        }
      }
    }
  }

  return 0;
}

/*---------------------------------------------------------------------------*/
/* Density make up for the electrons
 */
int ChargeMakeup(ChemEvln *Evln, Real dne)
{
  int i, j;
  Real ratio, negchargetot, de;
  Real *negcharge;
  Chemistry *Chem= Evln->Chem;

  negchargetot = 0.0;
  negcharge = (Real*)calloc(Chem->NGrain, sizeof(Real));
  for (j=0; j<Chem->NGrain; j++)
    negcharge[j] = 0.0;

  /* calculate the total negative charge */
  for (i=Chem->GrInd; i<Chem->ManInd; i++)
  {
    /* get the index of grain type */
    j = (i-Chem->GrInd)/(2*Chem->GrCharge+1);

    if (Chem->Species[i].charge < 0)
    {
      de = Evln->NumDen[i]*Chem->Species[i].charge;
      negcharge[j] += de;
      negchargetot += de;
    }
  }

  /* charge makeup ratio */
  ratio = -dne/negchargetot;

  if (ratio > 1.0)
    return -1;
  else
  {
    for (i=Chem->GrInd; i<Chem->ManInd; i++)
    {
      if (Chem->Species[i].charge < 0)
      {
        Evln->NumDen[i] *= (1.0-ratio);
      }
      if (Chem->Species[i].charge == 0)
      {
        /* get the index of grain type */
        j = (i-Chem->GrInd)/(2*Chem->GrCharge+1);
        Evln->NumDen[i] += ratio*negcharge[j];
      }
    }
  }

  free_1d_array(negcharge);

  return 0;
}

/*---------------------------------------------------------------------------*/
/* Copy array
 */
void cparray(Real *ary1, Real *ary2, int N)
{
  int i;
  for (i=0; i<N; i++)
    ary2[i]=ary1[i];

  return;
}

/* convert index from neu-man to full sp list */
int ConvertInd (int i)
{
  int p;
  if (Chem.NGrain>0)
  {
    if (i<Chem.N_Neu_f)
      p = i+1;
    else if (i<(Chem.N_Neu_f+Chem.N_Neu))
      p = i+1+2*Chem.N_Neu_f;
    else if (i<(Chem.N_Neu_f+Chem.N_Neu+Chem.N_Neu_s))
      p = i+1+2*Chem.N_Neu_f+Chem.N_Neu;
    else
      p = i-(Chem.N_Neu_f+Chem.N_Neu
          +Chem.N_Neu_s)+Chem.ManInd;
    return p;
  }else
  return i;
}

int IConvertInd (int p)
{
  int i;
  if (Chem.NGrain>0)
  {
    if (p>=Chem.ManInd)
      i = p-Chem.ManInd+(Chem.N_Neu_f+Chem.N_Neu+Chem.N_Neu_s); 
    else if (p>=Chem.SNeuInd)
      i = p-1-2*Chem.N_Neu_f-Chem.N_Neu;
    else if (p>=Chem.NeuInd)
      i = p-1-2*Chem.N_Neu_f;
    else
      i = p-1;
  }else{
  return p;
  }
}

int AdDesRate(ChemEvln *Evln)
{
  Chemistry *Chem = Evln->Chem;
  EquationTerm *EqTerm;
  int q,k,i;
  /* add ad/des reactions rate */
  if (Chem->NGrain>0){
    Chem->NNeuT = Chem->N_Neu_f+Chem->N_Neu+Chem->N_Neu_s;
    Chem->AdDesRate = (AdDesCoef*)calloc_1d_array(Chem->NNeuT,sizeof(AdDesCoef));

    for (k=0;k<Chem->NNeuT;k++){
      q = ConvertInd(k);
      for (i=0; i<Chem->Equations[q].NTerm; i++)
      {
        EqTerm = &(Chem->Equations[q].EqTerm[i]);
        if (EqTerm->type == 3)  // adsorption
           Chem->AdDesRate[k].ad= Evln->K[EqTerm->ind];
        if (EqTerm->type == 4)  // desorption
           Chem->AdDesRate[k].des = Evln->K[EqTerm->ind];
      }
     }
    return 1;
   }
   else
     return -1;
}

