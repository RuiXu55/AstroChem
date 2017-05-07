/*=============================================================================
 * FILE: evolve.c
 * PURPOSE: Contains functions to evolve the number densities of all species
 *   based on the chemistry model.
 * CONTAINS PUBLIC FUNCTIONS:
 *   EleMakeup() - density makeup for charge/element conservation
 * REFERENCES:
 *   Bai, X.-N. & Goodman, J., 2009, ApJ, 701, 737
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
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_spgmr.h>

int EleMakeup( int verbose);
int EleMakeup_sub(int q, Real dn);
int ChargeMakeup ( Real dne);

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*============================================================================*/
int evolve(Real tend, Real dttry, Real abstol)
{
  realtype t,reltol=1.e-5;
  int flag, status,verbose,i;
  N_Vector numden,dndt,vrtol;
  void* cvode_mem;
  Chemistry *Chem = Evln.Chem;
  numden = cvode_mem = NULL;

  dndt = N_VNew_Serial(Chem->Ntot);
  numden = N_VNew_Serial(Chem->Ntot);
  vrtol = N_VNew_Serial(Chem->Ntot);

  /* initialize number density for calculation */
  for(i=0;i<Chem->Ntot;i++){
    NV_Ith_S(numden,i) = Evln.NumDen[i]; 
    NV_Ith_S(dndt,i) = 0.0;
    NV_Ith_S(vrtol,i) = 1.e-2/Evln.DenScale[i];
  }

  /* init CVode */ 
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeInit(cvode_mem,f,0.0,numden);
  if(check_flag(&flag,"CVodeInit", 1)) return(1);

  //flag = CVodeSVtolerances(cvode_mem, abstol, vrtol);
  //if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1); 
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  //flag = CVSpgmr(cvode_mem,PREC_LEFT,Chem->Ntot);
  flag = CVSpgmr(cvode_mem,PREC_RIGHT,Chem->Ntot);
  flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
    if(check_flag(&flag, "CVSpilsSetGSType", 1)) return(1);
  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  //flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
  //if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
  flag = CVBandPrecInit(cvode_mem,Chem->Ntot,Chem->Ntot,Chem->Ntot); //N, mu, ml
  if(check_flag(&flag,"CVBandPrecInit", 0)) return(1); 
  flag = CVodeSetMaxNumSteps(cvode_mem, 500000);

  clock_t c0, c1; /* Timing the code */
  c0 = clock();
  ath_pout(0,"\n Chemical evolution started...\n");
  verbose = 0;
  Evln.t = dttry;
  while(Evln.t<tend)
  {
    //coeff_adj(&Evln);
    /* copy species number density to cvode to evolve */
    for(i=0;i<Chem->Ntot;i++)
      NV_Ith_S(numden,i) = Evln.NumDen[i];
    flag = CVode(cvode_mem,Evln.t, numden, &t, CV_NORMAL);
    Evln.t *= 1.2;
    //Evln.t = MIN(1.2*Evln.t, 1e4*OneYear+Evln.t);

    /* copy species # density back and impose conservation */
    for(i=0;i<Chem->Ntot;i++)
      Evln.NumDen[i] = NV_Ith_S(numden,i);
    ath_pout(0,"evolution time (yr) = %e\n",Evln.t/OneYear);
    status = EleMakeup(verbose);

    /* ends if evolution time is too large */
    c1 = clock();
    if (((c1-c0)/CLOCKS_PER_SEC > 300.0) || (status < 0))
      break;
  }

  /* finalize and return the status */
  N_VDestroy_Serial(dndt);
  N_VDestroy_Serial(numden);
  CVodeFree(&cvode_mem);
  ath_pout(0,"Evolution completed at t=%e yr, with Abn(e-)=%e.\n",
     Evln.t/OneYear, Evln.NumDen[0]*Evln.Abn_Den);
  return(0);
}

/*----------------------------------------------------------------------------*/
/* user provided routine for calculating the Jacobi matrix
 */

static int f(realtype t, N_Vector numden, N_Vector dndt, void *user_data)
{
  int i, j, k, p;
  Real sum,rate;
  EquationTerm *EqTerm;
  for (k=0; k<Chem.Ntot; k++)
  {
    sum  = 0.0;
    for (i=0; i<Chem.Equations[k].NTerm; i++)
    {
      EqTerm = &(Chem.Equations[k].EqTerm[i]);
      rate = Evln.K[EqTerm->ind] * EqTerm->dir;
      for (j=0; j<EqTerm->N; j++)
      {
        p = EqTerm->lab[j];
        rate *= NV_Ith_S(numden,p);
      }
      sum += rate;
    }
   NV_Ith_S(dndt,k) = sum;
  }
  return(0);
}

/*---------------------------------------------------------------------------*/
/* Make up the element number density to enforce conservation laws
 */
int EleMakeup(int verbose)
{
  int i, j, k, l, status=0;
  Real den, denmax, disp, frac;

  Real *EleNumDen;        /* Number density of each element */
  Real ChargeDen = 0.0;   /* Total charge number density excluding electron */

  Chemistry *Chem = Evln.Chem;
  Real    *NumDen = Evln.NumDen;

  /* Initialization */
  EleNumDen = (Real*)calloc_1d_array((Chem->N_Ele+Chem->NGrain), sizeof(Real));

  for (i=0; i<Chem->N_Ele + Chem->NGrain; i++) {
    EleNumDen[i] = 0.0;
  }

  /* For those with negative density, set them to zero */
  for (i=0; i<Chem->Ntot; i++) {
    if (Evln.NumDen[i]< 0.0)
    {
      ath_pout(verbose, "Warning: At t=%e yr, [%s] = %e < 0!\n",
                            Evln.t/OneYear, Chem->Species[i].name,Evln.NumDen[i]);
      NumDen[i] = 0.0;
  }

  /* Calculate the elemental density */
  for (i=0; i<Chem->Ntot; i++)
  {
    for (j=0; j<Chem->N_Ele+Chem->NGrain; j++)
    {
     if (Chem->Species[i].composition[j] > 0)
        EleNumDen[j] += Evln.NumDen[i]*Chem->Species[i].composition[j];
    }
  }

   /* Make up for the element densities */
   for (i=0; i<Chem->N_Ele; i++)
   {
     /* Calculate the discrepency */
     disp = (EleNumDen[i] - Chem->Elements[i].abundance/Evln.Abn_Den);
     ath_pout(verbose,"Discrepancy for %3s : %e over %e\n",
     Chem->Elements[i].name, disp, Chem->Elements[i].abundance/Evln.Abn_Den);

    /* if abundance is smaller than the true value, then increase
     * the number densities of its single-element species
     */
     if (disp<0.0){
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
        if( NumDen[l]<0.0) NumDen[l] = 0.0;
      }
     }
     /* if abundance is larger than the true value, then reduce
      * 307      * the number densities of all neutral species containing this element
      * 308      */
     else
     {
       status = EleMakeup_sub(i, disp);
      }
      
      if (status < 0)
        return status;
      }
/* Make up for the grain densities */
  for (i=Chem->N_Ele; i<Chem->N_Ele+Chem->NGrain; i++)
  {
    /* Calculate the discrepency */
    disp = (EleNumDen[i] - Chem->Elements[i].abundance/Evln.Abn_Den);

    ath_pout(verbose,"Discrepancy for %3s : %e over %e\n",
      Chem->Elements[i].name, disp, Chem->Elements[i].abundance/Evln.Abn_Den);

    /* Density make up */
    frac = disp / EleNumDen[i];

    for (j=Chem->GrInd; j<Chem->Ntot; j++)
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

    status = ChargeMakeup( -ChargeDen);
  }

  free_1d_array(EleNumDen);

  return status;
}
}

/*---------------------------------------------------------------------------*/
/* Density make up for single element species
 * q:  name of the element
 * dn: density makeup
 */
int EleMakeup_sub(int q, Real dn)
{
  int i, j, k=0;
  Real den, frac, dni;
  Chemistry *Chem = Evln.Chem;
  SpeciesInfo *Species = Chem->Species;

/* Find the neutral species containing this element
 * and have the largest number density
 */
  den = 0.0;

  for (i=0; i<Chem->Ntot; i++)
  {
    if ((Species[i].composition[q] > 0) && (Species[i].charge == 0))
    {
      den += Evln.NumDen[i] * Species[i].composition[q];
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
      dni = Evln.NumDen[i] * frac;

      Evln.NumDen[i] *= (1.0-frac);

      for (j=0; j<Chem->N_Ele+Chem->NGrain; j++)
      {
        if (Species[i].composition[j] > 0)
        {
          if (j != q)
          {
            k = Chem->Elements[j].single[0];
            Evln.NumDen[k] +=
                  dni*Species[i].composition[j]/Species[k].composition[j];
          }
        }
      }
    }
  }

  return 0;
}


/*============================================================================*/
/*------------------------------ PRIVATE FUNCTIONS ---------------------------*/
/*---------------------------------------------------------------------------*/
/* Density make up for the electrons
 */
int ChargeMakeup(Real dne)
{
  int i, j;
  Real ratio, negchargetot, de;
  Real *negcharge;
  Chemistry *Chem= Evln.Chem;

  negchargetot = 0.0;
  negcharge = (Real*)calloc(Chem->NGrain, sizeof(Real));
  for (j=0; j<Chem->NGrain; j++)
    negcharge[j] = 0.0;

  /* calculate the total negative charge */
  for (i=Chem->GrInd; i<Chem->Ntot; i++)
  {
    /* get the index of grain type */
    j = (i-Chem->GrInd)/(2*Chem->GrCharge+1);

    if (Chem->Species[i].charge < 0)
    {
      de = Evln.NumDen[i]*Chem->Species[i].charge;
      negcharge[j] += de;
      negchargetot += de;
    }
  }

  /* charge makeup ratio */
  ratio = -dne/negchargetot;

  if (ratio > 1.0){
    ath_perr(0,"can't makeup charge density!\n");
    return -1;
    exit(0);
  }
  else
  {
    for (i=Chem->GrInd; i<Chem->Ntot; i++)
    {
      if (Chem->Species[i].charge < 0)
      {
        Evln.NumDen[i] *= (1.0-ratio);
      }
      if (Chem->Species[i].charge == 0)
      {
        /* get the index of grain type */
        j = (i-Chem->GrInd)/(2*Chem->GrCharge+1);
        Evln.NumDen[i] += ratio*negcharge[j];
      }
    }
  }

  free_1d_array(negcharge);

  return 0;
}

/* Check flag for CVode Setup */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }
  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }
  return(0);
}

