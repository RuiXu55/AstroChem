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

/* inverse matrix */
void inverse(Real* A, int N)
{
  int IPIV[N];
  int LWORK = N*N;
  Real WORK[LWORK];
  int INFO;
  dgetrf_(&N, &N, A, &N, IPIV, &INFO);
  dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);
}

/*============================================================================*/
int evolve(Real tend, Real dttry, Real abstol)
{

  realtype t,t1,dt,reltol=1.e-6;
  Real sum,rate, Tn;
  int flag, status,verbose,i,j,k,k1,p,p1,q, Ns,Nsp,Nsp1;
  N_Vector numden,dndt;
  Real *numden1;
  EquationTerm *EqTerm;
  void* cvode_mem;
  Chemistry *Chem = Evln.Chem;
  cvode_mem = NULL;
  Real ad, des,tmp,tmp1;
  AdDes ds;


  // Ns=# of Neutral species, Nsp is non-mantle species
  Ns = Chem->N_Neu_f + Chem->N_Neu+Chem->N_Neu_s;
  Nsp = Chem->Ntot-Ns*Chem->NGrain;

  dndt = N_VNew_Serial(Nsp);
  numden = N_VNew_Serial(Nsp);

  /* initialize number density for calculation */
  for(i=0;i<Nsp;i++){
    NV_Ith_S(numden,i) = Evln.NumDen[i]; 
    NV_Ith_S(dndt,i) = 0.0;
  }
  /* init CVode */ 
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeInit(cvode_mem,f,0.0,numden);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  flag = CVSpgmr(cvode_mem,PREC_LEFT,Nsp);
  flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
  flag = CVBandPrecInit(cvode_mem,Nsp,Nsp,Nsp); //N, mu, ml
  flag = CVodeSetMaxNumSteps(cvode_mem, 500000);

  // Mantle and their neutral counterpart
  if (Chem->NGrain>0){
    Nsp1    = (Chem->NGrain+1)*Ns;
    numden1 = (Real*)calloc_1d_array(Nsp1,sizeof(Real));

    ds.inv = (AdDesCoef*)calloc_1d_array(Ns,sizeof(AdDesCoef));
    AdDesCoef dat, invd;

    for (k=0;k<Ns;k++){
      q = ConvertInd(k);
      for (i=0; i<Chem->Equations[q].NTerm; i++)
      {
        EqTerm = &(Chem->Equations[q].EqTerm[i]);
        if (EqTerm->type == 3)  // adsorption
           ad= Evln.K[EqTerm->ind];
        if (EqTerm->type == 4)  // desorption
           des = Evln.K[EqTerm->ind];
      }
      ds.inv[k].d1 = ad;
      ds.inv[k].d2 = des*0;
      /*
      tmp  = 1.+ ad*dt-ad*des*dt*dt/(1.+des*dt);
      tmp1 = 1.+des*dt-ad*des*dt*dt/(1.+ad*dt);
      ds.inv[k].d1 = 1./tmp;
      ds.inv[k].d2 = des*dt/(1.+des*dt)/tmp;
      ds.inv[k].d3 = ad*dt/(1.+ad*dt)/tmp1;
      ds.inv[k].d4 = 1./tmp1;
      */
      /*
      dat.d1 = 1.0 + dt*ad;
      dat.d2 = -ad*dt;
      dat.d3 = -des*dt;
      dat.d4 = 1.0+ dt*des;

      Real coef = 1.0/(dat.d1*dat.d4-dat.d2*dat.d3);
      ds.inv[k].d1 = coef*dat.d4;
      ds.inv[k].d2 = -coef*dat.d2;
      ds.inv[k].d3 = -coef*dat.d3;
      ds.inv[k].d4 = coef*dat.d1;
      */ 
    }
  }


  clock_t c0, c1; /* Timing the code */
  c0 = clock();
  verbose = 0;
  Evln.t = dttry;
  dt = dttry; // year

  ath_pout(0,"\n Chemical evolution started...\n");

  while(Evln.t<tend)
  {
    //coeff_adj(&Evln);
    for(i=0;i<Nsp;i++)
      NV_Ith_S(numden,i) = Evln.NumDen[i];
    //flag = CVode(cvode_mem, Evln.t, numden, &t, CV_NORMAL);
    for(i=0;i<Nsp;i++)
      Evln.NumDen[i] = NV_Ith_S(numden,i);

    //for(i=0;i<Chem->Ntot;i++)
     // ath_pout(0,"A:sp=%s, #=%10e\n",Chem->Species[i].name,Evln.NumDen[i]);

    // update ad/des separately 
    if (Chem->NGrain>0)
    {
      for (k=0;k<Ns;k++){
	p   = ConvertInd(k);
	p1  = ConvertInd(k+Ns);
	Tn  = Evln.NumDen[p]+Evln.NumDen[p1];
	numden1[k]  = (Evln.NumDen[p]-Tn*ds.inv[k].d2
          /(ds.inv[k].d1+ds.inv[k].d2))*
	  exp(-(ds.inv[k].d1+ds.inv[k].d2)*dt)+
	  Tn*ds.inv[k].d2/(ds.inv[k].d1+ds.inv[k].d2);
	numden1[k+Ns] = Tn-numden1[k];
	ath_pout(0,"sp1=%s,sp2=%s\n",Chem->Species[p].name,Chem->Species[p1].name);
	ath_pout(0,"Tn=%10e, #1=%10e, #2=%10e\n",Tn,numden1[k],numden1[k+Ns]);
      }

      for (k=0;k<Nsp1;k++){
        p = ConvertInd(k);
        Evln.NumDen[p] = numden1[k];
	ath_pout(0,"sp=%s\n",Chem->Species[p].name);
      }
    }
    dt = MIN(1.5*Evln.t, 1e2*OneYear+Evln.t)-Evln.t;
    Evln.t = MIN(1.5*Evln.t, 1e2*OneYear+Evln.t);
    //for(i=0;i<Chem->Ntot;i++)
    // ath_pout(0,"B:sp=%s, #=%10e\n",Chem->Species[i].name,Evln.NumDen[i]);

    status = EleMakeup(verbose);

    ath_pout(0,"evolution time (yr) = %e\n",Evln.t/OneYear);

    /* ends if evolution time is too large */
    c1 = clock();
    if (((c1-c0)/CLOCKS_PER_SEC > 300.0) || (status < 0))
      break;
  }

  Real *EleNumDen;        /* Number density of each element */
  /* Initialization */
  EleNumDen = (Real*)calloc_1d_array((Chem->N_Ele+Chem->NGrain), sizeof(Real));

  for (i=0; i<Chem->N_Ele + Chem->NGrain; i++) {
    EleNumDen[i] = 0.0;
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
  for (i=0; i<Chem->N_Ele + Chem->NGrain; i++) {
    ath_pout(0,"sp=%s, #=%10e\n",Chem->Elements[i].name, EleNumDen[i]);
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
  int i, j, k, p,Nsp;
  Real sum,rate;
  EquationTerm *EqTerm;
  Nsp = Chem.Ntot-(Chem.N_Neu_f +
    Chem.N_Neu + Chem.N_Neu_s)*Chem.NGrain;
  for (k=0; k<Nsp; k++)
  {
    sum  = 0.0;
    for (i=0; i<Chem.Equations[k].NTerm; i++)
    {
      EqTerm = &(Chem.Equations[k].EqTerm[i]);
      if ((EqTerm->type) != 3 && (EqTerm->type !=4)){
        rate = Evln.K[EqTerm->ind] * EqTerm->dir;
        for (j=0; j<EqTerm->N; j++)
        {
          p = EqTerm->lab[j];
          rate *= NV_Ith_S(numden,p);
        }
        sum += rate;
      }
    }
   NV_Ith_S(dndt,k) = sum;
  }
  return(0);
}

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
      Evln.NumDen[i] = 0.0;
    }
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
         den += Evln.NumDen[l]*Chem->Species[l].composition[i];
       }
       frac = disp/MAX(TINY_NUMBER,den);

       for (j=0; j<Chem->Elements[i].numsig; j++)
       {
         l = Chem->Elements[i].single[j];
         Evln.NumDen[l] *= (1.0 - frac);
       }
     }

     /* if abundance is larger than the true value, then reduce
      * the number densities of all neutral species containing this element
      */
     else
     {
       status = EleMakeup_sub(i, disp);
     }
      
     if (status < 0)
       return status;
   } // loop for i 0-Chem.Ele

  /* Make up for the grain densities */
  for (i=Chem->N_Ele; i<Chem->N_Ele+Chem->NGrain; i++)
  {
    /* Calculate the discrepency */
    disp = (EleNumDen[i] - Chem->Elements[i].abundance/Evln.Abn_Den);

    ath_pout(verbose,"Discrepancy for %3s : %e over %e\n",
      Chem->Elements[i].name, disp, Chem->Elements[i].abundance/Evln.Abn_Den);

    /* Density make up */
    frac = disp / MAX(TINY_NUMBER,EleNumDen[i]);

    //for (j=Chem->GrInd; j<Chem->Ntot; j++)
    for (j=Chem->GrInd; j<Chem->ManInd; j++)
    {
      if (Chem->Species[j].composition[i] > 0)
        Evln.NumDen[j] *= (1.0-frac);
    }
  }

  /* Calculate the charge density */
  for (i=0; i<Chem->Ntot; i++)
  {
    if ((i != 0) && (Chem->Species[i].charge != 0))
      ChargeDen += Evln.NumDen[i] * Chem->Species[i].charge;
  }

  /* Make up for the charge density */

  if (ChargeDen >= 0.0)
  {
    Evln.NumDen[0] = ChargeDen;
  }
  else
  {
    Evln.NumDen[0] = 0.0;
    status = ChargeMakeup( -ChargeDen);
  }

  free_1d_array(EleNumDen);

  return status;
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

  frac = dn / MAX(TINY_NUMBER,den);

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
  //for (i=Chem->GrInd; i<Chem->Ntot; i++)
  for (i=Chem->GrInd; i<Chem->ManInd; i++)
  {
    /* get the index of grain type */
    j = (i-Chem->GrInd)/(2*Chem->GrCharge+1);
    if(j>=Chem->NGrain)
      ath_perr(0,"Grain type index error in evolve.c!\n");

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
    for (i=Chem->GrInd; i<Chem->ManInd; i++)
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

