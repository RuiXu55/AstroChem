#include "../header/copyright.h"
/*=============================================================================
 * FILE: coeff.c
 *
 * PURPOSE: Calculate the rate coefficients for all reactions. There are 6
 *   types of reactions labeled by "rtype", ranging from 0 to 5:
 *
 *     0: ionization reactions;
 *     1: gas-phase reactions;
 *     2: ion/electron + grain reactions;
 *     3: neutral-grain reactions (adsorption);
 *     4: desorption reactions;
 *     5: grain-grain reactions;
 *
 *   Different types of reactions are calculated by different functions.
 *
 *   We use cgs unit for all rate coefficients.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  - IonizationCoeff() - set the ionization rate coefficient
 *  - CalCoeff()        - calculate rate coefficients for all other reactions
 *  - ChemCoeff()       - get the rate coefficients for gas-phase reactions
 *  - IonGrCoeff()      - get the rate coefficients for ion-grain reactions
 *  - NeuGrCoeff()      - get the rate coefficients for neutral-grain reactions
 *  - DesorpCoeff()     - get the rate coefficients for desorption reactions
 *  - GrGrCoeff()       - get the rate coefficients for grain-grain reactions
 *  - EleStickCoeff()   - electron-grain sticking coefficient
 *
 * REFERENCES:
 *  Bai, X.-N. & Goodman, J., 2009, ApJ, 701, 737
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

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   NeuStickCoeff() - neutral-grain sticking coefficient
 *   PAdsorb()       - electron-grain absoption probability
 *   CSfun()         - electron-grain collision cross section
 *   Rootx()         - rootfinding for cross section/sticking coeff calculation
 *============================================================================*/

Real NeuStickCoeff(Real D, Real T0);
Real PAdsorb(Real size, Real epsi, Real nu);
Real PEscape(Real Delta, Real epsi, Real nu);
Real CSfun(Real epsi, Real nu);
Real Rootx(Real epsi, Real nu);

/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*------------------------------------------------------------------------------
 * Set the ionization rate coefficient
 */
void IonizationCoeff(ChemEvln *Evln, Real zeta_eff, Real Av, int verbose)
{
  int i;
  Chemistry *Chem = Evln->Chem;

  Evln->zeta_eff = zeta_eff;

  for (i=0; i<Chem->NReaction; i++)
  {
    if (Chem->Reactions[i].rtype == 0)
    { /* ionization reaction */
      Evln->K[i] = zeta_eff * Chem->Reactions[i].coeff[0].gamma;

      if (verbose == 0) {
        PrintReaction(Chem,i,Evln->K[i]);
      }
    }
    if (Chem->Reactions[i].rtype == 10)
    { /* photoionization reaction */
      Evln->K[i] = 10000.0*Chem->Reactions[i].coeff[0].alpha
             *exp(-Chem->Reactions[i].coeff[0].gamma*Av);

      if (verbose == 0) {
        PrintReaction(Chem,i,Evln->K[i]);
      }
    }
  }

  return;
}

void IonizationCoeff1(ChemEvln *Evln, Real zeta_eff, Real Av, Real G,int verbose)
{
  int i;
  Chemistry *Chem = Evln->Chem;

  Evln->zeta_eff = zeta_eff;

  for (i=0; i<Chem->NReaction; i++)
  {
    if (Chem->Reactions[i].rtype == 0)
    { /* ionization reaction */
      Evln->K[i] = zeta_eff * Chem->Reactions[i].coeff[0].gamma;

      if (verbose == 0) {
        PrintReaction(Chem,i,Evln->K[i]);
      }
    }
    if (Chem->Reactions[i].rtype == 10)
    { /* photoionization reaction */
      Evln->K[i] = G*Chem->Reactions[i].coeff[0].alpha
             *exp(-Chem->Reactions[i].coeff[0].gamma*Av);

      if (verbose == 0) {
        PrintReaction(Chem,i,Evln->K[i]);
      }
    }
  }

  return;
}
/*------------------------------------------------------------------------------
 * Calculate all other rate coefficients
 */
void CalCoeff(ChemEvln *Evln, Real T, int verbose)
{
  int i, type;
  Real K;
  Chemistry *Chem = Evln->Chem;

  Evln->T = T;

  /* Read coeffients from reaction */
  for (i=0; i<Chem->NReaction; i++)
  {
    type = Chem->Reactions[i].rtype;
    switch (type)
    {
      /* Ionization Reaction (K is set elsewhere) */
      case 0: break;

      /* Chemical Reaction */
      case 1: K = ChemCoeff(Chem->Reactions[i].coeff, T, 
                            Chem->Reactions[i].NumTRange);
              break;
      /* Photo-reactions */
      case 10: break;

      /* Charge+grain Reaction */
      case 2: K = IonGrCoeff(Chem, i, T);
              break;
      /* Neutral+grain */
      case 3: K = NeuGrCoeff(Chem, Evln, i);
              break;
      /* Desorption */
      case 4: K = DesorpCoeff(Chem->Reactions[i].coeff, T);
              break;
      /* Grain+Grain reaction */
      case 5: K = GrGrCoeff(Chem, i, T);
              break;
      /* Grain-surface reaction */
      case 6: K = GrSurfCoeff(Evln, i, T);
              break;
      default : ath_error("[coefficients]: reaction type should be 0-6!\n");
    }

    Evln->rate_adj[i] = 1.0;

    if ((type != 0) && (type != 10))
    {
      Evln->K[i] = K;

      if (verbose == 0) { 
        PrintReaction(Chem,i,K);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------
 * Dynamical rate adjustment
 */
void coeff_adj(ChemEvln *Evln)
{
  int i,i1,i2;
  Chemistry *Chem=Evln->Chem;

  if (Chem->NGrain>0)
  {
    /* update the GrAvail array */
    GrAvailFac(Evln);

    /* Read coeffients from reaction */
    for (i=0; i<Chem->NReaction; i++)
    {
      /* Adjust rate coefficients for desorption */
      if (Chem->Reactions[i].rtype == 4)
      {
        i1 = Chem->Reactions[i].reactant[0];
        Evln->rate_adj[i] = Evln->GrAvail[i1];
      }
      /* Adjust rate coefficients for grain-surface reactions */
      if (Chem->Reactions[i].rtype == 6)
      {
        i1 = Chem->Reactions[i].reactant[0];
        i2 = Chem->Reactions[i].reactant[1];
        Evln->rate_adj[i] = Evln->GrAvail[i1]*Evln->GrAvail[i2];
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------
 * Rate for chemical reactions (type 1)
 */
Real ChemCoeff(Coefficient *coeff, Real T, int NumTRange)
{
  int     i = NumTRange - 1;
  Real coef = 0.0;
  Real coef1,coef2;

  if (T > coeff[i].Tmax)	/* T too high */
  {
    coef = coeff[i].alpha*pow(T/300.0, coeff[i].beta)*exp(-coeff[i].gamma/T);
    //coef = coeff[0].alpha * pow(coeff[0].Tmax/300.0, coeff[0].beta)
    //       * exp(-coeff[0].gamma/coeff[0].Tmax);
  }
  else if (T < coeff[0].Tmin)	/* T too low */
  {
    coef = coeff[0].alpha*pow(T/300.0, coeff[0].beta)*exp(-coeff[0].gamma/T);
    //coef = coeff[0].alpha * pow(coeff[0].Tmin/300.0, coeff[0].beta)
    //     * exp(-coeff[0].gamma/coeff[0].Tmin);
  }

  else{

    while (i > 0)
      if (T < coeff[i].Tmin) i--;
      else break;

    if (T <= coeff[i].Tmax)
      coef = coeff[i].alpha*pow(T/300.0, coeff[i].beta)*exp(-coeff[i].gamma/T);
    else if ((T-coeff[i].Tmax) < (coeff[i+1].Tmin-T))
      coef = coeff[i].alpha * pow(coeff[i].Tmax/300.0, coeff[i].beta)
                            * exp(-coeff[i].gamma/coeff[i].Tmax);
    else
      coef = coeff[i].alpha * pow(coeff[i+1].Tmin/300.0, coeff[i].beta)
                            * exp(-coeff[i].gamma/coeff[i+1].Tmin);
  }

  return coef;
}

/*-----------------------------------------------------------------------------
 * Rate coeffieicnt for ion + grain reactions (type 2)
 * Reference: Draine & Sutin, 1987, ApJ, 320, 803
 */
Real IonGrCoeff(Chemistry *Chem, int i, Real T)
{
  Real q;	/* sign of charge */
  int r1, r2, Z;
  Real tau, nu, Jt, temp1, temp2, size;
  Real s;	/* sticking coefficient */
  Coefficient *coeff = Chem->Reactions[i].coeff;

  r1 = Chem->Reactions[i].reactant[0];
  r2 = Chem->Reactions[i].reactant[1];

  q = Chem->Species[r1].charge;

  Z = Chem->Species[r2].charge;      /* charge of grain */
  size = Chem->Species[r2].gsize;    /* GrainSize */

  tau = 17.96*(T/300.0)*size;
  nu = (Real)(Z)/(Real)(q);

  /* Calculating the cross section */
  if (Z == 0)
    Jt = 1.0+sqrt(1.5708/tau);
  else if (nu < 0.0)
    Jt = (1.0-nu/tau)*(1.0+sqrt(2.0/(tau-2.0*nu)));
  else
  {
    temp1 = 1.0+sqrt(1.0/(4.0*tau+3.0*nu));
    temp2 = nu/(1.0+1.0/sqrt(nu));
    Jt = SQR(temp1)*exp(-temp2/tau);
  }

  /* calculating sticking coefficient */
  if (r1 != 0)  /* Ions */
    s = 1.0;
  else          /* electrons */
    s = EleStickCoeff(size, Z, T);
  return 7.893e-3 * s * sqrt(T/300.0/coeff[0].alpha)
                  * SQR(size) * Jt * coeff[0].gamma;
	
}

/*-----------------------------------------------------------------------------
 * Rate coeffieicnt for neutral + grain reactions (type 3)
 */
Real NeuGrCoeff(Chemistry *Chem, ChemEvln *Evln, int i)
{/* Note: grainlab is -k-1 where k is grainkind */

  int grainlab = -1-Chem->Reactions[i].reactant[1];

  Coefficient *coeff = Chem->Reactions[i].coeff;

  Real s = NeuStickCoeff(coeff[0].beta, Evln->T);
  Real size = Chem->GrSize[grainlab];
  Real ngr = Chem->Elements[Chem->N_Ele+grainlab].abundance/Evln->Abn_Den;
  return 8.57e-3 * s * ngr * sqrt(Evln->T/300.0/coeff[0].alpha) * SQR(size);
}

/*-----------------------------------------------------------------------------
 * Desorption reaction (type 4)
 */
Real DesorpCoeff(Coefficient *coeff, Real T)
{
  return 2.74e12*sqrt(coeff[0].beta/300.0/coeff[0].alpha)*exp(-coeff[0].beta/T);
}

/*-----------------------------------------------------------------------------
 * Grain+Grain reaction (type 5)
 */
Real GrGrCoeff(Chemistry *Chem, int i, Real T)
{
  int m, n, r1, r2;
  Real size, mu;
  
  r1 = Chem->Reactions[i].reactant[0];
  r2 = Chem->Reactions[i].reactant[1];

  m = Chem->Species[r1].charge;
  n = Chem->Species[r2].charge;
  size = Chem->Species[r1].gsize + Chem->Species[r2].gsize;

  mu = Chem->Species[r1].mass * Chem->Species[r2].mass
            / (Chem->Species[r1].mass+Chem->Species[r2].mass);

  return 7.89e-3*size*size*sqrt(T/300.0/mu)*(1-0.0556*m*n/(size*T/300.0));
}

/*-----------------------------------------------------------------------------
 *  Grain-surface reaction (type 6)
 *  Reference: Hasegawa, Herbst & Leung, 1992, ApJS, 82, 167
 */
Real GrSurfCoeff(ChemEvln *Evln, int i, Real T)
{
  int ind1,ind2,grind,k;
  Chemistry *Chem=Evln->Chem;
  Real Ea,Eb1,Eb2,grsize,gr_abn,Ns,ngr;
  Real m1,m2,mu,nu01,nu02,thop1,thop2,kappa,ratio;

  ind1 = Chem->Reactions[i].reactant[0];
  ind2 = Chem->Reactions[i].reactant[1];
  grind= (int)(Chem->Reactions[i].coeff[0].beta);

  Eb1 = Chem->Species[ind1].Eb; /* in units of Kelvin */
  Eb2 = Chem->Species[ind2].Eb;
  Ea  = Chem->Reactions[i].coeff[0].alpha;
  ratio = Chem->Reactions[i].coeff[0].gamma; /* branching ratio */

  grsize = Chem->GrSize[grind]; /* micron */
  k = Chem->N_Ele+grind;
  gr_abn = Chem->Elements[k].abundance;
  Ns     = 1.8e8*SQR(grsize); /* total number of surface sites */
  ngr    = gr_abn/Evln->Abn_Den;

  m1 = Chem->Species[ind1].mass; /* in units of proton mass */
  m2 = Chem->Species[ind1].mass;
  mu = m1*m2/(m1+m2); /* reduced mass */

  nu01 = 1.58e11*sqrt(Eb1/m1); /* oscillation frequency (s^(-1)) */
  nu02 = 1.58e11*sqrt(Eb2/m2);

  thop1 = exp(0.3*Eb1/T)/nu01; /* hopping time */
  thop2 = exp(0.3*Eb2/T)/nu02;

  kappa = exp(-0.41*(mu*Ea));

  /* effective reaction rate coefficient */
  return kappa*(thop1+thop2)/(thop1*thop2*Ns*ngr)*ratio;
}

/*-----------------------------------------------------------------------------
 *  Grain-surface reaction: (dynamical) adjustment factor
 */
void GrAvailFac(ChemEvln *Evln)
{
  int i,k,ind;
  Chemistry *Chem = Evln->Chem;
  Real grsize,Ns,Nmol,fac;

  /* calculate the total mantle species number density for each grain type */
  for (i=0;i<Chem->Ntot;i++)
  {
    /* if it is a mantle species, return grain index, otherwise negative */
    ind = -(int)(Chem->Species[i].gsize)-1;

    /* here use GrAvail as a temperory array to save data */
    if (ind>=0)
    {
      k = ind+Chem->N_Ele;

      grsize = Chem->GrSize[ind]; /* micron */
      Ns     = 1.8e8*SQR(grsize);   /* total number of surface sites */

      /* number of mantle molecules per grain */
      Nmol = Evln->NumDen[i]/(Chem->Elements[k].abundance/Evln->Abn_Den);

      Evln->GrAvail[i] = 1.0;

      fac = Nmol/Ns;
      if (fac > 1.0)
      {
        Evln->GrAvail[i] = 1.0/fac;
        ath_pout(0,"GrAvail of species %s is %e.\n",Chem->Species[i].name,1.0/fac);
      }
    }
  }
}

/*-----------------------------------------------------------------------------
 * Calculating the sticking coefficient of e-+gr(Z) at temperature T0
 */
Real EleStickCoeff(Real size, int Z, Real T0)
{
  int m, n, i, j;
  Real x[10], a[10];
  Real valP, wei, temp, v0, vmin, vmax, v, dv, epsi;

  /* Numerical integration from 0 to 10 (Gauss-Legendre method) */
  valP = 0.0;	wei = 0.0;

  dv = 0.05*pow(100.0/T0,0.75); /* step size */
  m  = (int)(4.0/dv);       /* number of loops */

  n = 5;	/* 5th order G-L method, 9th order accuracy */

  x[0] = -0.9061798459;	x[1] = -0.5384693101;
  x[2] = 0;		x[3] = -x[1];	x[4] = -x[0];
  a[0] = 0.2369268851;	a[1] = 0.4786286705;
  a[2] = 0.5688888889;	a[3] = a[1];	a[4] = a[0];

  /* find the point where cross section is marginally 0 */
  if (CSfun(dv*T0/16.71*size, (Real)(-Z)) > 0.0)
    v0 = 0.0;
  else {
    vmin = dv;
    vmax = (Real)(-Z)*16.71/size/T0;
    while ((vmax-vmin) > 0.05*dv)
    {
      v0   = 0.5*(vmin+vmax);
      temp = CSfun(v0*T0/16.71*size, (Real)(-Z));
      if (temp > 0.0)
        vmax = v0;
      else
        vmin = v0;
    }
    v0   = 0.5*(vmin+vmax);
  }

  /* start the calculation */
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
    {
      v = v0 + dv*i+0.5*dv*(x[j]+1.0);	/* v=T/T0 */

      epsi = v*T0/16.71*size;           /* epsi = E*a/e^2 */

      temp = a[j]*v*exp(-v+v0)*CSfun(epsi, (Real)(-Z));
      valP += temp*PAdsorb(size, epsi, (Real)(-Z));

      wei += temp;
    }
  }

  if (wei == 0.0) /* Also, valP = 0 */
    wei = 1.0;

	//return 0.6;
  return valP/wei;
}

/*============================================================================*/
/*----------------------------- Private Functions ----------------------------*/

/*-----------------------------------------------------------------------------
 * Calculating the sticking coefficient of X+gr(Z) at temperature T0
 */
Real NeuStickCoeff(Real D, Real T0)
{/* D: binding energy */
  int m, n, i, j;
  Real x[10], a[10];
  Real valP, wei, temp, v;

  /* Numerical integration from 0 to 10 (Gauss-Radau method) */
  valP = 0.0;	wei = 0.0;

  m = 10;	/* accuracy to exp(-10); */

  n = 5;	/* 5th order G-L method, 9th order accuracy */

  x[0] = -0.9061798459;	x[1] = -0.5384693101;
  x[2] = 0;		x[3] = -x[1];	x[4] = -x[0];
  a[0] = 0.2369268851;	a[1] = 0.4786286705;
  a[2] = 0.5688888889;	a[3] = a[1];	a[4] = a[0];

  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
    {
      v = i+0.5*(x[j]+1.0);	/* v=T/T0; */
      temp = a[j]*v*exp(-v);
      valP += temp*exp(-SQR(v*T0)/(D*46.4));
      wei += temp;
    }
  }

  if (wei == 0.0) /* Also, valP = 0 */
    wei = 1.0;

  return valP/wei;
}

/*-----------------------------------------------------------------------------
 * The probability of a electron being adsorbed by the grain (a=size) with
 * initial energy E
 */
Real PAdsorb(Real size, Real epsi, Real nu) /* Note: E in unit of Kelvin */
{
  Real D = 3.0*694.3*size;      /* Binding energy, 1 eV, in unit of e^2/a */
  Real mr = 1.0/1836.0/12.0;    /* mass ratio of e- to C (grain made of C) */
  Real dEu = 4.0*mr*(epsi-nu+D);/* upper limit of energy transfer */
  Real Td = 420.0/16.71*size;   /* T_Debye for grain lattice (graphite) */
  Real alpha = 1.5*dEu/Td;      /* probability of inelastic collision */
  Real delE = Td/18.0;          /* average energy transfer at each collision */
  Real Ei = epsi;               /* remaining energy after the ith collision */
  Real bi = 0.0;                /* escape probability */
  Real PE = 1.0;                /* absorption probability */

  if (alpha > 1.0)
    alpha = 1.0;  /* probability <= 1 */

  if (dEu <= 0.0)
    return 1.0;  /* energy is too small, can not escape! */

  while ((Ei > 0.0) && (bi < 1.0))
  {
    bi = PEscape(D, Ei, nu);
    PE *= (1.0-bi)/(1.0-bi+bi/alpha);
    Ei -= delE;
  }

  return PE;
}

/* Newton Iteration in the PAdsorb calculation */
Real PEscape(Real Delta, Real epsi, Real nu)
{
  Real x, sin2, E1, E2;

  x = Rootx(epsi, nu);

  /* calculate the critical angle */
  E1 = epsi*SQR(x) + 0.5/(SQR(x)-1.0) - nu*x;
  E2 = epsi- nu + Delta;

  if ((E1 <= 0.0) || (E2 <= 0))
    return 0.0; /* energy is too small to escape */

  sin2 = E1/E2;

  if (sin2 >= 1.0)
    return 1.0;
  else
    return 1.0-sqrt(1.0-sin2);
}

/*-----------------------------------------------------------------------------
 * cross section function (epsi=E*a/e^2, nu=Q(gr)/q(ion/e))
 */
Real CSfun(Real epsi, Real nu)
{
/* Find the root related to the extreme (max) of the potential, 
 * then calculate reduced cross section
 */
  Real x,beta2;

  x = Rootx(epsi, nu);

  beta2 = SQR(x)+(0.5/(SQR(x)-1.0)-nu*x)/(epsi+TINY_NUMBER);

  if (beta2<0.0)	beta2 = 0.0;

  return beta2;
}

/*----------------------------------------------------------------------------- 
 * Newton iteration to obtain the critical radius
 * Parameters:
 *   epsi - E*a/e^2
 *   nu   - Z/q
 *
 * Solves equation (B11) of Nishi et al. (1991) but with grain charge included.
 * The equation is equivalent to equation (2.10) of Draine & Sutin (1987).
 *
 * Needed for calculating the ion/electron-grain cross section as well as the
 * electron-grain sticking coefficient.
 */
Real Rootx(Real epsi, Real nu)
{
  Real x, res, resp;

  /* find the starting point for the iteration */
  x = 2.0;
  res = x/SQR(SQR(x)-1.0) - 2.0*x*epsi + nu;
  while (res < 0.0)
  {
    x = 1.0 + 0.5*(x-1.0);
    res = x/SQR(SQR(x)-1.0) - 2.0*x*epsi + nu;
  }

  /* perform Newton's iteration */
  while (res > 1.0e-8)
  {
    resp = -(3.0*SQR(x)+1.0)/pow((SQR(x)-1.0),3.0)-2.0*epsi;
    x -= res/resp;
    res = x/SQR(SQR(x)-1.0) - 2.0*x*epsi + nu;
  }

  return x;
}

#endif /* CHEMISTRY */
