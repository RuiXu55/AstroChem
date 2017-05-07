#include "../copyright.h"
/*=============================================================================
 * FILE: stifbs.c
 *
 * PURPOSE: Contains solvers of stiff set of differential equations stifbs.c
 *   adopted from Numerical Recipes, modified so that the array indices range
 *   from 0..N-1 (rather than 1..N), and with the chemistry structure inserted
 *   as a parameter.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  - stifbs()       - solve stiff set of ODEs
 *
 * REFERENCES:
 *   Press, W., Teukolsky, S., Vetterling, W., Flannery, B., 1992, Numerical
 *     Recipes in C, Chapter 16.6, Cambridge University Press
 *
==============================================================================*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../defs.h"
#include "chemistry.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef CHEMISTRY

#define NRANSI
#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

Real **d,*x;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   simpr()    - semi-implicit ODE integrator (2nd order)
 *   pzextr()   - extropolation and error estimate
 *============================================================================*/

void simpr(ChemEvln *Evln, Real *y, Real *dydx, Real *dfdx, Real **dfdy,
                          int n,Real xs, Real htot, int nstep, Real *yout);
void pzextr(int iest, Real xest, Real *yest, Real *yz, Real *dy, int nv);


/*============================================================================*/
/*---------------------------- Public Functions ------------------------------*/

/*----------------------------------------------------------------------------*/
/* Numerical Recipes routine stifb: use the semi-implicit extropolation method
 * to solve stiff set of ordinary differential equations
 * Input parameters:
 *   y:     number density;               dydx:  reaction rates;
 *   nv:    Ntot (total # of species);    xx:    time (s)
 *   htry:  trial time step (s);          eps:   error level (%)
 *   yscal: density variation scale;      hdid:  time step did
 *   hnext: suggested next time step
 */
int stifbs(ChemEvln *Evln, Real *y, Real *dydx, int nv, Real *xx,
                Real htry, Real eps, Real *yscal, Real *hdid, Real *hnext)
{
  int i,iq,k,kk,km;
  static int first=1,kmax,kopt,nvold = -1;
  static Real epsold = -1.0,xnew;
  Real eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
  Real *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
  static Real a[IMAXX];
  static Real alf[KMAXX][KMAXX];
  static int nseq[IMAXX]={2,6,10,14,22,34,50,70};
  int reduct,exitflag=0;

/* initialization */

  d    = (Real**)calloc_2d_array(nv,KMAXX,sizeof(Real));
  dfdx = (Real*) calloc_1d_array(nv, sizeof(Real));
  dfdy = (Real**)calloc_2d_array(nv,nv, sizeof(Real));
  err  = (Real*) calloc_1d_array(KMAXX, sizeof(Real));
  x    = (Real*) calloc_1d_array(KMAXX, sizeof(Real));
  yerr = (Real*) calloc_1d_array(nv, sizeof(Real));
  ysav = (Real*) calloc_1d_array(nv, sizeof(Real));
  yseq = (Real*) calloc_1d_array(nv, sizeof(Real));

  if(eps != epsold || nv != nvold)
  {
    *hnext = xnew = -1.0e29;

    eps1 = SAFE1*eps;

    a[0] = nseq[0]+1;
    for (k=0;k<KMAXX;k++) {
      a[k+1]=a[k]+nseq[k+1];
    }

    for (iq=1;iq<KMAXX;iq++)
    {
      for (k=0;k<iq;k++)
        alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/((a[iq+1]-a[0]+1.0)*(2*k+3))));
    }

    epsold=eps;
    nvold=nv;
    a[0] += nv;

    for (k=0;k<KMAXX;k++)
      a[k+1]=a[k]+nseq[k+1];

    for (kopt=1;kopt<KMAXX-1;kopt++)
      if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;

    kmax=kopt;
  }

  h=htry;

  for (i=0;i<nv;i++)
    ysav[i]=y[i];

  jacobi(Evln, y, dfdy);

  if (*xx != xnew || h != (*hnext)) {
    first=1;
    kopt=kmax;
  }

  reduct=0;

  for (;;)
  {/* outer loop: adjust the total time step h */

    for (k=0;k<=kmax;k++)
    { /* inner loop: the step size sequence */

      xnew=(*xx)+h;
      if (xnew == (*xx))
      {
        ath_perr(0,"[stifbs]: step size underflow in stifbs\n");
        return -1;
      }

      simpr(Evln,ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq);
      xest=SQR(h/nseq[k]);

      pzextr(k,xest,yseq,y,yerr,nv);

      if (k != 0) {
        errmax=TINY;

        for (i=0;i<nv;i++) {
          errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
        }

        errmax /= eps;

        km=k-1;
        err[km]=pow(errmax/SAFE1,1.0/(2*km+3));
      }

      if (k != 0 && (k >= kopt-1 || first)) {
        if (errmax < 1.0) {
          exitflag=1;
          break;
        }
        if (k == kmax || k == kopt+1) {
          red=SAFE2/err[km];
          break;
        }
        else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
          red=1.0/err[km];
          break;
        }
        else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
          red=alf[km][kmax-1]*SAFE2/err[km];
          break;
        }
        else if (alf[km][kopt] < err[km]) {
          red=alf[km][kopt-1]/err[km];
          break;
        }
      }
    }

    if (exitflag) break;
    red=MIN(red,REDMIN);
    red=MAX(red,REDMAX);
    h *= red;
    reduct=1;
  }

  *xx=xnew;
  *hdid=h;
  first=0;
  wrkmin=1.0e35;

  for (kk=0;kk<=km;kk++) {
    fact=MAX(err[kk],SCALMX);
    work=fact*a[kk+1];

    if (work < wrkmin) {
      scale=fact;
      wrkmin=work;
      kopt=kk+1;
    }
  }

  *hnext=h/scale;

  if (kopt >= k && kopt != kmax && !reduct) {
    fact=MAX(scale/alf[kopt-1][kopt],SCALMX);
    if (a[kopt+1]*fact <= wrkmin) {
      *hnext=h/fact;
      kopt++;
    }
  }

  free_1d_array(yseq);
  free_1d_array(ysav);
  free_1d_array(yerr);
  free_1d_array(x);
  free_1d_array(err);
  free_2d_array(dfdy);
  free_2d_array(dfdx);
  free_2d_array(d);

  return 0;
}

/*============================================================================*/
/*---------------------------- Private Functions -----------------------------*/

/*----------------------------------------------------------------------------*/
/* Semi-implicit method to integrate stiff ordinary differential equations
 */
void simpr(ChemEvln *Evln, Real *y, Real *dydx, Real *dfdx, Real **dfdy,
                           int n, Real xs, Real htot, int nstep, Real *yout)
{
  int i,j,nn,*indx;
  Real d,h,x,**a,*del,*ytemp;

  indx=(int*)calloc_1d_array(n,sizeof(int));
  a=(Real**)calloc_2d_array(n,n,sizeof(Real));
  del=(Real*)calloc_1d_array(n,sizeof(Real));
  ytemp=(Real*)calloc_1d_array(n,sizeof(Real));

  h=htot/nstep;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) a[i][j] = -h*dfdy[i][j];
      ++a[i][i];
  }

  ludcmp(a,n,indx,&d);

  for (i=0;i<n;i++)
    yout[i]=h*(dydx[i]+h*dfdx[i]);

  lubksb(a,n,indx,yout);

  for (i=0;i<n;i++)
    ytemp[i]=y[i]+(del[i]=yout[i]);

  x=xs+h;

  /* user defined function to calculate the derivatives */
  derivs(Evln,ytemp,yout);

  for (nn=2;nn<=nstep;nn++)
  {
    for (i=0;i<n;i++)
        yout[i]=h*yout[i]-del[i];

    lubksb(a,n,indx,yout);
    for (i=0;i<n;i++)
      ytemp[i] += (del[i] += 2.0*yout[i]);

    x += h;
    (*derivs)(Evln,ytemp,yout);
  }

  for (i=0;i<n;i++)
    yout[i]=h*yout[i]-del[i];

  lubksb(a,n,indx,yout);

  for (i=0;i<n;i++)
    yout[i] += ytemp[i];

  free_1d_array(ytemp);
  free_1d_array(del);
  free_2d_array(a);
  free_1d_array(indx);

  return;
}

/*----------------------------------------------------------------------------*/
/* Extropolation and error estimate
 */
void pzextr(int iest, Real xest, Real *yest, Real *yz, Real *dy, int nv)
{
  int k1,j;
  Real q,f2,f1,delta,*c;

  c=(Real*)calloc_1d_array(nv, sizeof(Real));

  x[iest]=xest;

  for (j=0;j<nv;j++)
    dy[j]=yz[j]=yest[j];

  if (iest == 0) {
    for (j=0;j<nv;j++)
      d[j][0]=yest[j];
  }
  else
  {
    for (j=0;j<nv;j++)
      c[j]=yest[j];

    for (k1=0;k1<iest;k1++)
    {
      delta=1.0/(x[iest-k1-1]-xest);
      f1=xest*delta;
      f2=x[iest-k1-1]*delta;

      for (j=0;j<nv;j++) {
        q=d[j][k1];
        d[j][k1]=dy[j];
        delta=c[j]-q;
        dy[j]=f1*delta;
        c[j]=f2*delta;
        yz[j] += dy[j];
      }
    }
    for (j=0;j<nv;j++)
      d[j][iest]=dy[j];
  }

  free_1d_array(c);

  return;
}

#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
#undef NRANSI

#endif
