/*=============================================================================
 * FILE: stifkr.c
 *
 * PURPOSE: Contains solvers of stiff set of differential equations stifkr.c
 *   adopted from Numerical Recipes, modified so that the array indices range
 *   from 0..N-1 (rather than 1..N), and with the chemistry structure inserted
 *   as a parameter. It uses the Kaps-Rantrop algorithm (see NR for details).
 *   This integrator is an alternative for the stifbs, and works better for
 *   large set of equations (e.g., the complex network).
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  - stifkr()       - solve stiff set of ODEs
 *
 * REFERENCES:
 *   Press, W., Teukolsky, S., Vetterling, W., Flannery, B., 1992, Numerical
 *     Recipes in C, Chapter 16.6, Cambridge University Press
 *
==============================================================================*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../header/prototypes.h"
#include "../header/defs.h"
#include "../header/chemistry.h"

#ifdef CHEMISTRY

#define SAFETY 0.9
#define GROW 1.5
#define PGROW -0.25
#define SHRNK 0.5
#define PSHRNK (-1.0/3.0)
#define ERRCON 0.1296
#define MAXTRY 500
#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)
#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)
#define TINY 1.0e-20;

int stifkr(ChemEvln *Evln, Real *y, Real *dydx, int n, Real *x,
                Real htry, Real eps, Real *yscal, Real *hdid, Real *hnext)
{
  int i,j,jtry,*indx;
  Real d,errmax,h,xsav,**a,*dfdx,**dfdy,*dysav,*err;
  Real *g1,*g2,*g3,*g4,*ysav;

  indx  = (int*)calloc_1d_array(n, sizeof(int));
  a     = (Real**)calloc_2d_array(n,n, sizeof(Real));
  dfdx  = (Real*) calloc_1d_array(n, sizeof(Real));
  dfdy  = (Real**)calloc_2d_array(n,n, sizeof(Real));
  dysav = (Real*) calloc_1d_array(n, sizeof(Real));
  err   = (Real*) calloc_1d_array(n, sizeof(Real));
  g1    = (Real*) calloc_1d_array(n, sizeof(Real));
  g2    = (Real*) calloc_1d_array(n, sizeof(Real));
  g3    = (Real*) calloc_1d_array(n, sizeof(Real));
  g4    = (Real*) calloc_1d_array(n, sizeof(Real));
  ysav  = (Real*) calloc_1d_array(n, sizeof(Real));

  xsav=(*x);
  for (i=0;i<n;i++) {
    ysav[i]=y[i];
    dysav[i]=dydx[i];
  }

  jacobi(Evln, ysav, dfdy);

  h=htry;

  for (jtry=0;jtry<MAXTRY;jtry++) {
    for (i=0;i<n;i++) {

      for (j=0;j<n;j++)
        a[i][j] = -dfdy[i][j];

      a[i][i] += 1.0/(GAM*h);
    }

    ludcmp(a,n,indx,&d);

    for (i=0;i<n;i++)
      g1[i]=dysav[i]+h*C1X*dfdx[i];

    lubksb(a,n,indx,g1);

    for (i=0;i<n;i++)
      y[i]=ysav[i]+A21*g1[i];

    *x=xsav+A2X*h;

    derivs(Evln, y, dydx);

    for (i=0;i<n;i++)
      g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;

    lubksb(a,n,indx,g2);

    for (i=0;i<n;i++)
      y[i]=ysav[i]+A31*g1[i]+A32*g2[i];

    *x=xsav+A3X*h;

    derivs(Evln, y, dydx);

    for (i=0;i<n;i++)
      g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;

    lubksb(a,n,indx,g3);

    for (i=0;i<n;i++)
      g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;

    lubksb(a,n,indx,g4);

    for (i=0;i<n;i++) {
      y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];

      err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
    }

    *x=xsav+h;

    if (*x == xsav){
      ath_perr(0, "[stifkr]: stepsize not significant in stifkr\n");

      return -1;
    }

    errmax=0.0;
    for (i=0;i<n;i++)
      errmax=MAX(errmax,fabs(err[i]/yscal[i]));

    errmax /= eps;
    if (errmax <= 1.0) {

      *hdid=h;
      *hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);

      free_1d_array(ysav);
      free_1d_array(g4);
      free_1d_array(g3);
      free_1d_array(g2);
      free_1d_array(g1);
      free_1d_array(err);
      free_1d_array(dysav);
      free_1d_array(dfdx);
      free_1d_array(indx);
      free_2d_array(dfdy);
      free_2d_array(a);

      return 0;
    }
    else
    {
       *hnext=SAFETY*h*pow(errmax,PSHRNK);
       h=(h >= 0.0 ? MAX(*hnext,SHRNK*h) : MIN(*hnext,SHRNK*h));
    }
  }

  *hdid=h;
  *hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);

  free_1d_array(ysav);
  free_1d_array(g4);
  free_1d_array(g3);
  free_1d_array(g2);
  free_1d_array(g1);
  free_1d_array(err);
  free_1d_array(dysav);
  free_1d_array(dfdx);
  free_1d_array(indx);
  free_2d_array(dfdy);
  free_2d_array(a);

  ath_perr(0,"[stifkr]: Maximum number of trial steps exceeded...\n");

  return -1;
}

#undef TINY
#undef SAFETY
#undef GROW
#undef PGROW
#undef SHRNK
#undef PSHRNK
#undef ERRCON
#undef MAXTRY
#undef GAM
#undef A21
#undef A31
#undef A32
#undef C21
#undef C31
#undef C32
#undef C41
#undef C42
#undef C43
#undef B1
#undef B2
#undef B3
#undef B4
#undef E1
#undef E2
#undef E3
#undef E4
#undef C1X
#undef C2X
#undef C3X
#undef C4X
#undef A2X
#undef A3X
#undef NRANSI

#endif /* CHEMISTRY */
