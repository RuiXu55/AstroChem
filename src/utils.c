/*============================================================================*/
/*! \file utils.c
 *  \brief A variety of useful utility functions.
 *
 * PURPOSE: A variety of useful utility functions.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - ath_strdup()     - not supplied by fancy ANSI C, but ok in C89 
 * - ath_gcd()        - computes greatest common divisor by Euler's method
 * - ath_big_endian() - run-time detection of endianism of the host cpu
 * - ath_bswap()      - fast byte swapping routine
 * - ath_error()      - fatal error routine
 * - minmax1()        - fast Min/Max for a 1d array using registers
 * - minmax2()        - fast Min/Max for a 2d array using registers
 * - minmax3()        - fast Min/Max for a 3d array using registers
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "header/defs.h"
#include "header/chemistry.h"
#include "header/prototypes.h"
#include "header/chemproto.h"

/*----------------------------------------------------------------------------*/
/*! \fn char *ath_strdup(const char *in)
 *  \brief This is really strdup(), but strdup is not available in 
 *   ANSI  (-pendantic or -ansi will leave it undefined in gcc)
 *   much like allocate.
 */

char *ath_strdup(const char *in)
{
  char *out = (char *)malloc((1+strlen(in))*sizeof(char));
  if(out == NULL) {
    ath_perr(-1,"ath_strdup: failed to alloc %d\n",(int)(1+strlen(in)));
    return NULL; /* malloc failed */
  }
  return strcpy(out,in);
}

/*----------------------------------------------------------------------------*/
/*! \fn int ath_gcd(int a, int b)
 *  \brief Calculate the Greatest Common Divisor by Euler's method
 */

int ath_gcd(int a, int b)
{
  int c;
  if(b>a) {c=a; a=b; b=c;} 
  while((c=a%b)) {a=b; b=c;}
  return b;
}

/*----------------------------------------------------------------------------*/
/*! \fn int ath_big_endian(void)
 *  \brief Return 1 if the machine is big endian (e.g. Sun, PowerPC)
 * return 0 if not (e.g. Intel)
 */

int ath_big_endian(void)
{
  short int n = 1;
  char *ep = (char *)&n;

  return (*ep == 0); /* Returns 1 on a big endian machine */
}

/*----------------------------------------------------------------------------*/
/*! \fn void ath_bswap(void *vdat, int len, int cnt)
 *  \brief Swap bytes, code stolen from NEMO  
 */
 
void ath_bswap(void *vdat, int len, int cnt)
{
  char tmp, *dat = (char *) vdat;
  int k;
 
  if (len==1)
    return;
  else if (len==2)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
      dat += 2;
    }
  else if (len==4)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
      tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
      dat += 4;
    }
  else if (len==8)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
      tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
      tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
      tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
      dat += 8;
    }
  else {  /* the general SLOOOOOOOOOW case */
    for(k=0; k<len/2; k++) {
      tmp = dat[k];
      dat[k] = dat[len-1-k];
      dat[len-1-k] = tmp;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void ath_error(char *fmt, ...)
 *  \brief Terminate execution and output error message
 *  Uses variable-length argument lists provided in <stdarg.h>
 */

void ath_error(char *fmt, ...)
{
  va_list ap;
   FILE *atherr = atherr_fp();

  fprintf(atherr,"### Fatal error: ");   /* prefix */
  va_start(ap, fmt);              /* ap starts with string 'fmt' */
  vfprintf(atherr, fmt, ap);      /* print out on atherr */
  fflush(atherr);                 /* flush it NOW */
  va_end(ap);                     /* end varargs */

#ifdef MPI_PARALLEL
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif

  exit(EXIT_FAILURE);
}

#if defined(PARTICLES) || defined(CHEMISTRY)
/*----------------------------------------------------------------------------*/
/*! \fn void ludcmp(Real **a, int n, int *indx, Real *d)
 *  \brief LU decomposition from Numerical Recipes
 *
 * Using Crout's method with partial pivoting
 * a is the input matrix, and is returned with LU decomposition readily made,
 * n is the matrix size, indx records the history of row permutation,
 * whereas d =1(-1) for even(odd) number of permutations.
 */
void ludcmp(Real **a, int n, int *indx, Real *d)
{
  int i,imax,j,k;
  Real big,dum,sum,temp;
  Real *rowscale;  /* the implicit scaling of each row */

  rowscale = (Real*)calloc_1d_array(n, sizeof(Real));
  *d=1.0;  /* No row interchanges yet */

  for (i=0;i<n;i++)
  { /* Loop over rows to get the implicit scaling information */
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) ath_error("[LUdecomp]:Input matrix is singular!");
    rowscale[i]=1.0/big;  /* Save the scaling */
  }

  for (j=0;j<n;j++) { /* Loop over columns of Crout's method */
    /* Calculate the upper block */
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    /* Calculate the lower block (first step) */
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      /* search for the largest pivot element */
      if ( (dum=rowscale[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    /* row interchange */
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      rowscale[imax]=rowscale[j];
    }
    indx[j]=imax; /* record row interchange history */
    /* Calculate the lower block (second step) */
    if (a[j][j] == 0.0) a[j][j]=TINY_NUMBER;
    dum=1.0/(a[j][j]);
    for (i=j+1;i<n;i++) a[i][j] *= dum;
  }
  free(rowscale);
}

/*----------------------------------------------------------------------------*/
/*! \fn void lubksb(Real **a, int n, int *indx, Real b[])
 *  \brief Backward substitution (from numerical recipies)
 *
 * a is the input matrix done with LU decomposition, n is the matrix size
 * indx id the history of row permutation
 * b is the vector on the right (AX=b), and is returned with the solution
 */
void lubksb(Real **a, int n, int *indx, Real *b)
{
  int i,ii=-1,ip,j;
  Real sum;
  /* Solve L*y=b */
  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  /* Solve U*x=y */
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void InverseMatrix(Real **a, int n, Real **b)
 *  \brief Inverse matrix solver
 *
 * a: input matrix; n: matrix size, b: return matrix
 * Note: the input matrix will be DESTROYED
 */
void InverseMatrix(Real **a, int n, Real **b)
{
  int i,j,*indx;
  Real *col,d;

  indx = (int*)calloc_1d_array(n, sizeof(int));
  col = (Real*)calloc_1d_array(n, sizeof(Real));

  ludcmp(a,n,indx,&d);

  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a, n, indx, col);
    for (i=0; i<n; i++)    b[i][j] = col[i];
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c)
 *  \brief Matrix multiplication: a(m*n) * b(n*l) = c(m*l) */
void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c)
{
  int i, j, k;
  for (i=0; i<m; i++)
    for (j=0; j<l; j++)
    {
      c[i][j] = 0.0;
      for (k=0; k<n; k++) c[i][j] += a[i][k] * b[k][j];
    }
}

/*--------------------------------------------------------------------------- */
/*! \fn Real Erf(Real z)
 *  \brief Error function  */
Real Erf(Real z)
{
  /* arrays of the error function y=erf(x) */
  static double x[101]={
        0.000000e+000,  3.783387e-003,  7.709914e-003,  1.178500e-002,  1.601425e-002,  2.040352e-002,
        2.495885e-002,  2.968653e-002,  3.459307e-002,  3.968525e-002,  4.497008e-002,  5.045486e-002,
        5.614715e-002,  6.205480e-002,  6.818596e-002,  7.454909e-002,  8.115295e-002,  8.800667e-002,
        9.511969e-002,  1.025018e-001,  1.101632e-001,  1.181145e-001,  1.263667e-001,  1.349310e-001,
        1.438193e-001,  1.530440e-001,  1.626176e-001,  1.725534e-001,  1.828652e-001,  1.935671e-001,
        2.046738e-001,  2.162008e-001,  2.281639e-001,  2.405796e-001,  2.534651e-001,  2.668380e-001,
        2.807169e-001,  2.951209e-001,  3.100699e-001,  3.255844e-001,  3.416859e-001,  3.583966e-001,
        3.757395e-001,  3.937386e-001,  4.124186e-001,  4.318054e-001,  4.519256e-001,  4.728071e-001,
        4.944786e-001,  5.169701e-001,  5.403124e-001,  5.645379e-001,  5.896800e-001,  6.157732e-001,
        6.428537e-001,  6.709587e-001,  7.001271e-001,  7.303990e-001,  7.618162e-001,  7.944220e-001,
        8.282614e-001,  8.633812e-001,  8.998296e-001,  9.376570e-001,  9.769156e-001,  1.017659e+000,
        1.059945e+000,  1.103830e+000,  1.149376e+000,  1.196644e+000,  1.245701e+000,  1.296614e+000,
        1.349454e+000,  1.404292e+000,  1.461205e+000,  1.520272e+000,  1.581573e+000,  1.645193e+000,
        1.711221e+000,  1.779746e+000,  1.850864e+000,  1.924673e+000,  2.001274e+000,  2.080774e+000,
        2.163281e+000,  2.248909e+000,  2.337778e+000,  2.430008e+000,  2.525728e+000,  2.625070e+000,
        2.728170e+000,  2.835170e+000,  2.946219e+000,  3.061469e+000,  3.181080e+000,  3.305216e+000,
        3.434048e+000,  3.567755e+000,  3.706521e+000,  3.850536e+000,  4.000000e+000
  };
  static double y[101]={
        0.00000000e+000,        4.26907434e-003,        8.69953340e-003,        1.32973284e-002,        1.80686067e-002,        2.30197153e-002,
        2.81572033e-002,        3.34878242e-002,        3.90185379e-002,        4.47565113e-002,        5.07091186e-002,        5.68839404e-002,
        6.32887618e-002,        6.99315688e-002,        7.68205444e-002,        8.39640613e-002,        9.13706742e-002,        9.90491090e-002,
        1.07008250e-001,        1.15257124e-001,        1.23804883e-001,        1.32660778e-001,        1.41834139e-001,        1.51334337e-001,
        1.61170754e-001,        1.71352743e-001,        1.81889576e-001,        1.92790394e-001,        2.04064148e-001,        2.15719527e-001,
        2.27764884e-001,        2.40208149e-001,        2.53056730e-001,        2.66317410e-001,        2.79996226e-001,        2.94098338e-001,
        3.08627885e-001,        3.23587825e-001,        3.38979770e-001,        3.54803790e-001,        3.71058224e-001,        3.87739454e-001,
        4.04841688e-001,        4.22356710e-001,        4.40273635e-001,        4.58578645e-001,        4.77254725e-001,        4.96281391e-001,
        5.15634428e-001,        5.35285634e-001,        5.55202571e-001,        5.75348359e-001,        5.95681482e-001,        6.16155658e-001,
        6.36719759e-001,        6.57317799e-001,        6.77889021e-001,        6.98368078e-001,        7.18685336e-001,        7.38767318e-001,
        7.58537287e-001,        7.77916009e-001,        7.96822665e-001,        8.15175962e-001,        8.32895397e-001,        8.49902691e-001,
        8.66123358e-001,        8.81488386e-001,        8.95935967e-001,        9.09413237e-001,        9.21877939e-001,        9.33299942e-001,
        9.43662512e-001,        9.52963249e-001,        9.61214608e-001,        9.68443923e-001,        9.74692870e-001,        9.80016358e-001,
        9.84480847e-001,        9.88162149e-001,        9.91142807e-001,        9.93509200e-001,        9.95348535e-001,        9.96745927e-001,
        9.97781755e-001,        9.98529475e-001,        9.99054014e-001,        9.99410828e-001,        9.99645625e-001,        9.99794704e-001,
        9.99885782e-001,        9.99939162e-001,        9.99969080e-001,        9.99985060e-001,        9.99993164e-001,        9.99997050e-001,
        9.99998805e-001,        9.99999548e-001,        9.99999841e-001,        9.99999948e-001,        9.99999985e-001
  };
  double z1, g;
  int il, iu, im;
  z1 = fabs(z);
  /* look for the location of z1 in the x array */
  il=0; iu=100;
 while(iu-il>1) {
    im = (iu+il)/2;
    if (x[im] > z1) iu = im;
    else il = im;
  }
  /* linear interpolation */
  g = (y[iu]*(z1-x[il])+y[il]*(x[iu]-z1))/(x[iu]-x[il]);

  g = MIN(g,1.0);
  if (z<0.0) g = -g;
  return g;
}

#endif /* PARTICLES or CHEMISTRY */
