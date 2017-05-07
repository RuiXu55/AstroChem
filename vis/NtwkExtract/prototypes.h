#ifndef PROTOTYPES_H
#define PROTOTYPES_H 
#include "copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions from the /src directory.      */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "defs.h"

/*----------------------------------------------------------------------------*/
/* main.c */
int athena_main(int argc, char *argv[]);

/*----------------------------------------------------------------------------*/
/* ath_array.c */
void*   calloc_1d_array(                      size_t nc, size_t size);
void**  calloc_2d_array(           size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);
void free_1d_array(void *array);
void free_2d_array(void *array);
void free_3d_array(void *array);

/*----------------------------------------------------------------------------*/
/* ath_log.c */
void ath_log_set_level(const int out, const int err);
void ath_log_open(const char *basename, const int lazy, const char *mode);
void ath_log_close(void);
FILE *athout_fp(void);
FILE *atherr_fp(void);
void ath_flush_out(void);
void ath_flush_err(void);
int ath_perr(const int level, const char *fmt, ...);
int ath_pout(const int level, const char *fmt, ...);

/*----------------------------------------------------------------------------*/
/* par.c */
void   par_open(char *filename);
void   par_cmdline(int argc, char *argv[]);
int    par_exist(char *block, char *name);

char  *par_gets(char *block, char *name);
int    par_geti(char *block, char *name);
double par_getd(char *block, char *name);

char  *par_gets_def(char *block, char *name, char   *def);
int    par_geti_def(char *block, char *name, int    def);
double par_getd_def(char *block, char *name, double def);

void   par_sets(char *block, char *name, char *sval, char *comment);
void   par_seti(char *block, char *name, char *fmt, int ival, char *comment);
void   par_setd(char *block, char *name, char *fmt, double dval, char *comment);

void   par_dump(int mode, FILE *fp);
void   par_close(void);

/*----------------------------------------------------------------------------*/
/* utils.c */
char *ath_strdup(const char *in);
int ath_gcd(int a, int b);
int ath_big_endian(void);
void ath_bswap(void *vdat, int sizeof_len, int cnt);
void ath_error(char *fmt, ...);
void ludcmp(Real **a, int n, int *indx, Real *d);
void lubksb(Real **a, int n, int *indx, Real *b);
void InverseMatrix(Real **a, int n, Real **b);
void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c);
Real Erf(Real z);
#endif /* PROTOTYPES_H */
