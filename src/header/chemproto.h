#ifndef CHEMISTRY_PROTOTYPES_H
#define CHEMISTRY_PROTOTYPES_H 
#include "copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the src/chemistry
 *   directory.
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "defs.h"
#include "chemistry.h"

/*----------------------------------------------------------------------------*/
#ifdef CHEMISTRY

/*----------------------------------------------------------------------------*/
/* coeff.c */
void IonizationCoeff(ChemEvln *Evln, Real zeta_eff, Real Av, Real G, int verbose);
void CalCoeff(ChemEvln *Evln, Real T, int verbose);
void coeff_adj(ChemEvln *Evln);

Real ChemCoeff(Coefficient *coeff, Real T, int NumTRange);
Real IonGrCoeff(Chemistry *Chem, int i, Real T);
Real NeuGrCoeff(Chemistry *Chem, ChemEvln *Evln, int i);
Real DesorpCoeff(Coefficient *coeff, Real T);
Real GrGrCoeff(Chemistry *Chem, int i, Real T);
Real GrSurfCoeff(ChemEvln *Evln, int i, Real T);
void GrAvailFac(ChemEvln *Evln);
Real EleStickCoeff(Real size, int Z, Real T0);

/*----------------------------------------------------------------------------*/
/* density.c */
void init_numberden(ChemEvln *Evln, Real rho, int verbose);
void reset_numberden(ChemEvln *Evln, Real rho_new, int verbose);
void denscale(ChemEvln *Evln, int verbose);

/*----------------------------------------------------------------------------*/
/* diffusivity.c */
void Cal_NIMHD(ChemEvln *Evln);
void Cal_recomb(ChemEvln *Evln, Real Bmin, Real Bmax, int  nB,   Real Dt,
                                           Real *t_O, Real *t_H, Real *t_A);

/*----------------------------------------------------------------------------*/
/* disk.c */
void init_disk(Nebula *Disk);
Real Ionization_disk(Nebula *Disk, Real radius, Real z);
Real Ionization_disk1(Nebula *Disk, Real radius, Real sigmaz);

Real SurfDen_disk (Nebula *Disk, Real radius);
Real Temp_disk    (Nebula *Disk, Real radius);
Real Height_disk  (Nebula *Disk, Real radius);
Real Cs_disk      (Nebula *Disk, Real radius);
Real Omega_disk   (Nebula *Disk, Real radius);
Real eta0_disk    (Nebula *Disk, Real radius);
Real Rho_disk     (Nebula *Disk, Real radius, Real z);
Real SurfDenZ_disk(Nebula *Disk, Real radius, Real z);

/*----------------------------------------------------------------------------*/
/* evolve.c */
int evolve(Real te, Real dttry, Real err);
//int evolve(ChemEvln *Evln, Real te, Real *dttry, Real err);
void jacobi(ChemEvln *Evln, Real *numden, Real **jacob);
void derivs(ChemEvln *Evln, Real *numden, Real *drv);
//int EleMakeup(N_Vector &numden, int verbose);
int ConvertInd (int n);
int IConvertInd (int m);

/*----------------------------------------------------------------------------*/
/* init_chemistry.c */
void init_chemistry (Chemistry *Chem);
void init_chemevln  (Chemistry *Chem, ChemEvln *Evln);
void dup_chemevln(ChemEvln *Evln, ChemEvln *Evln_new);
void final_chemistry(Chemistry *Chem);
void final_chemevln (ChemEvln *Evln);

/*----------------------------------------------------------------------------*/
/* init_equations.c */
void init_equations(Chemistry *Chem);

/*----------------------------------------------------------------------------*/
/* init_reactions.c */
void init_reactions(Chemistry *Chem);

void InsertReactionInit(Chemistry *Chem);
int  CheckReaction(Chemistry *Chem, ReactionInfo Reaction);
int  FindInverse(Chemistry *Chem);
void PrintReaction(Chemistry *Chem, int i, Real K);

/*----------------------------------------------------------------------------*/
/* init_species.c */
void init_species(Chemistry *Chem);
int  FindSpecies(Chemistry *Chem, char name[NL_SPE]);

/*----------------------------------------------------------------------------*/
/* output_chemistry.c */
void init_chemout(Chemistry *Chem, ChemOutput *ChemOut, int mode,
                                   char *bname, Real bvalue,  char *id);
void final_chemout(ChemOutput *ChemOut);
void output_nspecies(ChemEvln *Evln, ChemOutput *ChemOut,
                           char *pname,   Real value);
void read_nspecies(ChemEvln *Evln, ChemOutput *ChemOut, int nskip);
void output_etaB(ChemEvln *Evln, ChemOutput *ChemOut, char *pname, Real value,
                     Real Omega, int nB);
void output_etafit(ChemEvln *Evln, ChemOutput *ChemOut,
                                     char *pname, Real value);
void output_recomb(ChemEvln *Evln, ChemOutput *ChemOut, char *pname, Real value,
                     Real Omega, Real Bmin, Real Bmax, int nB);
void ChemSet_allspecies(Chemistry *Chem, ChemOutput *ChemOut);
void ChemSet_allgrain(Chemistry *Chem, ChemOutput *ChemOut);
void ChemSet_allgas(Chemistry *Chem, ChemOutput *ChemOut);
void ChemSet_selected(Chemistry *Chem, ChemOutput *ChemOut);

/*----------------------------------------------------------------------------*/
/* stifbs.c */
int stifbs(ChemEvln *Evln, Real *y, Real *dydx, int nv, Real *xx,
               Real htry, Real eps, Real *yscal, Real *hdid, Real *hnext);

/*----------------------------------------------------------------------------*/
/* stifkr.c */
int stifkr(ChemEvln *Evln, Real *y, Real *dydx, int n, Real *x,
                Real htry, Real eps, Real *yscal, Real *hdid, Real *hnext);

#endif /* CHEMISTRY */

#endif /* CHEMISTRY_PROTOTYPES_H */
