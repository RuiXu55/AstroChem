#ifndef CHEMISTRY_H
#define CHEMISTRY_H
/*==============================================================================
 * FILE: chemistry.h
 *
 * PURPOSE: Structure definition and global variable declaration for chemical
 *          reaction networks.
 *============================================================================*/

#include <stdio.h>
#include "defs.h"

//#include "../config.h"

#ifdef CHEMISTRY

#define NL_ELE 5   /* maximum # of letters for an element */
#define NL_SPE 20  /* maximum # of letters for a species  */
#define NM_SIG 20  /* maximum # of single element species */

#define OneYear 3.15576e7   /* one year in unit of second */
#define ChemErr 0.001  /* maximum allowable error in the chemical evolution */
#define MUN 2.34   /* mean molecular weight */

/*----------------------------------------------------------------------------*/
/***************************** Structure Definition ***************************/
/*----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Information about one element
 */
typedef struct ElementInfo_s {

  char name[NL_ELE];  /* name of the elements */
  Real mass;          /* in unit of m_p */
  Real abundance;     /* initial abundance */

  int numsig;         /* number of its single element species */
  int single[NM_SIG]; /* the label of all such species */

}ElementInfo;

/*-----------------------------------------------------------------------------
 * Information about one species
 */
typedef struct SpeciesInfo_s {

  char name[NL_SPE];  /* name of the species */
  Real mass;          /* in unit of m_p */
  int charge;         /* charge of the species */
  int numelem;        /* number of elements it contains */
  int *composition;   /* Chemical composition (0..NumElement) */
  short type;	      /* 0: electron; 1: neutral; 2: ion; 3: grain; 4: mantle */

  Real Eb;            /* binding energy with grains (for neutrals ONLY) */
  Real gsize;         /* GrainSize (for grains ONLY) */

}SpeciesInfo;

/*-----------------------------------------------------------------------------
 * Reaction coefficient
 */
typedef struct Coefficient_s {

/* K=alpha*pow(T/300.0, beta)*exp(-gamma/T) for chemical reactions */
/* Definition is different for reactions with grains */
  Real alpha;
  Real beta;
  Real gamma;

  Real Tmin;    /* Minimum temperature */
  Real Tmax;    /* Maximum temperature */

}Coefficient;

/*-----------------------------------------------------------------------------
 * Information about one reaction
 */
typedef struct ReactionInfo_s {

  /* Reaction type */
  int rtype;           /* 0: ionization; 1: chemical; >=2: grain related */

  /* Reaction properties */
  int reactant[2];     /* label of the reactants; */
  int product[4];      /* label of the products; */
  int NumTRange;       /* Number of temperature range */
  int inverse;         /* Label of its inverse reaction; =-1 if none */
  int use;             /* in use (1) or not (0) */

  /* Coefficients for calculating K */
  Coefficient *coeff;  /* have different meanings for different rtypes */

}ReactionInfo;

/*-----------------------------------------------------------------------------
 * Information about one reaction term in the reaction equation
 */
typedef struct EquationTerm_s {

  int ind;             /* Reaction index (from which one can get K) */
  int N;               /* Number of species involved (<=3) */
  int lab[3];          /* The label of the species */
  int dir;             /* Reactant (-1) or Product (1) */

}EquationTerm;

/*-----------------------------------------------------------------------------
 * All reactions concerned with one species
 */
typedef struct EquationInfo_s {

  int NTerm;           /* Number of reactions associated with this species */
  int EqTerm_size;     /* Size of the EqTerm array */

  EquationTerm *EqTerm;

}EquationInfo;


/*-----------------------------------------------------------------------------
 * Global information of the chemistry model (independent of cells)
 */
typedef struct Chemistry_s {

  int N_Ele;          /* # of elements in chemical reaction network */
  int N_Ele_tot;      /* # of all elements (including grains) */
  int N_Neu_f;        /* # of neutral species with +/- ionized counterpart */
  int N_Neu;          /* # of neutral species with + ionized counterpart */
  int N_Neu_s;        /* # of neutral species without inoized counterpart */
  int N_Ion_s;        /* # of charged species without neutral counterpart */
  int Ntot;           /* Total Number of species */

  int NGrain;         /* # of grain types */
  int GrCharge;       /* Maximum charge of a grain */
  Real GrDen;         /* Grain mass density, g/cm^3 */
  Real *GrSize;       /* Grain size, micron */
  Real *GrFrac;       /* Grain mass fraction */

  int NeuInd;         /* Starting index of normal neutrals */
  int SNeuInd;        /* Starting index of special neutrals (if exist)*/
  int SIonInd;        /* Starting index of special ions (if exist) */
  int GrInd;          /* Starting index of grains (if exist) */

  int NReaction;      /* Total Number of reactions */
  int Reaction_size;  /* Size of the reaction array */

  /* Array to store the information of all the elements */
  ElementInfo *Elements;     /* 0..N_Ele+NGrain-1 */

  /* Array to store the information of all the species */
  SpeciesInfo *Species;      /* 0..Ntot-1 */

  /* Array of all reactions */
  ReactionInfo *Reactions;   /* 0..Reaction_size-1 */

  /* Array of evolution equations of all species */
  EquationInfo *Equations;   /* 0..Ntot-1 */

}Chemistry;

Chemistry Chem;

/*-----------------------------------------------------------------------------
 * Structure for a complete chemistry evolution model
 */
typedef struct ChemEvln_s {

  /* Parent chemistry model */
  Chemistry *Chem;

  /* MHD quantities */
  Real rho, T, B;

  /* ionization rate */
  Real zeta_eff;

  /* Current time */
  Real t;

  /* Reaction rate coefficient array */
  Real *K;                   /* 0..NReaction-1 */

  Real *rate_adj;  /* dynamical adjustment to the rate coefficient */

  /* Number density of all species */
  Real *NumDen;              /* 0..Ntot-1 */

  /* Fraction of grain mantle species available for surface reactions */
  Real *GrAvail;       /* 0..Ntot-1 */

  /* Abundance to number density ratio */
  Real Abn_Den;              /* abundance / number density (1/n_H) */

  /* Density scale of each species */
  Real *DenScale;            /* 0..Ntot-1 */

  /* Magnetic diffusivity */
  Real eta_O, eta_H, eta_A;

}ChemEvln;

ChemEvln Evln;

/*-----------------------------------------------------------------------------
 * Output parameters
 */
typedef struct ChemOutput_s{

  /* name of the output file */
  char outbase[30];
  char outid[30];

  /* background information */
  char bname[30];
  Real bvalue;

  /* output number label */
  int lab;

  /* Four output modes: */
  int mode;  /* 1: number density; 2: diffusivity;
                3: recombination time; 4: diffusivity fitting parameters */

  /* For Mode 1: A list of selected species */
  int nsp;
  int *ind;

}ChemOutput;


/*-----------------------------------------------------------------------------
 * Disk Model parameters
 */
typedef struct Nebula_s{

  Real Mstar;		/* Stellar mass (M_sun) */

  Real Sigma0;		/* Disk surface density at 1AU (g/cm^2) */
  Real pS;		/* density profile power law index */

  Real Tdisk;		/* Disk temperature at 1AU (K) */
  Real pT;		/* temperature profile power law index */

  Real CR_rate;		/* Cosmic ray ionization rate (1e-17 /s) */
  Real RD_rate;		/* ionization rate by radioactive decay (1e-19 /s) */

  Real Lx;		/* X-ray flux (10^29 erg/s) */
  Real Tx;		/* X-ray source temperature (keV) */

}Nebula;


#endif /* CHEMISTRY */

#endif /* CHEMISTRY_H */
