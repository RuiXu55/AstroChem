#ifndef CHEMISTRY_H
#define CHEMISTRY_H
#include "../copyright.h"
/*==============================================================================
 * FILE: chemistry.h
 *
 * PURPOSE: Structure definition and global variable declaration for chemical
 *          reaction networks.
 *============================================================================*/

#include <stdio.h>
#include "../defs.h"

//#include "../config.h"

#ifdef CHEMISTRY

#define NL_ELE 5   /* maximum # of letters for an element */
#define NL_SPE 15  /* maximum # of letters for a species  */
#define NM_SIG 20  /* maximum # of single element species */

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

  int nreact;
  int nprod;

  Real Eb;            /* binding energy with grains (for neutrals ONLY) */

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
  char rtype[20];

  /* Reaction properties */
  int reactant[2];     /* label of the reactants; */
  int product[4];      /* label of the products; */

  /* Coefficients for calculating K */
  Coefficient *coeff;  /* have different meanings for different rtypes */

}ReactionInfo;

/*-----------------------------------------------------------------------------
 * Global information of the chemistry model (independent of cells)
 */
typedef struct Chemistry_s {

  int N_Ele;          /* # of elements in chemical reaction network */
  int N_Neu_f;	      /* # of neutral species with +/- ionized counterpart */
  int N_Neu;          /* # of neutral species with + ionized counterpart */
  int N_Neu_s;        /* # of neutral species without an inoized counterpart */
  int N_Ion_s;        /* # of charged species without neutral counterpart */
  int Ntot;           /* Total Number of species */

  int NeuInd;
  int SNeuInd;        /* Starting index of special neutrals (if exist)*/
  int SIonInd;        /* Starting index of special ions (if exist) */

  int NReaction;      /* Total Number of reactions */
  int Reaction_size;  /* Size of the reaction array */

  /* Array to store the information of all the elements */
  ElementInfo *Elements;     /* 0..N_Ele-1 */

  /* Array to store the information of all the species */
  SpeciesInfo *Species;      /* 0..Ntot-1 */

  /* Array of all reactions */
  ReactionInfo *Reactions;   /* 0..Reaction_size-1 */

}Chemistry;

#endif /* CHEMISTRY */

#endif /* CHEMISTRY_H */
