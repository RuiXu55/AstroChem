#include "../copyright.h"
/*=============================================================================
 * FILE: init_species.c
 *
 * PURPOSE: Initialize chemical species. It reads the elemental and species
 *   information from standard input file "SpeciesNames.txt" in the same
 *   directory as the executable.  Below is a sample of the input file format:
 *
 *  # Number of Elements
 *      2
 *  # Number of Grain Types
 *      1
 *  # Maximum Grain Charge
 *      2
 *  # Element       Mass (m_p)      Abundance (H=1)
 *      H           1               1.0
 *      Mg          24              1.25e-8
 *  # Grain mass density (g/cm^3)
 *      3
 *  # Gr-Size       MassRatio
 *      0.1         0.01
 *  # Number of Neutral Species (with ion counterpart)
 *      2
 *  # Number of Neutral Species (without ion counterpart)
 *      0
 *  # Number of Ion Species (without neutral counterpart)
 *      0
 *  # Species       E_B (K)
 *      H2          450.0
 *      Mg          5300.0
 *
 *  This code further construct all grain related species, such as neutral and
 *  charged grains, mantle species, etc. For details on species construction,
 *  see Bai & Goodman (2009).
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_species()
 *   FindSpecies()
 *
 * REFERENCES:
 *   Bai, X.-N. & Goodman, J., 2009, ApJ, 701, 737
 *
 * History:
 *   Written by  Xuening Bai      Nov. 2010
==============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "chemistry.h"
#include "../prototypes.h"
#include "prototypes.h"

#ifdef CHEMISTRY  /* endif at the end of the file */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   Analyze()    - analyze the composition of a chemical species
 *   FindElem()   - find an element from the species
 *============================================================================*/
void Analyze(Chemistry *Chem, int i);
int  FindElem(Chemistry *Chem, char name[NL_SPE], int p, int *l);
void OutputSpecies(Chemistry *Chem);

/*============================================================================*/
/*---------------------------- Public Functions ------------------------------*/

/*----------------------------------------------------------------------------*/
/* The First step of initializing the chemical model
 * Read and analyze chemical species names, add grain and mantle species
 */
void init_species(Chemistry *Chem)
{
  FILE *fp;
  int i, j, k, p, n;
  char line[MAXLEN],mantle[8];
  Real sumgas, grtot, ratio;

  fp = fopen("SpeciesNames.txt","r");

  if(fp == NULL){
    ath_perr(-1,"[construct_species]: Unable to open the species name file!\n");
    return;
  }

/*---------------------- Read and construct Elements -------------------------*/

/* Elements info */

  fgets(line,MAXLEN,fp);
  fscanf(fp, "%d\n", &(Chem->N_Ele));
  if (Chem->N_Ele <= 0)
    ath_error("[init_species]: Number of element must be positive!\n");

  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);

/* Initiate element arrays */

  Chem->Elements = (ElementInfo*)calloc_1d_array(Chem->N_Ele,
                                                 sizeof(ElementInfo));

/* Read elements info */

  sumgas = 0.0;

  fgets(line,MAXLEN,fp);
  for (i=0; i<Chem->N_Ele; i++)
  {
    fscanf(fp, "%s", Chem->Elements[i].name);
    fscanf(fp, "%lf", &(Chem->Elements[i].mass));
    fscanf(fp, "%lf\n", &(Chem->Elements[i].abundance));

    sumgas += Chem->Elements[i].mass * Chem->Elements[i].abundance;

    Chem->Elements[i].numsig = 0;
    for (j=0; j<NM_SIG; j++)
      Chem->Elements[i].single[j] = 0;

    if ((Chem->Elements[i].mass<=0.0) || (Chem->Elements[i].abundance<=0.0))
      ath_error("[init_species]: mass and abundance of %s must be positive!\n",
                 Chem->Elements[i].name);
  }

  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);
  fgets(line,MAXLEN,fp);
/*------------------------ Read and construct Species ------------------------*/

/* Number of species */

  fgets(line,MAXLEN,fp); /* Neutral species with +/- ion counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Neu_f));
  if (Chem->N_Neu_f < 0)
    ath_error("[init_species]: Number of neutral species must be positive!\n");

  fgets(line,MAXLEN,fp); /* Neutral species with + ion counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Neu));
  if (Chem->N_Neu < 0)
    ath_error("[init_species]: Number of neutral species must be positive!\n");

  fgets(line,MAXLEN,fp); /* Neutral species without ion counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Neu_s));
  if (Chem->N_Neu_s < 0)
    ath_error("[init_species]: Number of special neutral species must >=0 !\n");

  fgets(line,MAXLEN,fp); /* Ion species without neutral counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Ion_s));
  if (Chem->N_Ion_s < 0)
    ath_error("[init_species]: Number of special ion species must be >=0!\n");

  Chem->Ntot = 3*Chem->N_Neu_f + 2*Chem->N_Neu + Chem->N_Neu_s + Chem->N_Ion_s + 1;

/* Initialize arrays */

  Chem->Species  = (SpeciesInfo*)calloc_1d_array(Chem->Ntot,
                   sizeof(SpeciesInfo));

  for (i=0; i<Chem->Ntot; i++) {
    Chem->Species[i].composition = (int*)calloc(Chem->N_Ele, sizeof(int));
    Chem->Species[i].nreact = 0;
    Chem->Species[i].nprod  = 0;
  }

/* The first species: electron */
  sprintf(Chem->Species[0].name, "e-");

  Chem->Species[0].mass = 1.0/1836.0;
  Chem->Species[0].charge = -1;

  Chem->Species[0].numelem = 0;    /* Electron does not contain any element */
  for (j=0; j<Chem->N_Ele; j++)
    Chem->Species[0].composition[j] = 0;

  Chem->Species[0].Eb = 0.0;      /* Not used for electron */

/* Read and construct Neutral species with +/- ionized counterpart */

  fgets(line,MAXLEN,fp);
  for (i=1; i<=Chem->N_Neu_f; i++) 
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    fscanf(fp, "%lf\n", &(Chem->Species[i].Eb));    /* Binding energy */

    if (Chem->Species[i].Eb < 0)
      ath_error("[init_species]: Binding energy must be non-negative!\n");

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);
    
    /* Construct positively ionized counterpart */
    n = Chem->N_Neu_f + i;
  
    strcpy(Chem->Species[n].name, Chem->Species[i].name);
    strcat(Chem->Species[n].name, "+");

    Chem->Species[n].mass = Chem->Species[i].mass;
    Chem->Species[n].charge = 1;
  
    Chem->Species[n].numelem = Chem->Species[i].numelem;
    for (j=0; j<Chem->N_Ele; j++)
      Chem->Species[n].composition[j] = Chem->Species[i].composition[j];

    Chem->Species[n].Eb = Chem->Species[i].Eb;

    /* Construct negatively ionized counterpart */
    n = 2*Chem->N_Neu_f + i;

    strcpy(Chem->Species[n].name, Chem->Species[i].name);
    strcat(Chem->Species[n].name, "-");

    Chem->Species[n].mass = Chem->Species[i].mass;
    Chem->Species[n].charge = -1; 

    Chem->Species[n].numelem = Chem->Species[i].numelem;
    for (j=0; j<Chem->N_Ele; j++)
      Chem->Species[n].composition[j] = Chem->Species[i].composition[j];
    
    Chem->Species[n].Eb = Chem->Species[i].Eb;
  }

/* Read and construct species for Neutral species with an ionized counterpart */

  p = 3*Chem->N_Neu_f+1;
  Chem->NeuInd = p; /* Starting index of special neutrals */

  for (i=p; i<p+Chem->N_Neu; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    fscanf(fp, "%lf\n", &(Chem->Species[i].Eb));    /* Binding energy */

    if (Chem->Species[i].Eb < 0)
      ath_error("[init_species]: Binding energy must be non-negative!\n");

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);

    /* Construct ionized counterpart */
    n = Chem->N_Neu + i;

    strcpy(Chem->Species[n].name, Chem->Species[i].name);
    strcat(Chem->Species[n].name, "+");

    Chem->Species[n].mass = Chem->Species[i].mass;
    Chem->Species[n].charge = 1;

    Chem->Species[n].numelem = Chem->Species[i].numelem;
    for (j=0; j<Chem->N_Ele; j++)
      Chem->Species[n].composition[j] = Chem->Species[i].composition[j];

    Chem->Species[n].Eb = Chem->Species[i].Eb;
  }

/* Read and construct species for Neutral species without ionized counterpart */

  p += 2*Chem->N_Neu;
  Chem->SNeuInd = p; /* Starting index of special neutrals */

  for (i=p; i<p+Chem->N_Neu_s; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    fscanf(fp, "%lf", &(Chem->Species[i].Eb));

    if (Chem->Species[i].Eb < 0)
      ath_error("[init_species]: Binding energy must be non-negative!\n");

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);
  }

/* Read and construct ionized species without a neutral counterpart */

  p += Chem->N_Neu_s;
  Chem->SIonInd = p; /* Starting index of special ions */

  for (i=p; i<p+Chem->N_Ion_s; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    Chem->Species[i].Eb = 0.0;        /* Not used */

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);
  }

  fclose(fp);

/*----------- Output All Species -------------*/

  OutputSpecies(Chem);

  return;
}

/*---------------------------------------------------------------------------*/
/* Find the label of a chemical species in the catalog
 */
int FindSpecies(Chemistry *Chem, char name[NL_SPE])
{
  int i, k;

  for (i=0; i<Chem->Ntot; i++)
  {
    k = strcmp(name,Chem->Species[i].name);
    if (k == 0)
      return i;
  }

  if (strcmp(name,"0") == 0)
    return -1;
  else
    return -10;   /* not found! */
}


/*============================================================================*/
/*------------------------------ PRIVATE FUNCTIONS ---------------------------*/

/*----------------------------------------------------------------------------*/
/* Analyze the composition of a species i and do the related operations
 */
void Analyze(Chemistry *Chem, int i)
{/* Warning: This function ONLY analyze chemical species, no grain! */

  int j, k, n, l=0;
  int p = 0;
  Real mass = 0.0;

/* Find out the compositon of the species */
  for (j=0; j<Chem->N_Ele; j++)
    Chem->Species[i].composition[j] = 0;

  while ((p<(NL_SPE-1)) && (Chem->Species[i].name[p] != '\0')
                        && (Chem->Species[i].name[p] != '+'))
  {
    k = FindElem(Chem, Chem->Species[i].name, p, &l);
    p += l;
    n = 1;

    if ((Chem->Species[i].name[p]>'1') && (Chem->Species[i].name[p]<='9'))
    {
       n = Chem->Species[i].name[p]-'0';
       p += 1;
    }

    Chem->Species[i].composition[k] += n;
  }

  if (Chem->Species[i].name[p] == '+')
    Chem->Species[i].charge = 1;
  else
    Chem->Species[i].charge = 0;/* Note: we don't have negative charged ions! */

  if (p >= (NL_SPE-1))
    ath_perr(-1, "[analyze_species]: Too many letters in  %s!\n",
                 Chem->Species[i].name);

/* Calculate the mass and number of elements */
  k = 0;
  for (j=0; j<Chem->N_Ele; j++)
  if (Chem->Species[i].composition[j] != 0)
  {
    mass += Chem->Species[i].composition[j]*Chem->Elements[j].mass;
    k += 1;
    l = j;
  }

  Chem->Species[i].mass = mass;
  Chem->Species[i].numelem = k;

/* If it is a single neutral element species, tell the array Chem->Elements[] */
  if ((k == 1) && (Chem->Species[i].charge == 0))
  {
    n = Chem->Elements[l].numsig;
    Chem->Elements[l].numsig += 1;
    Chem->Elements[l].single[n] = i;      /* Neutral single element */
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* Find the element names from the p_th letter of a species
 * l: the length of the element name (1 or 2)
 * Warning: Only find chemical elements, NOT grain!
 */
int FindElem(Chemistry *Chem, char name[NL_SPE], int p, int *l)
{
  int i = 0;

  /* If the second letter is not capital, then this element has 2 letters */
  if ((name[p+1]>='a') && (name[p+1]<='z'))
  {
    *l = 2;
    while ((i<Chem->N_Ele) && !((Chem->Elements[i].name[0] == name[p])
                           && (Chem->Elements[i].name[1] == name[p+1]))) 
      i++;
  }
  else    /* this element has 1 letter */
  {
    *l = 1;
    while ((i<Chem->N_Ele) && ((Chem->Elements[i].name[0] != name[p])
                           || (Chem->Elements[i].name[1] != '\0')))
      i++;
  }

  if (i == Chem->N_Ele)
    ath_perr(-1, "[find_element]: Element not found in %s!\n", name);

  return i;
}

/*---------------------------------------------------------------------------*/
/* Output species info
 */
void OutputSpecies(Chemistry *Chem)
{
  int i,j,k;
  SpeciesInfo *sp;

  ath_pout(0,"\n");
  ath_pout(0,"Total Number of Elements:          %d\n",Chem->N_Ele);

  ath_pout(0,"Number of Full Neutral Species:  %d\n",Chem->N_Neu_f);
  ath_pout(0,"Number of Normal Neutral Species:  %d\n",Chem->N_Neu);
  ath_pout(0,"Number of Special Neutral Species: %d\n",Chem->N_Neu_s);
  ath_pout(0,"Number of Special Ion Species:     %d\n",Chem->N_Ion_s);
  ath_pout(0,"Total Number of Species:           %d\n\n",Chem->Ntot);

  if (Chem->N_Neu_s > 0)
    ath_pout(0,"Starting index for Special Neutrals: %d:\n", Chem->SNeuInd);

  if (Chem->N_Ion_s > 0)
    ath_pout(0,"Starting index for Special Ions:     %d:\n", Chem->SIonInd);

  ath_pout(0,"\n");
  ath_pout(0,"List of species:\n");
  ath_pout(0,"\n");
  for (i=0; i<Chem->Ntot; i++)
  {
    sp = &(Chem->Species[i]);
    ath_pout(0,"%4d. %7s: M=%9.2e m_p, Q=%2de, N_Ele=%2d, Eb=%6.1f", i, sp->name,
                                sp->mass, sp->charge, sp->numelem, sp->Eb);
    ath_pout(0,"  Comp:");
    for (j=0; j<Chem->N_Ele; j++)
    {
      k = sp->composition[j];
      if (k > 0)
        ath_pout(0,"  %3s:%2d",Chem->Elements[j].name, k);
    }
    ath_pout(0,"\n");
  }
  ath_pout(0,"\n");

  return;
}

#endif /* CHEMISTRY */
