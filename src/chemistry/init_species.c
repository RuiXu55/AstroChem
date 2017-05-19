#include "../header/copyright.h"
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
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/chemproto.h"
#include "../header/prototypes.h"

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
  int i, j, k, p, n,TmpInd;
  char line[MAXLEN],mantle[8],fname[20];
  Real sumgas, grtot, ratio;
  ath_pout(0,"Read species begin!\n");
  sprintf(fname,"%s",par_gets("job","read_species"));
  fp = fopen(fname,"r");

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
  fscanf(fp, "%d\n", &(Chem->NGrain));
  if (Chem->NGrain < 0)
    ath_error("[init_species]: Number of grain types must be non-negative!\n");

  fgets(line,MAXLEN,fp);
  fscanf(fp, "%d\n", &(Chem->GrCharge));
  if (Chem->GrCharge <= 0)
    ath_error("[init_species]: Number of grain charges must be positive!\n");

/* Initiate element arrays */

  Chem->N_Ele_tot= Chem->N_Ele+Chem->NGrain;

  Chem->Elements = (ElementInfo*)calloc_1d_array(Chem->N_Ele_tot,
                                                 sizeof(ElementInfo));

  if (Chem->NGrain > 0) {
    Chem->GrSize   = (Real*)calloc_1d_array(Chem->NGrain, sizeof(Real));

    Chem->GrFrac   = (Real*)calloc_1d_array(Chem->NGrain, sizeof(Real));
  }

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

/* Read grain info */

  grtot = 0.0;

  fgets(line,MAXLEN,fp);
  fscanf(fp, "%lf\n", &(Chem->GrDen));
  if (Chem->GrCharge < 0)
    ath_error("[init_species]: Grain mass density must be non-negative!\n");

  fgets(line,MAXLEN,fp);
  for (i=0; i<Chem->NGrain; i++)
  {
    k = Chem->N_Ele+i;  /* label of grain element */

    sprintf(Chem->Elements[k].name, "gr%d", i+1);/* Grain as an element */

    fscanf(fp, "%lf", &(Chem->GrSize[i]));     /* Size of this grain */
    fscanf(fp, "%lf\n", &(Chem->GrFrac[i]));   /* Mass fraction of this grain */

    /* Grain mass (m_p) */
    Chem->Elements[k].mass = 2.505e12*Chem->GrDen*pow(Chem->GrSize[i],3.0);

    grtot += Chem->GrFrac[i];

    Chem->Elements[k].numsig = 0;
    for (j=0; j<NM_SIG; j++)
      Chem->Elements[k].single[j] = 0;

    if ((Chem->GrSize[i]<=0.0) || (Chem->GrFrac[i]<=0.0) || (grtot>=1.0))
      ath_error("[init_species]: Grain size and mass fraction must be >0!\n");
  }

/* Calculate grain abundance */

  ratio = sumgas / (1.0-grtot);

  for (i=0; i<Chem->NGrain; i++)
  {
    k = Chem->N_Ele+i;  /* label of grain element */

    Chem->Elements[k].abundance
         = Chem->GrFrac[i]/Chem->Elements[k].mass*ratio;
    ath_pout(0,"grain abundance=%1.12f",Chem->Elements[k].abundance);
  }

/*------------------------ Read and construct Species ------------------------*/

/* Number of species */

  fgets(line,MAXLEN,fp); /* Neutral species with ion counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Neu_f));
  if (Chem->N_Neu_f < 0)
    ath_error("[init_species]: Number of full neutral species must be >=0!\n");

  fgets(line,MAXLEN,fp); /* Neutral species with ion counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Neu));
  if (Chem->N_Neu < 0)
    ath_error("[init_species]: Number of neutral species must be >=0!\n");

  fgets(line,MAXLEN,fp); /* Neutral species without ion counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Neu_s));
  if (Chem->N_Neu_s < 0)
    ath_error("[init_species]: Number of special neutral species must >=0 !\n");

  fgets(line,MAXLEN,fp); /* Ion species without neutral counterpart */
  fscanf(fp, "%d\n", &(Chem->N_Ion_s));
  if (Chem->N_Ion_s < 0)
    ath_error("[init_species]: Number of special ion species must be >=0!\n");

  Chem->Ntot = (3+Chem->NGrain)*Chem->N_Neu_f + (2+Chem->NGrain)*Chem->N_Neu
             + (1+Chem->NGrain)*Chem->N_Neu_s
             + Chem->N_Ion_s + Chem->NGrain*(2*Chem->GrCharge+1) + 1;

/* Initialize arrays */

  Chem->Species  = (SpeciesInfo*)calloc_1d_array(Chem->Ntot,
                   sizeof(SpeciesInfo));

  for (i=0; i<Chem->Ntot; i++) {
    Chem->Species[i].composition = (int*)calloc(Chem->N_Ele_tot, sizeof(int));
  }

/* The first species: electron */
  sprintf(Chem->Species[0].name, "e-");

  Chem->Species[0].mass = 1.0/1836.0;
  Chem->Species[0].charge = -1;
  Chem->Species[0].type = 0; /* electron */

  Chem->Species[0].numelem = 0;    /* Electron does not contain any element */
  for (j=0; j<Chem->N_Ele_tot; j++)
    Chem->Species[0].composition[j] = 0;

  Chem->Species[0].Eb = 0.0;      /* Not used for electron */
  Chem->Species[0].gsize = 0.0;   /* Not used for electron */

/* Read and construct Neutral species with +/- ionized counterpart */

  fgets(line,MAXLEN,fp);
  for (i=1; i<=Chem->N_Neu_f; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    fscanf(fp, "%lf\n", &(Chem->Species[i].Eb));    /* Binding energy */
    Chem->Species[i].gsize = 0.0;                   /* Not used */
    Chem->Species[i].type  = 1;  /* neutral */

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
    Chem->Species[n].type   = 2;  /* ion */

    Chem->Species[n].numelem = Chem->Species[i].numelem;
    for (j=0; j<Chem->N_Ele_tot; j++)
      Chem->Species[n].composition[j] = Chem->Species[i].composition[j];

    Chem->Species[n].Eb = Chem->Species[i].Eb;
    Chem->Species[n].gsize = 0.0;     /* Not used */

    /* Construct negatively ionized counterpart */
    n = 2*Chem->N_Neu_f + i;

    strcpy(Chem->Species[n].name, Chem->Species[i].name);
    strcat(Chem->Species[n].name, "-");

    Chem->Species[n].mass = Chem->Species[i].mass;
    Chem->Species[n].charge = -1;
    Chem->Species[n].type   = 2; /* ion */

    Chem->Species[n].numelem = Chem->Species[i].numelem;
    for (j=0; j<Chem->N_Ele_tot; j++)
      Chem->Species[n].composition[j] = Chem->Species[i].composition[j];

    Chem->Species[n].Eb = Chem->Species[i].Eb;
    Chem->Species[n].gsize = 0.0;     /* Not used */

		ath_pout(0,"Neu_f= %s\n",Chem->Species[i].name);
  }

/* Read and construct Neutral species with + ionized counterpart */

  p = 3*Chem->N_Neu_f+1;
  Chem->NeuInd = p; /* Starting index of normal neutrals */

  for (i=p; i<p+Chem->N_Neu; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    fscanf(fp, "%lf\n", &(Chem->Species[i].Eb));    /* Binding energy */
		ath_pout(0,"N_Neu= %s\n",Chem->Species[i].name);
    Chem->Species[i].gsize = 0.0;                   /* Not used */
    Chem->Species[i].type   = 1;

    if (Chem->Species[i].Eb < 0)
      ath_error("[init_species]: Binding energy must be non-negative!\n");

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);

    /* Construct ionized counterpart */
    n = i + Chem->N_Neu;

    strcpy(Chem->Species[n].name, Chem->Species[i].name);
    strcat(Chem->Species[n].name, "+");

    Chem->Species[n].mass = Chem->Species[i].mass;
    Chem->Species[n].charge = 1;
    Chem->Species[n].type   = 2;

    Chem->Species[n].numelem = Chem->Species[i].numelem;
    for (j=0; j<Chem->N_Ele_tot; j++)
      Chem->Species[n].composition[j] = Chem->Species[i].composition[j];

    Chem->Species[n].Eb = Chem->Species[i].Eb;
    Chem->Species[n].gsize = 0.0;     /* Not used */
  }

/* Read and construct species for Neutral species without ionized counterpart */

  p = Chem->NeuInd + 2*Chem->N_Neu;
  Chem->SNeuInd = p; /* Starting index of special neutrals */

  for (i=p; i<p+Chem->N_Neu_s; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    fscanf(fp, "%lf", &(Chem->Species[i].Eb));
		ath_pout(0,"Neu_s= %s\n",Chem->Species[i].name);
    Chem->Species[i].gsize = 0.0;                   /* Not used */
    Chem->Species[n].type   = 1;

    if (Chem->Species[i].Eb < 0)
      ath_error("[init_species]: Binding energy must be non-negative!\n");

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);

  }

/* Read and construct ionized species without a neutral counterpart */

  p = Chem->SNeuInd + Chem->N_Neu_s;
  Chem->SIonInd = p; /* Starting index of special ions */

  for (i=p; i<p+Chem->N_Ion_s; i++)
  {
    fscanf(fp, "%s", Chem->Species[i].name);
    Chem->Species[i].Eb = 0.0;        /* Not used */
    Chem->Species[n].gsize = 0.0;     /* Not used */

    /* Analyze the name of this species to obtain its compositon (IMPORTANT!) */
    Analyze(Chem, i);
    Chem->Species[n].type   = 2;
  }

/* Construct grains */
  p = p + Chem->N_Ion_s - 1;
  Chem->GrInd = p + 1;  /* Starting index of grains (if exist) */

  for (k=0; k<Chem->NGrain; k++)
  {
    n = p+k*(2*Chem->GrCharge+1);

    for (i=1; i<Chem->GrCharge; i++)
      sprintf(Chem->Species[n+i].name, "gr%d(%d-)", k+1, Chem->GrCharge-i+1);

    if (Chem->GrCharge >= 1)
      sprintf(Chem->Species[n+Chem->GrCharge].name, "gr%d(-)", k+1);

    sprintf(Chem->Species[n+Chem->GrCharge+1].name, "gr%d", k+1);

    if (Chem->GrCharge >= 1)
      sprintf(Chem->Species[n+Chem->GrCharge+2].name, "gr%d(+)", k+1);

    for (i=1; i<Chem->GrCharge; i++)
      sprintf(Chem->Species[n+Chem->GrCharge+2+i].name, "gr%d(%d+)", k+1, i+1);

    for (i=n+1; i<=n+(2*Chem->GrCharge+1); i++)
    {
      Chem->Species[i].mass = Chem->Elements[Chem->N_Ele+k].mass;
      Chem->Species[i].charge = i-n-Chem->GrCharge-1;
      Chem->Species[n].type   = 3;

      Chem->Species[i].numelem = 1;
      for (j=0; j<Chem->N_Ele_tot; j++)
        Chem->Species[i].composition[j] = 0;

      Chem->Species[i].composition[Chem->N_Ele+k] = 1;
      Chem->Species[i].gsize = Chem->GrSize[k];
      Chem->Species[i].Eb = 0.0;     /* Not used */
    }

    Chem->Elements[Chem->N_Ele+k].numsig = 1;

    /* neutral grain as the single element species */
    Chem->Elements[Chem->N_Ele+k].single[0] = n+Chem->GrCharge+1;
  }

  /* Construct mantle species */
  p = p + Chem->NGrain*(2*Chem->GrCharge+1);
	for (i=1;i<=Chem->N_Neu_f + Chem->N_Neu + Chem->N_Neu_s; i++)
  {
		/* Find counterpart gas phase species index */
		if(i<=Chem->N_Neu_f)
			TmpInd = i;
		else if (i<=Chem->N_Neu)
			TmpInd = 2*Chem->N_Neu_f + i;
		else
			TmpInd = 2*Chem->N_Neu_f + Chem->N_Neu + i;

    for (k=0; k<Chem->NGrain; k++)
    {
      n = i+p+k*(Chem->N_Neu_f+Chem->N_Neu+Chem->N_Neu_s);

			sprintf(mantle, "[m%d]", k+1);
			strcpy(Chem->Species[n].name, Chem->Species[TmpInd].name);
			strcat(Chem->Species[n].name, mantle);

			Chem->Species[n].mass = Chem->Species[TmpInd].mass;
			Chem->Species[n].charge = Chem->Species[TmpInd].charge;
			Chem->Species[n].type   = 4;

			Chem->Species[n].numelem = Chem->Species[TmpInd].numelem;
			for (j=0; j<Chem->N_Ele_tot; j++)
				 Chem->Species[n].composition[j] = Chem->Species[TmpInd].composition[j];

			Chem->Species[n].Eb = Chem->Species[TmpInd].Eb;
			Chem->Species[n].gsize = (Real)(-k-1); // N.B.: indicator of associated grain type 
    }
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

  return -10;   /* not found! */
}


/*============================================================================*/
/*------------------------------ PRIVATE FUNCTIONS ---------------------------*/

/*----------------------------------------------------------------------------*/
/* Analyze the composition of a species i and do the related operations
 */
void Analyze(Chemistry *Chem, int i)
{/* Warning: This function ONLY analyze chemical species, no grain! */

  int f, N_Neu, j, k, n, l=0;
  int p = 0;
  Real mass = 0.0;

/* Find out the compositon of the species */
  for (j=0; j<Chem->N_Ele_tot; j++)
    Chem->Species[i].composition[j] = 0;

  while ((p<(NL_SPE-1)) && (Chem->Species[i].name[p] != '\0')
                        && (Chem->Species[i].name[p] != '+')
                        && (Chem->Species[i].name[p] != '-'))
  {
		/* k is element index, l is species length*/
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
  else if (Chem->Species[i].name[p] == '-')
    Chem->Species[i].charge = -1;
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

/* If it is a neutral single-element species, tell the array Chem->Elements[] */
  if ((k == 1) && (Chem->Species[i].charge == 0))
  {
    n = Chem->Elements[l].numsig;
		// numsig is number of its single elements
    Chem->Elements[l].numsig += 1;         //+Chem->NGrain;
    Chem->Elements[l].single[n] = i;      // Neutral single element
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
  ath_pout(0,"Total Number of Grain Types:       %d\n",Chem->NGrain);
  ath_pout(0,"Max Number of Grain Charge:        %d\n\n",Chem->GrCharge);

  ath_pout(0,"List of elements and their single-element species:\n");
  for (i=0; i<Chem->N_Ele; i++)
  {
    ath_pout(0,"%7s:", Chem->Elements[i].name);
    for (j=0; j<Chem->Elements[i].numsig; j++)
    {
      k = Chem->Elements[i].single[j];
      ath_pout(0,"  %5s,",Chem->Species[k].name);
    }
    ath_pout(0,"\n");
  }

  ath_pout(0,"Number of full Neutral Species:    %d\n",Chem->N_Neu_f);
  ath_pout(0,"Number of Normal Neutral Species:  %d\n",Chem->N_Neu);
  ath_pout(0,"Number of Special Neutral Species: %d\n",Chem->N_Neu_s);
  ath_pout(0,"Number of Special Ion Species:     %d\n",Chem->N_Ion_s);
  ath_pout(0,"Total Number of Species:           %d\n\n",Chem->Ntot);

  if (Chem->N_Neu_f > 0)
    ath_pout(0,"Starting index for Normal Neutrals: %d:\n", Chem->NeuInd);

  if (Chem->N_Neu_s > 0)
    ath_pout(0,"Starting index for Special Neutrals: %d:\n", Chem->SNeuInd);

  if (Chem->N_Ion_s > 0)
    ath_pout(0,"Starting index for Special Ions:     %d:\n", Chem->SIonInd);

  if (Chem->NGrain > 0)
    ath_pout(0,"Starting index for Grains:           %d:\n", Chem->GrInd);

  ath_pout(0,"\n");
  ath_pout(0,"List of species:\n");
  ath_pout(0,"\n");
  for (i=0; i<Chem->Ntot; i++)
  {
    sp = &(Chem->Species[i]);
    ath_pout(0,"%4d. %7s: M=%9.2e m_p, Q=%2de, N_Ele=%2d, Eb=%6.1f", i, sp->name,
                                sp->mass, sp->charge, sp->numelem, sp->Eb);
    ath_pout(0,"  Comp:");
    for (j=0; j<Chem->N_Ele_tot; j++)
    {
      k = sp->composition[j];
      if (k > 0)
        ath_pout(0,"  %3s:%2d",Chem->Elements[j].name, k);
    }
    ath_pout(0,"\n");
  }
  ath_pout(0,"\n");
 for(i=0;i<Chem->Ntot;i++)
	ath_pout(0,"%s\n",Chem->Species[i].name);

  return;
}

#endif /* CHEMISTRY */
