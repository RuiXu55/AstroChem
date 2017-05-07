/****************************************************************************
 * Code to calculate chemical reactions
 * Last updated November 2010
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "header/chemistry.h"
#include "header/prototypes.h"
#include "header/chemproto.h"

int main(int argc, char *argv[])
{
  int i,j,k,nt,verbose,niter;
  char id[20];
  char *athinput = NULL;
  Real tend,dttry,atol,rho,Tg,zeta;
  Real sigmax,sigmin,r,nG;
  Real z,sigmaz,G,Tgmin;
  //Chemistry Chem;
  //ChemEvln  Evln;
  ChemOutput ChemOut;
  Nebula Disk;
/*--- Step 1. ----------------------------------------------------------------*/
/* Check for command line options and respond.  See comments in usage()
 * for description of options.  */

  for (i=1; i<argc; i++) {
/* If argv[i] is a 2 character string of the form "-?" then: */
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                                /* -i <file>   */
        /* specify input file */
	athinput = argv[++i];
	break;
      default:
        ath_error(0, "%s -i <file>\n", argv[0]);
        exit(EXIT_FAILURE);
        break;
      }
    }
  }

  /* Print usage message if no input file specified */
  if (athinput == NULL) {
    ath_error(0, "%s -i <file> \n", argv[0]);
    exit(EXIT_FAILURE);
  }

/*--- Step 2. ----------------------------------------------------------------*/
/* Initialization and basic info */

  par_open(athinput);   /* opens AND reads */
  par_cmdline(argc,argv);

  par_dump(0,stdout); /* Dump a copy of the parsed information to athout */

/*--- Step 3. ----------------------------------------------------------------*/
/* initialization */
  ath_pout(0,"begin init chemistry!\n");
  init_chemistry(&Chem);
  ath_pout(0,"begin init ChemEvln\n");
  init_chemevln (&Chem, &Evln);

/*--- Step 4. ----------------------------------------------------------------*/
/* main calculation */

  tend     = par_getd("problem","te")*OneYear;
  dttry    = par_getd("problem","dt0") * OneYear; /* initial trial time step */
  atol     = par_getd("problem","atol");

  r        = par_getd("problem","r");
  sigmax   = par_getd("problem","sigmax");
  sigmin   = par_getd("problem","sigmin");
  Real sig = par_getd("problem","sig");
  Real depth = par_getd("problem","depth");
  niter    = par_getd("problem","niter");

  /* Disk property */
  init_disk(&Disk); 
  /* mode+bname+bvalue+outid */
  init_chemout(&Chem,&ChemOut,1,"R",r,"0");
  
  if(r<2.){
   G = 1e7;
   }
   else if (r<20.){
   G = 1e5;}
   else{
   G = 1e3;
   }

  Tgmin = Temp_disk(&Disk,r);
  for(i=0;i<niter;i++)
  {
    sigmaz = pow(10.,sigmax-i*(sigmax-sigmin)/(niter-1.));
    nG     = G*exp(-pow((sigmaz/depth),2));
    Tg     = Tgmin+i*(10.*Tgmin-Tgmin)/(niter-1.);
    zeta   = Ionization_disk1(&Disk,r,sigmaz);
    rho    = sigmaz/Height_disk(&Disk,r)/AU; 
    if (sigmaz<sig){
    verbose  = 0;
    /* choose all species for the output */
    ChemSet_allspecies(&Chem, &ChemOut);
    /* initialize the number density with single-element species */
    init_numberden(&Evln, rho, verbose);
    /* Ionization with G */
    IonizationCoeff1(&Evln,zeta,0.,nG,verbose);
    CalCoeff       (&Evln, Tg, verbose);       /* all other reactions */
    /* evolve the network from 0 to tend */

    Evln.t = 0.0;
    evolve(tend, dttry,atol);
    /* chemical network reduction */
    //species_reduction(&Evln);
    select_reaction(&Evln);
    /* output the number densities */
    output_nspecies(&Evln, &ChemOut, "G", G);	
    break;
    }
  }
/*--- Step 5. ----------------------------------------------------------------*/
/* finalization */
final_chemout(&ChemOut);
final_chemevln (&Evln);
final_chemistry(&Chem);
par_close();
return EXIT_SUCCESS;

}

