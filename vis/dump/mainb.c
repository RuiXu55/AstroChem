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
  int i,j,k,nt,verbose;
  char id[20];
  char *athinput = NULL;
  Real tend,dttry,atol,rho, Tg, r,zeta_eff,zs,ze;
  Real G;
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
  printf("begin init chemistry");
  init_chemistry(&Chem);
  printf("begin init ChemEvln");
  init_chemevln (&Chem, &Evln);

/*--- Step 4. ----------------------------------------------------------------*/
/* main calculation */

r     = par_getd("problem","r");
tend  = par_getd("problem","te")*OneYear;
dttry = par_getd("problem","dt0") * OneYear; /* initial trial time step */
atol  = par_getd("problem","atol");
zs    = par_getd("problem","zstart");
ze    = par_getd("problem","zend");
G     = par_getd("disk","G");


/* Disk property */
init_disk(&Disk); 
init_chemout(&Chem,&ChemOut,1,"R",r,"0");

for(k=zs;k<ze;k++){
  zeta_eff = Ionization_disk(&Disk,r,k/10.);
  Tg = Temp_disk(&Disk,r);	  // the temperature at 1AU
  ath_pout(0,"Tg = %e \n",Tg);
  rho = Rho_disk(&Disk,r,k/10.0);	 //radius + height
  sprintf(id,"%d",k);
  verbose = k;
  /* choose all species for the output */
  ChemSet_allspecies(&Chem, &ChemOut);
  /* initialize the number density with single-element species */
  //if(k==0){
  init_numberden(&Evln, rho, verbose);
  //}
  //else{
  //  reset_numberden(&Evln,rho,verbose);
  //}

  /* calculate the rate coefficients for all reactions */
  //IonizationCoeff(&Evln, zeta_eff, 0.0, verbose); /* ionization reations */
  /* Ionization with G */
  //IonizationCoeff1(&Evln, zeta_eff, 0.0,G, verbose);

  CalCoeff       (&Evln, Tg, verbose);       /* all other reactions */

  /* evolve the network from 0 to tend */
  Evln.t = 0.0;
  evolve(tend, dttry,atol);
  /* chemical network reduction */
  //species_reduction(&Evln);
  select_reaction(&Evln);

  /* output the number densities */
  output_nspecies(&Evln, &ChemOut, "z", k/10.0);	
  //output_etaB(&Evln, &ChemOut, "rho",rho,rho,nB);
}

/*--- Step 5. ----------------------------------------------------------------*/
/* finalization */
final_chemout(&ChemOut);
final_chemevln (&Evln);
final_chemistry(&Chem);
par_close();
return EXIT_SUCCESS;

}

