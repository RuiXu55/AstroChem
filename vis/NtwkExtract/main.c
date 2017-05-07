/****************************************************************************
 * Code to calculate chemical reactions
 * Last updated November 2010
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "chemistry/chemistry.h"
#include "chemistry/prototypes.h"
#include "prototypes.h"

int main(int argc, char *argv[])
{
  int i, k;
  char id[10];
  char *athinput = NULL;
  Chemistry Chem;

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
/* extract reactions */

  init_species(&Chem);

  reprocess_list(athinput);

  extract_reactions(&Chem);

  final_chemistry(&Chem);

  return EXIT_SUCCESS;
}

