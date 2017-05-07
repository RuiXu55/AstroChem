#include "../copyright.h"
/*=============================================================================
 * FILE: reprocess_list.c
 *
 * PURPOSE: Initialize all chemical reactions. It reads ionization reactions
 *   and other chemical reactions from standard input file "ChemReactions.txt"
 *   in the same directory as the executable. Below is a sample input format:
 *
==============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../defs.h"
#include "chemistry.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef CHEMISTRY

/*============================================================================*/
/*----------------------------- Public Functions -----------------------------*/

/*----------------------------------------------------------------------------*/
/* Second step of initializing the chemical model
 * Read all chemical reactions, construct all grain related reactions
 */
void reprocess_list(char *list)
{
  FILE *frp,*fwp;
  int i, j, k, n, len;
  char line[MAXLEN],newline[MAXLEN];

  frp = fopen(list,"r");
  fwp = fopen("AllReactions.txt","w");

  if(frp == NULL){
    ath_perr(-1,"[reprocess_list]: Unable to open the input file!\n");
    return;
  }
  if(fwp == NULL){
    ath_perr(-1,"[reprocess_list]: Unable to open the AllReactions.txt!\n");
    return;
  }

  fscanf(frp, "%d\n", &n);	/* Number of reactions */
  fprintf(fwp, "%d\n", n);
  if (n < 0)
    ath_error("[reprocess_list]: Number of reactions must >=0!\n");

  for (i=0; i<n; i++) {

    fgets(line,MAXLEN,frp);
    len = strlen(line);

    j=0; k=0;
    while (k<len-1) {
      if (line[k] == ':'){
        if (line[k+1] == ':')
        {
          newline[j]=' ';
          newline[j+1]='0';
          j += 2;
        }
        else
        {
          newline[j]=' ';
          j += 1;
        }
      }
      else {
       newline[j] = line[k];
       j++;
      }
      k++;
    }
    newline[j]='\0';

    fprintf(fwp,"%s\n",newline);
  }

  fclose(frp);
  fclose(fwp);
	
  return;
}

#endif /* CHEMISTRY */

