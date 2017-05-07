#include "../header/copyright.h"
/*=============================================================================
 *FILE: reduction.c
 * PURPOSE: Contains functions to reduce the full chemical network
 * CONTAINS PUBLIC FUNCTIONS:
 *  select_reaction.c - reduct chemical network
 * REFERENCES:
 *   D. Wiebe, et al. 2003, A&A, 399, 197-210
 * History:
 *   Written by  Rui Xu         Sept. 2014
 ==============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "../header/defs.h"
#include "../header/chemistry.h"
#include "../header/prototypes.h"
#include "../header/chemproto.h"
          
/*============================================================================
 * private function:
   public function:
   select_reaction(ChemEvln *Evln)
   choose the most important formation and destruction channel for all
   species that can compensate 90% abundance of the species.
 =========================================================================*/

void select_reaction(ChemEvln *Evln)
{
  int i,j,k,l,m,p,num,ind;
  Real rate,*des,*form;
  EquationTerm *EqTerm;
  Chemistry *Chem = Evln->Chem;
  int Ntot = Chem->Ntot;
  FILE *fp;
  char fname[20];

  Real limit = par_getd("problem","limit");
  printf("limit: %e\n ",limit);
  des   = (Real*)calloc_1d_array((Ntot),sizeof(Real));
  form  = (Real*)calloc_1d_array((Ntot),sizeof(Real));
  for(i=0;i<Ntot;i++){
    des[i] = form[i] = 0.;
  }

  /* calculate the G[i] and L[i]  for all species in the network */
  for(i=0;i<Ntot;i++){
    for (k=0;k<Chem->Equations[i].NTerm;k++)
    {
      EqTerm = &(Chem->Equations[i].EqTerm[k]);
      /* destruction channel */
      if(EqTerm->dir<0){
        rate = Evln->K[EqTerm->ind];
        for (l=0;l<EqTerm->N;l++)
        {
          p = EqTerm->lab[l];
          rate *= Evln->NumDen[p];
        }
        des [i] += rate;
      }
      /* formation channel */
      if(EqTerm->dir>0){
        rate = Evln->K[EqTerm->ind];
        for (l=0;l<EqTerm->N;l++)
        {
          p = EqTerm->lab[l];
          rate *= Evln->NumDen[p];
        }
        form[i] += rate;
      }
    }/* end iter over reactions for same species*/
  }/* iteration over species */

  sprintf(fname,"%s",par_gets("job","saver"));
  fp = fopen(fname,"w");
  for (i=0;i<Ntot;i++)
  {	
    //printf("%d th species: %7s\n ",i,Chem->Species[i].name);
    for (k=0;k<Chem->Equations[i].NTerm;k++)
    {
      EqTerm = &(Chem->Equations[i].EqTerm[k]);
      /* absolute value for all the reactions */
      rate = Evln->K[EqTerm->ind];
      for (l=0;l<EqTerm->N;l++)
      {
        p = EqTerm->lab[l];
        rate *= Evln->NumDen[p];
      }
      /* ith species is a product
       * formation reaction */
      if(EqTerm->dir>0)
      {
        /* print reactions if the reaction 
         * rate is larger than threshold*/
        printf("rate: %e form: %e limit:%e\n ",rate,form[i],limit);
        if(rate/form[i]>limit){
          printf("ratio:%e limit:%e\n",rate/form[i],limit);
          ind   = EqTerm->ind;
          /* first reactant */
          p = Chem->Reactions[ind].reactant[0];
          fprintf(fp,"%7s ",Chem->Species[p].name);
          /* second reactant */
          p = Chem->Reactions[ind].reactant[1];
          if(p>=0)
            fprintf(fp,"%7s ",Chem->Species[p].name);
          else
            fprintf(fp,"%7s ", "0");
          /* first product */
          p = Chem->Reactions[ind].product[0];
          fprintf(fp,"%7s ",Chem->Species[p].name);
          /* second product */
          p = Chem->Reactions[ind].product[1];
          if(p>=0)
            fprintf(fp,"%7s ", Chem->Species[p].name);
          else
            fprintf(fp,"%7s ","0");
          /* third product */
          p = Chem->Reactions[ind].product[2];
          if(p>=0)
            fprintf(fp,"%7s ", Chem->Species[p].name);
          else
            fprintf(fp,"%7s ","0");
          /* fourth product */
          p = Chem->Reactions[ind].product[3];
          if(p>=0)
            fprintf(fp,"%7s ",Chem->Species[p].name);
          else
            fprintf(fp,"%7s","0");
          fprintf(fp,"\n");
        } /* end print reactions */
      }else{ 
      /* ith species is a reactant/destruction reaction
       * print reactions if the reaction 
       * rate is larger than threshold */
      printf("rate: %e des: %e limit:%e\n ",rate,des[i],limit);
      if(rate/des[i]>limit){
        printf("ratio:%e limit:%e\n",rate/des[i],limit);
        ind   = EqTerm->ind;
        /* first reactant */
        p = Chem->Reactions[ind].reactant[0];
        fprintf(fp,"%7s ",Chem->Species[p].name);
        /* second reactant */
        p = Chem->Reactions[ind].reactant[1];
        if(p>=0)
          fprintf(fp,"%7s ",Chem->Species[p].name);
        else
          fprintf(fp,"%7s ", "0");
        /* first product */
        p = Chem->Reactions[ind].product[0];
        fprintf(fp,"%7s ",Chem->Species[p].name);
        /* second product */
        p = Chem->Reactions[ind].product[1];
        if(p>=0)
          fprintf(fp,"%7s ", Chem->Species[p].name);
        else
          fprintf(fp,"%7s ","0");
        /* third product */
        p = Chem->Reactions[ind].product[2];
        if(p>=0)
          fprintf(fp,"%7s ", Chem->Species[p].name);
        else
          fprintf(fp,"%7s ","0");
        /* fourth product */
        p = Chem->Reactions[ind].product[3];
        if(p>=0)
          fprintf(fp,"%7s ",Chem->Species[p].name);
        else
          fprintf(fp,"%7s","0");
        fprintf(fp,"\n");
       } /*end if limit*/
     } /* end if reactant/product */
    }/* end k */
  }/* end i */
  printf(" Output file completed!\n");
  return ;
}








