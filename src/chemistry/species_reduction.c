#include "../header/copyright.h"
/*=============================================================================
 * FILE: species_reduction.c
 *
 * PURPOSE: Contains functions to reduce the full chemical network
 *          using species based reduction method.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   species_reduction.c - reduct chemical network
 *
 * REFERENCES:
 *   D. Wiebe, et al. 2003, A&A, 399, 197-210
 *
 * History:
 *   Written by  Rui Xu		Oct. 2014
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

void species_reduction(ChemEvln *Evln)
{

  int i,j,k,l,m,iter,label,p,q,*add;
  Real rate,MaxB,*sens,*df,*des,*form,*maxrat;
  EquationTerm *EqTerm;
  Chemistry *Chem = Evln->Chem;
  int Ntot= Chem->Ntot;
  FILE *fp; char fname[50];

  Real minsens = par_getd("problem","minsens");
  int maxiter  = par_getd("problem","maxiter");
  Real minratio= par_getd("problem","minratio");
  int red      = par_getd("problem","red");

  /* initialization */
  sens = (Real*)calloc_1d_array(Ntot,sizeof(Real));
  add = (int*)calloc_1d_array(Ntot,sizeof(int));
  df = (Real*)calloc_1d_array(Ntot,sizeof(Real));
  des = (Real*)calloc_1d_array(Ntot,sizeof(Real));
  form = (Real*)calloc_1d_array(Ntot,sizeof(Real));
  maxrat = (Real*)calloc_1d_array(Ntot,sizeof(Real));
  for(i=0;i<Ntot;i++){
    sens[i]=add[i]=df[i]=des[i]=form[i]=maxrat[i]= 0;
  }

  /* first choose important species */
  add[0] = 1;        //e-

  /* calculate the maxrat(j) for all species in the network */
  /* ignore desorption since its reaction rate are far more
   * larger than other reactions */
  for(i=0;i<Ntot;i++){
    for (k=0;k<Chem->Equations[i].NTerm;k++)
    {
      EqTerm = &(Chem->Equations[i].EqTerm[k]);
      if(Chem->Reactions[EqTerm->ind].rtype !=4){//exclude desorption
        /* destruction channel */
        if(EqTerm->dir<0){
          rate = Evln->K[EqTerm->ind];
          for (l=0;l<EqTerm->N;l++){
            p = EqTerm->lab[l];
            rate *= Evln->NumDen[p];
          }
          des[i] += rate;
        }
        /* formation channel */
        if(EqTerm->dir>0){
          rate = Evln->K[EqTerm->ind];
          for (l=0;l<EqTerm->N;l++){
            p = EqTerm->lab[l];
            rate *= Evln->NumDen[p];
          }
          form[i] += rate;
        }
      }// end rtype
    } // end Nterm
    maxrat[i] = MAX(des[i],form[i]);
  }//end of calculate maxrat[i]

  /* Add one species each iteration */
  iter = 0;
  while(iter<maxiter){
    ath_pout(0,"iteration=%d\n",iter);
    /* For each species, sum over the 
     * influence of all other species */
    for(i=0;i<Ntot;i++)/*loop over all species*/
    {
      sens[i] = 0.0;
      for(j=0;j<Ntot;j++){
      if(add[j]== 1){
        for (k=0;k<Chem->Equations[j].NTerm;k++){
          EqTerm = &(Chem->Equations[j].EqTerm[k]);
          if(Chem->Reactions[EqTerm->ind].rtype !=4){//exclude desorption
            /* compute the value of df/dn */
            /* One reactant */
            if(EqTerm->N==1){
              if(EqTerm->lab[0] == i){
                df[j] = Evln->K[EqTerm->ind]*EqTerm->dir;
                sens[i] +=pow(Evln->NumDen[i]*df[j]/maxrat[j],2);
            }}
            else if(EqTerm->N ==2){ /*Two reactant*/
              if(EqTerm->lab[0] ==i){
                q = EqTerm->lab[1];
                df[j] = Evln->K[EqTerm->ind]*Evln->NumDen[q]*EqTerm->dir;
                sens[i] += pow(Evln->NumDen[i]*df[j]/maxrat[j],2);
              }
              if(EqTerm->lab[1] ==i){
                q = EqTerm->lab[0];
                df[j] = Evln->K[EqTerm->ind]*Evln->NumDen[q]*EqTerm->dir;
                sens[i] += pow(Evln->NumDen[i]*df[j]/maxrat[j],2);
              }}
          }
       }}
    }} //loop for k && j && i;

    /* estimate sensi and determine 
     * which species need to be added */
    Real minsensin = 1e3;
    int minlabel;
    for(k=0;k<Ntot;k++){
      if(sens[k]<minsensin && add[k]>0.0){
         minsensin=sens[k] ; minlabel = k;
      }
    }
    if(red<1){ /* add one species one time*/ 
      MaxB = -1.0;
      for(k=0;k<Ntot;k++){
        if(sens[k]>MaxB && add[k]<1.0){
          MaxB = sens[k]; label = k;
        }
      }
      /* end of iteration when maxB<minsens*/
      if(MaxB/minsensin<minratio){
        ath_pout(0,"minin=%e, max=%e\n",minsensin,MaxB);
        ath_pout(0,"NO more important species!");
        break;
      }
      add[label] = 1.0;
      ath_pout(0,"Add species: %7s sens: %e\n",Chem->Species[label].name,sens[label]);
      iter += 1.;
    }
    else{
    /* add species larger than threshold */
      int  addnum = 0;
      for(k=0;k<Ntot;k++){
        if(sens[k]>minsens && add[k]<1.0){
          add[k] = 1.;
          addnum += 1;
          ath_pout(0,"Add species: %7s sens: %e\n",Chem->Species[k].name,sens[k]);
        }
      }
      iter += 1.;
      if (addnum<1)
        break;
    }

  }//while loop.


  /* output the all the species added in the network */
  sprintf(fname,"%s",par_gets("job","savep"));
  if((fp = fopen(fname,"w")) == NULL)
    ath_error("Error open file %s \n",fname);
  int out = par_getd("problem","outform");
  switch(out){
    case 0:
    for(i=0;i<Ntot;i++){
      if(add[i]>0)
      fprintf(fp,"%e %7s %d\n",sens[i],Chem->Species[i].name,add[i]);
    }
    break;
    default: 
    for(i=0;i<Ntot;i++){
      MaxB=0.0;
      for(j=1;j<Ntot;j++){
       if(sens[j]>MaxB){
         MaxB = sens[j];k=j;}
      }
      fprintf(fp,"%e  %7s %d \n",sens[k],Chem->Species[k].name,add[k]);
      sens[k]=0 ;
     }
    break;
  }
  fclose(fp);
} // end 



