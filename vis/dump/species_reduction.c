#include "../copyright.h"
/*=============================================================================
 * FILE: species_reduction.c
 *
 * PURPOSE: Contains functions to reduce the full chemical network
 *          using species based reduction method.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   reduction.c - reduct chemical network
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
#include "../defs.h"
#include "chemistry.h"
#include "prototypes.h"
#include "../prototypes.h"



void species_reduction(ChemEvln *Evln)
{

	int i,j,k,l,label,p,q,num,add_num;
	Real rate,tot_weight;
	EquationTerm *EqTerm;
	Chemistry *Chem = Evln->Chem;
	int Ntot	= Chem->Ntot;
	int NReaction 	= Chem->NReaction;
	FILE *fp;
	char fname[50];
        Real *B, *add_mark, *df, *g, *G, *L;
	B = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	add_mark = (int*)calloc_1d_array((Chem->Ntot),sizeof(int));
	df = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	g = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	G = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	L = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));

	/* initialization */
	for(i=0;i<Ntot;i++){
		B[i] = 0;
		add_mark[i]=0;
		f[i] = 0;
		G[i] = 0;
		L[i] = 0;
	}
	/* first we only choose eletron as important species */
	add_mark[0]=1;

	/* calculate the g(j) for all species in the network */
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
		L[i] += rate;
		}
		/* formation channel */
		if(EqTerm->dir>0){
		rate = Evln->K[EqTerm->ind];
		for (l=0;l<EqTerm->N;l++)
		{
			p = EqTerm->lab[l];
			rate *= Evln->NumDen[p];
		}
		G[i] += rate;
		}
	}
	if(G[i]>L[i]){
		g[i] = G[i];
	}
	else{
		g[i] = L[i];
	}
	}

	while(add_num>0)
	{
	    add_num =0;
	    /* for entire network we compute the sensitivity  B[i] = SIGMA_j[n(i)/g(j)][partial_f(j)/partial_n(i)]*/
	    for(i=0;i<Ntot;i++)
	    {
		f[i] =0;
		for(j=0;j<Ntot;j++){
		   for (k=0;k<Chem->Equations[i].NTerm;k++)
		   EqTerm = &(Chem->Equations[j].EqTerm[k]);
		   /* compute the value of df/dn */
		   if(EqTerm->N==1)
		   {
			if(EqTerm->lab[1] == i)
			{
			    df[j] += Evln->K[EqTerm->ind]*EqTerm->dir
			}
		   }
		   else
		   {
			if(EqTerm->lab[1] ==i)
			{
			    q = EqTerm->lab[2];
			    df[j] += Evln->K[EqTerm->ind]*Evln->NumDen[q]*EqTerm->dir;
			}
		   }

			
		}
		B[i] += power(Evln->NumDen[i]*df[j]/g[j],2);
	    }
	    /* estimate B[i] and determine how many necessary species need to be added */
	    for(i=0;i<Ntot;i++){
		if(add_mark[i]<1){
		   if(B[i]>1){
			add_mark[i] =1;
			add_num +=1;
		    }
		}
	     }		
	}//while loop

	/* output the all the species added in the network */
	sprintf(fname,"%s",par_gets("job","1_SpeciesNames"));
	if((fp = fopen(fname,"w")) == NULL)
		ath_error("Error open file %s \n",fname);
	for(i=0;i<Ntot;i++){
		if(add_mark[i]>0){
		fprintf(fp,"%7s \n",Chem->Species[species_ind[i]].name);
		}
	}
	fclose(fp);
}
