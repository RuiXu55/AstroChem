#include "../copyright.h"
/*=============================================================================
 * FILE: reduction.c
 *
 * PURPOSE: Contains functions to reduce the full chemical network
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   reduction.c - reduct chemical network
 *
 * REFERENCES:
 *   D. Wiebe, et al. 2003, A&A, 399, 197-210
 *
 * History:
 *   Written by  Rui Xu		Sept. 2014
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

/*============================================================================
private function:
	sort_reaction(ChemEvln *Evln, Real *omega_r, Real *omega_s,int num);
	give the first most important reactions and their index.
	print_reaction()
public function:
	reduction(ChemEvln *Evln);	
  ==========================================================================*/
void sort_reaction(ChemEvln *Evln, Real *omega_r, Real *omega_s, int num);
void save_reaction(ChemEvln *Evln, int *sort_ind, Real *sort_val, int num);
void save_species(ChemEvln *Evln, int *species_ind, Real *species_val, int num);




void reduction(ChemEvln *Evln)
{
	int i,j,k,l,label,p,num,iter;
	Real rate,tot_weight;
	EquationTerm *EqTerm;
	Chemistry *Chem = Evln->Chem;
	int Ntot	= Chem->Ntot;
	int NReaction 	= Chem->NReaction;
	FILE *fp,*fp1;
	char species_name[50];
        Real *omega_r, *old_omega_r, *omega_s, *Rtot;
	omega_r = (Real*)calloc_1d_array((Chem->NReaction),sizeof(Real));
	old_omega_r = (Real*)calloc_1d_array((Chem->NReaction),sizeof(Real));
	omega_s = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	Rtot	= (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	for(i=0;i<NReaction;i++){
		omega_r[i]=0;
	}
	for(i=0;i<Ntot;i++){
	omega_s[i] = 0;
	Rtot[i]	=0;
	}
	/* set the weight of electron to be 1, others remain to be zero */
	omega_s[0] 	= 1;
	printf("%s", Chem->Species[146].name);
	iter =1;
	/* update the weight of reactions and species each iteration step */
	while(1)
	{
		/* estimate the total rate containing species j */
		for(i=0;i<Ntot;i++){
		Rtot[i] =0;
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
			Rtot[i] += rate;
		}
		}//for loop

		for (j=0;j<NReaction;j++)
		{
			//Reactant
			omega_r[j]=0;
			rate = Evln->K[j];
			for(k=0;k<2;k++){
			if(Chem->Reactions[j].reactant[k]>=0)
			rate *=Evln->NumDen[Chem->Reactions[j].reactant[k]];
			}
			for(k=0;k<2;k++){
			if(Chem->Reactions[j].reactant[k]>=0)
			omega_r[j] += rate*omega_s[Chem->Reactions[j].reactant[k]]/Rtot[Chem->Reactions[j].reactant[k]];
			}
			//product
			for(k=0;k<4;k++){
			if(Chem->Reactions[j].product[k]>=0)
			omega_r[j] += rate*omega_s[Chem->Reactions[j].product[k]]/Rtot[Chem->Reactions[j].product[k]];
			}
			if(omega_r[j]<old_omega_r[j]){
			omega_r[j] = old_omega_r[j];
			}
		}
		/* normalize reaction weight */
		tot_weight =0;
		for (j=0;j<NReaction;j++){
		tot_weight += omega_r[j];
		}
		for (j=0;j<NReaction;j++){
		omega_r[j]=omega_r[j]/tot_weight;
		}
	/* update the weight of species */
	for (i=0;i<Ntot;i++)
	{
		//omega_s[i] =0;
		for(j=0;j<NReaction;j++)
		{
			if (i==Chem->Reactions[j].reactant[0] || i==Chem->Reactions[j].reactant[1] ||i==Chem->Reactions[j].product[0] || i==Chem->Reactions[j].product[1] || i==Chem->Reactions[j].product[2] || i==Chem->Reactions[j].product[3])
			{
		//		printf("%d\n",i); 
				if(omega_r[j]>omega_s[i])
				omega_s[i] = omega_r[j];
			}
		}
	}
	
	/* normalize the species weight */
	tot_weight =0;
	for (i=0;i<Ntot;i++){
	tot_weight += omega_s[i];
	}
	for(i=0;i<Ntot;i++){
	omega_s[i] = omega_s[i]/tot_weight;
	}

	/* copy the omega_r to the old omega_r */
	for (j=0;j<NReaction;j++){
	old_omega_r[j] = omega_r[j];
	}

	/* end of iteration : when all the species are involved in */
	label =1;
	for (i=0;i<Ntot;i++)
	{
		if(omega_s[i]<1e-100)
		label =0;
	}
	//if(label >0.5) break;
	if(iter>1000) break;
	printf("%d\n",iter);
	//printf("%d %d\n",Ntot,NReaction);
	iter +=1;
	num =0;
	for(i=0;i<NReaction;i++){
	if(omega_r[i]>1e-100)
	num +=1;
	}
	printf("%d\n",num);
	num=0;
	for(i=0;i<Ntot;i++){
	if(omega_s[i]>1e-100)
	num +=1;
	}
	printf("%d\n",num);
	}//while loop

	
	/* output file: output the reduced species and the reaction label*/
	sprintf(species_name, "%s-%s",
                            par_gets("job","outbase1"),par_gets("job","outid1"));
	if ((fp = fopen(species_name,"w")) == NULL)
   		 ath_error("[output_chem]: Error opening file %s...\n", species_name);
	for(j=0;j<NReaction;j++){
	fprintf(fp,"%e ",omega_r[j]);
	}
	/*for(i=0;i<Ntot;i++){
	fprintf(fp,"%8.4f ",omega_s[i]);
	}*/
	fclose(fp);
	/* sort the important reactions and save in 'reduce_reaction' file*/
	sort_reaction(Evln,omega_r,omega_s,NReaction);
	/* finalize and free space */
	free_1d_array(old_omega_r);
	free_1d_array(omega_r);
	free_1d_array(omega_s);
	free_1d_array(Rtot);
   return;
}


void sort_reaction(ChemEvln *Evln, Real *omega_r, Real *omega_s,int num)
{
	int i,j;
	int max_ind;
	char fname[20];
	Real max_val;
	Chemistry *Chem = Evln->Chem;
        int Ntot = Chem->Ntot;
	int NReaction = Chem->NReaction;
	int *sort_ind, *sort_species_ind;
	Real *sort_val, *sort_species_val;
	sort_ind = (int*)calloc_1d_array((Chem->NReaction),sizeof(int));
	sort_val = (Real*)calloc_1d_array((Chem->NReaction),sizeof(Real));
	sort_species_val = (Real*)calloc_1d_array((Chem->Ntot),sizeof(Real));
	sort_species_ind = (int*)calloc_1d_array((Chem->Ntot),sizeof(int));
	FILE *fp;
  for(i=0;i<num;i++)
  {
	max_val = 0;
	max_ind = 0;
	for(j=0;j<NReaction;j++){
		if(omega_r[j]>max_val){
		max_ind=j;
		max_val = omega_r[j];
		}
	}
	sort_val[i] = max_val;
	sort_ind[i] = max_ind;
	omega_r[max_ind] =0;
  }
	/* save the sorted reaction in 'reduce_reaction' file    */
	save_reaction(Evln, sort_ind, sort_val, num);

	/* sort species and save in 'reduce_species' file  */
	for(i=0;i<Chem->Ntot;i++){
	max_val =0 ;
	max_ind =0 ;
	for(j=0;j<Ntot;j++){
		if(omega_s[j]>max_val){
		max_ind =j;
		max_val =omega_s[j];
		}
	}
	sort_species_val[i] = max_val;
	sort_species_ind[i] = max_ind;
	omega_s[max_ind] =0;
	printf("%e %d %7s \n", sort_species_val[i],sort_species_ind[i],Chem->Species[max_ind].name);
	}
	save_species(Evln,sort_species_ind, sort_species_val,Ntot); 

}




void save_reaction(ChemEvln *Evln, int *max_ind, Real *max_val, int num)
{
	int i,j,p;
	char fname[20];
	FILE *fp;
	Chemistry *Chem = Evln->Chem;
	sprintf(fname,"%s",par_gets("job","r_react"));
	if((fp = fopen(fname,"w")) == NULL)
		ath_error("Error open file %s \n",fname);
	for(i=0;i<num;i++){
		fprintf(fp,"%e ",max_val[i]);
		p = Chem->Reactions[max_ind[i]].reactant[0];
		fprintf(fp,"%7s ",Chem->Species[p].name);
		p = Chem->Reactions[max_ind[i]].reactant[1];
		if(p>=0) fprintf(fp,"%7s ",Chem->Species[p].name);
		//product
		p = Chem->Reactions[max_ind[i]].product[0];
		fprintf(fp,"%7s ",Chem->Species[p].name);
		p = Chem->Reactions[max_ind[i]].product[1];
		if(p>=0) fprintf(fp,"%7s ", Chem->Species[p].name);
		p = Chem->Reactions[max_ind[i]].product[2];
		if(p>=0) fprintf(fp,"%7s ", Chem->Species[p].name);
		p = Chem->Reactions[max_ind[i]].product[3];
		if(p>=0) fprintf(fp,"%7s ",Chem->Species[p].name);
		fprintf(fp,"\n");
	}
	fclose(fp);
}


void save_species(ChemEvln *Evln, int *species_ind, Real *species_val, int num){

	int i,j,p;
	char fname[20];
	Chemistry *Chem = Evln->Chem;
	FILE *fp;
	sprintf(fname,"%s",par_gets("job","r_species"));
	if((fp = fopen(fname,"w")) == NULL)
		ath_error("Error open file %s \n",fname);
	for(i=0;i<num;i++){
		fprintf(fp,"%e  %7s \n",species_val[i],Chem->Species[species_ind[i]].name);
	}
	fclose(fp);

}
