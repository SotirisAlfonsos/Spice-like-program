/*		Circuit Simulation 2014 3 Solution of MNA System
 * 
 * Onoma programmatisti:	Alfonsos Swthrhs 959
 * 				Servou Katerina 852
 * 				Stauropoulos Antonis 1141
 * 				Vakerlis Christos 982
 * 
 * Hmeromhnia:		01/11/2014
 * 
 * Perigrafi programmatos:
 *	To programma exei ws skopo na diavasei ena netlist DC analyshs, toy opoiou to onoma dinetai ws deytero orisma sthn klhsh,
 *	na eisagei ta stoixeia se listes kai na ta ektypwsei.
 * 
 * Sinthikes/ipotheseis:	
 * 	To programma thewrei oti to netlist pou pairnei ws eisodo einai swsto kai opoiadhpote grammh pou den ksekinaei
 * 	apo V,I,R,L,C,D,M,Q den thn lamvanei ypopsin tis
 */

#define _XOPEN_SOURCE 500 /* Enable certain library functions (strdup) on linux.  See feature_test_macros(7) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Csparse/csparse/csparse.c"
#include "gnuplot_i/src/gnuplot_i.h"
#include "structs.h"
#include "functions.c"

//#include "transient_analysis.c"


int main (int argc, char* argv[]) {
  
	char line[MAX_SIZE_GRAMMIS];
	int counter;
	int flag_cholesky=0;
	int flag_biCG = 0;
	int flag_CG = 0;
	int flag_defaultItol=1;
	double ITOL=0.001;
	int flag_sparse = 0;
	int lines=0;
	double trans_h_fin[]={0,0};	//trans_h_fin[0] gia h , trans_h_fin[1] gia fin
	int trans_method=0;		//0 gia TRAPEZOIDAL , 1 gia backward euler
	
	tElements *head_elem = listElem_init();	//dhlwsh-arxikopoihsh listwn
	tDiodes *head_diod = listDiodes_init();
	tMOS *head_mos = listMOS_init();
	tBJT *head_bjt = listBJT_init();
	tPlot *head_plot = listPlot_init();
	dc_analysis dc_table={ "Nothing\0",0,0,0};
			
	tElements *curElem;			//voithitikoi deiktes gia prospelasi twn listwn
	tDiodes *curDiod;
	tMOS *curMos;
	tBJT *curBjt;	
	tPlot *curPlot;
	hashtable_t *hashtable;
	
	    //Anoigw to arxeio. An den exei anoiksei termatizei to programma.
	FILE * fp=fopen(argv[1],"r");	
	if(fp==NULL){	
		printf("To arxeio eisodou den anoikse\n");
		exit(1);
	}

	    //Diavazw to arxeio grammi-grammi kai se kathe epanalhpsh diavazw kai ena stoixeio analoga me to gramma pou arxizei. 		
	    //Spaw th grammh analoga me ta kena gia na parw tis times pou thelw.
	printf("/*** parsing ***/\n");
	while(fgets(line,MAX_SIZE_GRAMMIS,fp)!=NULL){	
		counter=1;
		lines++;
		char *token;
   		token = strtok(line, " ");
		if (token[0]=='I' || token[0]=='V' || token[0]=='i' || token[0]=='v' ) {
			while ( token != NULL ) {
			  if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
			  
			  if (token[0]!=NULL || token[0]!='\0')  {
			      head_elem = insertElem(head_elem, counter, token);
			      counter++;
			  }
			  else
			    printf("parse %s \n",token);
			  
			  token = strtok(NULL, "() \n");
			  
   			}
		}
		else if (token[0]=='R' || token[0]=='C' || token[0]=='L' || token[0]=='r' || token[0]=='c' || token[0]=='l') {
			while ( token != NULL ) {
			  if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
			  head_elem = insertElem(head_elem, counter, token);
			  counter++;
			  token = strtok(NULL, " \n");
   			}
		}		
		else if ( token[0]=='D' || token[0]=='d') {
			while ( token != NULL ) {
			  if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
			  head_diod = insertDiodes(head_diod, counter, token);
			  counter++;
			  token = strtok(NULL, " \n");
   			}
   			if (counter==5) 				//se periptwsi pou sto netlist den uparxei area, dinetai i timi 1
			    head_diod = insertDiodes(head_diod, 5, "1");
		}
		else if (token[0]=='M' || token[0]=='m') {
			while ( token != NULL ) {
			  if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
			  head_mos = insertMOS(head_mos, counter, token);
			  counter++;
			  token = strtok(NULL, " \n");
   			}
		}
		else if (token[0]=='Q' || token[0]=='q') {
			while ( token != NULL ) { 
				if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
				head_bjt = insertBJT(head_bjt, counter, token);
				counter++;
				token = strtok(NULL, " \n");
			}
			if (counter==6) 			   //se periptwsi pou sto netlist den uparxei area, dinetai i timi 1
			    head_bjt = insertBJT(head_bjt, 6, "1");
		}
		else if ((strcmp(token,".OPTIONS")==0) || (strcmp(token,".options")==0)) {
			token = strtok(NULL, " \n");
			while ( token!= NULL ) {
			    if (isspace(token[strlen(token)-1])!=0) {
			      token[strlen(token)-1] = '\0';
			    }
			    
			    if((strcmp(token,"METHOD")==0) || (strcmp(token,"method")==0)) {
			      trans_method=0;		//TRAPEZOIDAL RULE
			    }
			    else if( ((strcmp(token,"BE")==0) || (strcmp(token,"be")==0)) && (counter==2) ){ 
			      trans_method=1;		//backward euler rule
			    }
			    else if ((strcmp(token,"SPD")==0) || (strcmp(token,"spd")==0)) {
				flag_cholesky = 1; 
			    }else if (  (flag_cholesky==1) && ((strcmp(token,"ITER")==0) || (strcmp(token,"iter")==0))  ) {
				flag_CG = 1;
			    }else if (  (flag_cholesky==0) && ((strcmp(token,"ITER")==0) || (strcmp(token,"iter")==0))  ) {
				flag_biCG = 1;
			    }else if ((strcmp(token,"ITOL")==0) || (strcmp(token,"itol")==0)){
				flag_defaultItol=0; 
			    }else if ((strcmp(token,"SPARSE")==0) || (strcmp(token,"sparse")==0)){
				 flag_sparse=1;
			    }else if (counter==2 && flag_defaultItol==0){
				ITOL = atof(token);
			    }

			    counter++;
			    token = strtok(NULL, " \n");
   			}	
		}
		else if ( (strcmp(token,".PLOT")==0) || (strcmp(token,".plot")==0) || (strcmp(token,".PRINT")==0) || (strcmp(token,".print")==0) ){		//create plot list
			token = strtok(NULL, " \n");
			while ( token != NULL ) { 
				if(token[1]=='('){
				  if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
				  //printf("olo: %s \tisaxthes: %s\n",token,&token[2]);
				  token[strlen(token)-1] = '\0';
				  
				  head_plot = insertPlot(head_plot, &token[2],token[0]);
				}
				token = strtok(NULL, " \n");
			}

		}
		else if ((strcmp(token,".DC")==0) || (strcmp(token,".dc")==0)){
			token = strtok(NULL, " \n");
			while ( token != NULL ) { 
				if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
				
				if (counter==1)
				    dc_table.input_var = strdup(token);
				else if (counter==2)
				    dc_table.start_value = atof(token);
				else if (counter==3)
				    dc_table.end_value = atof(token);
				else if (counter==4)
				    dc_table.increment = atof(token);
				   
				counter++;
				token = strtok(NULL, " \n");
			}
		}
		else if( (strcmp(token,".TRAN")==0) || (strcmp(token,".tran")==0) ) {
		      token = strtok(NULL, " \n");
		      while ( token != NULL ) { 
			      if (isspace(token[strlen(token)-1])!=0)  token[strlen(token)-1] = '\0';
			      
			      if (counter==1){
				  trans_h_fin[0]=atof(token);
				  printf("trans_h_fin[0]: %g:\n",trans_h_fin[0]);
			      }  
			      else if (counter==2){
				  trans_h_fin[1]=atof(token);
			      	  printf("trans_h_fin[1]: %g:\n",trans_h_fin[1]);
				
			      }
			      counter++;
			      token = strtok(NULL, " \n");
		      }
		}
		

	}
	
		// fclose , print , free 
		
	fclose(fp);	

	    //ektupwsi twn stoixeiwn twn listwn
	if(lines<1000){
	    for (curElem = head_elem; curElem!=NULL; curElem=curElem->next) {  
	      if (curElem->value > 0.000001)
		printf("%c %s %s %s %lf \n", curElem->type, curElem->name, curElem->posNode, curElem->negNode, curElem->value);
	      else
		printf("%c %s %s %s %.18lf \n", curElem->type, curElem->name, curElem->posNode, curElem->negNode, curElem->value);
	    }
	    for (curDiod = head_diod; curDiod!=NULL; curDiod=curDiod->next) {  
	      if (curDiod->area > 0.000001)
		printf("%s %s %s %s %lf \n", curDiod->name, curDiod->posNode, curDiod->negNode, curDiod->modName, curDiod->area);
	      else
		printf("%s %s %s %s %.18lf \n", curDiod->name, curDiod->posNode, curDiod->negNode, curDiod->modName, curDiod->area);
	    }
	    for (curMos = head_mos; curMos!=NULL; curMos=curMos->next) {  
	      if (curMos->L > 0.000001 && curMos->W > 0.000001)
		printf("%s %s %s %s %s %s %lf %lf \n", curMos->name, curMos->D, curMos->G, curMos->S, curMos->B, curMos->modName, curMos->L, curMos->W);
	      else
		printf("%s %s %s %s %s %s %.18lf %.18lf \n", curMos->name, curMos->D, curMos->G, curMos->S, curMos->B, curMos->modName, curMos->L, curMos->W);
	    }
	    for (curBjt = head_bjt; curBjt!=NULL; curBjt=curBjt->next) { 
	      if (curBjt->area > 0.000001)
		printf("%s %s %s %s %s %lf \n", curBjt->name, curBjt->C, curBjt->B, curBjt->E, curBjt->modName, curBjt->area); 
	      else
		printf("%s %s %s %s %s %.18lf \n", curBjt->name, curBjt->C, curBjt->B, curBjt->E, curBjt->modName, curBjt->area);    
	    }
	    for (curPlot = head_plot; curPlot!=NULL; curPlot=curPlot->next) {
	      if (curPlot == head_plot)
		  printf("\nPLOT:");
	      printf("\t%c(%s)",curPlot->type,curPlot->node);
	    }
	}
	printf("\nflag_cholesky: %d\n",flag_cholesky);
	printf("DC: %s %lf %lf %lf\n",dc_table.input_var, dc_table.start_value, dc_table.end_value, dc_table.increment);
	
	
	    // dimiourgia hashtable , dimiourgia kai epilisi Ax=b sistimatos
	hashtable = sumNodes(head_elem, flag_cholesky, &dc_table, head_plot, flag_biCG, flag_CG, ITOL, flag_sparse, trans_h_fin, trans_method,lines);
	
	
	    //apodesmeusi tis desmeumenis mnimis

	for (curDiod = head_diod; curDiod!=NULL; curDiod=head_diod) {  
	  head_diod = head_diod->next;
	  free(curDiod);
	}  
	for (curMos = head_mos; curMos!=NULL; curMos=head_mos) {  
	  head_mos = head_mos->next;
	  free(curMos);
	}  
	for (curBjt = head_bjt; curBjt!=NULL; curBjt=head_bjt) { 
	  head_bjt = head_bjt->next;
	  free(curBjt);
	}  
	
	return 0;
}
	//    ./cs_parser NetList_Cholesky_Transient.txt
	//     gcc cs_parser.c -o cs_parser
	//     gcc cs_parser.c -o cs_parser -lgsl -lgslcblas -lm
	//     gcc cs_parser.c -o cs_parser -lgsl -lgslcblas -lm gnuplot_i/gnuplot_i.o
	
	/* 
	 * TIPOTA      			-> LU 
	 * .OPTIONS SPD			-> cholesky
	 * .OPTIONS ITER		-> Bi_CG
	 * .OPTIONS SPD ITER 		-> CG
	 * .OPTIONS SPARSE 		-> LU sparse
	 * .OPTIONS SPD SPARSE		-> cholesky sparse
	 * .OPTIONS ITER SPARSE		-> Bi_Cg sparse
	 * .OPTIONS ITER SPARSE 		-> CG sparse
	 * 
	 * .TRAN h fin
	 * meta apo piges:
	 * EXP (i1 i2 td1 tc1 td2 tc20)
	 * SIN (i1 ia fr td df ph)
	 * PULSE (i1 i2 td tr tf pw per)
	 * PWL (t1 i1) ... (tn in)
	 * 
	 * .OPTIONS ITOL 0.6
	   .OPTIONS  spd iter sparse
	   .DC v3bb 1.8 2.8 1
	   .PLOT v(n12_15316200_16398000)  v(n14_24316200_27396000) v(n12_15316200_16398000)  v(n14_24316200_27396000) v(n12_15316200_16398000)  v(n14_24316200_27396000) 
	   .DC
	  *.TRAN 0.1 10
	  *.OPTIONS METHOD be
	  *.OPTIONS SPARSE
	 * 
	 * 
	 * 
	 */ 

	