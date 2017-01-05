/*		Circuit Simulation 2014 5 Solution of MNA System
 * 
 * Onoma programmatisti:	Alfonsos Swthrhs 959
 * 				Servou Katerina 852
 * 				Stauropoulos Antonis 1141
 * 				Vakerlis Christos 982
 * 
 * Hmeromhnia:		28/11/2014
 */

int cs_LU (cs *C,double *B, int size, dc_analysis *dc_table, tPlot *head_plot, int thesiDC, tarrayX *X, dcIt dcInodes){								//14/11
     
    int i,j;
    css *S;
    csn *N;
    double *x;
    double inc;
    double *B_new;
    tPlot *cur; 
    FILE *fp;
    
    fp=fopen("results.txt","w");

    S = (css *) cs_calloc(1, sizeof(css)); 
    N = (csn *) cs_calloc(1, sizeof(csn)); 
    
    x = (double *)calloc(size,sizeof(double)); 
    if (x == NULL) {
      printf("Memory allocation failed. Exiting...\n");
      exit(1);
    }
     
    B_new = (double *)calloc(size,sizeof(double)); 
    if (B_new == NULL) {
      printf("Memory allocation failed. Exiting...\n");
      exit(1);
    }

     for (i = 0; i < size; i++){
	B_new[i] =  B[i];
     }
      
    S = cs_sqr(2, C, 0);    //printf("4.1\n");
    N = cs_lu(C, S, 1);    //printf("4.2\n");

 //   cs_spfree(C);    
    
    fprintf(fp,"Oi luseis tis eksiswsis me LU gia araious pinakes einai : \n");
        
    cs_ipvec(N->pinv, B_new, x, size);   // printf("4.5\n"); 
    cs_lsolve(N->L, x);    // printf("5\n");
    cs_usolve(N->U, x);    // printf("6\n");
    cs_ipvec(S->q, x, B_new, size);  
	   
    fprintf(fp,"---------------------------------------------------------------------\n");
    printf("/*** solved LU , starts print and DC sweep (if asked ***/\n");

    for (i=0; i<size; i++) {		
         fprintf(fp, "x[%d] = %g\n", i,B_new[i]);
    }
    
    if (thesiDC != -1){
	for (inc=dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	  
	  B[ thesiDC ] = inc;
	  
	  for (i = 0; i < size; i++){
	      B_new[i] =  B[i];
	  }
	  
	  cs_ipvec(N->pinv, B_new, x, size);
	  cs_lsolve(N->L, x);
	  cs_usolve(N->U, x);
	  cs_ipvec(S->q, x, B_new, size);
	  
	  fprintf(fp,"---------------------------------------------------------------------\n");
	    
	  fprintf(fp,"for %s:  %lf\t", dc_table->input_var, inc);
	  for(cur = head_plot; cur != NULL; cur = cur->next){
	    for(i = 0; i < size; i++){
	      if ( strcmp(X[i].name,cur->node) == 0 ){

		for (j = 0; j < size; j++){
		  if (S->q[j] == i) break;
		}
		  
		fprintf(fp," %c(%3s):  %g ",cur->type, X[i].name, x[j] );
	      }
	    }
	  }
	  fprintf(fp,"\n");
	}
    }else if (dcInodes.node1 != -1){				
	for (inc = dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc += dc_table->increment){
	    
	    if (dcInodes.flag==0){
		B[ dcInodes.node1 ] = -inc;
		if (dcInodes.node2 != -1)
		    B[dcInodes.node2] = inc;	
	    }else if (dcInodes.flag==1) {
		B[dcInodes.node1 ] = inc;
		if (dcInodes.node2 != -1)
		    B[ dcInodes.node2 ] = -inc;
	    }

	    for (i = 0; i < size; i++){
		B_new[i] =  B[i];
	    }
	    
	    cs_ipvec(N->pinv, B_new, x, size);
	    cs_lsolve(N->L, x);
	    cs_usolve(N->U, x);
	    cs_ipvec(S->q, x, B_new, size);
	    
	    fprintf(fp,"---------------------------------------------------------------------\n");
	      
	    fprintf(fp, "for %s:  %lf\t", dc_table->input_var, inc);
	    for(cur = head_plot; cur != NULL; cur = cur->next){
	      for(i = 0; i < size; i++){
		if ( strcmp(X[i].name,cur->node) == 0 ){
		  
		  for (j = 0; j < size; j++){
		      if (S->q[j] == i) break;
		  }
		  
		  fprintf(fp," %c(%3s): %g \t",cur->type, X[i].name, x[ j ]);
		}
	      }
	    }
	    fprintf(fp,"\n");
	  }
    }
    
    fprintf(fp,"---------------------------------------------------------------------\n");

    fclose(fp);
    
    cs_sfree(S);
    cs_nfree(N);
    cs_free(x);
    cs_free(B_new);
    
    return(0);
}


int cs_Cholesky (cs *C,double *B, int size, dc_analysis *dc_table, tPlot *head_plot, int thesiDC, tarrayX *X, dcIt dcInodes){								//14/11
    
    int i;
    css *S;
    csn *N;
    double *x;
    double inc;
    tPlot *cur;
    FILE *fp;
    double * B_new ; 
    
    fp=fopen("results.txt","w");
    
    S = cs_calloc (1, sizeof (css));
    N = cs_calloc (1, sizeof (csn)); 
    
    x = (double *)calloc(size,sizeof(double));
    if (x == NULL) {
      printf("Memory allocation failed. Exiting...\n");
      exit(1);
    }
    
    B_new = (double *)calloc(size,sizeof(double));
    if (B_new == NULL) {
      printf("Memory allocation failed. Exiting...\n");
      exit(1);
    }      
    
    for (i = 0; i < size; i++){
      B_new[i] =  B[i];
    }
    
    S = cs_schol(1, C);
    N = cs_chol(C, S);

    cs_spfree(C);	//free C
    
    fprintf(fp,"Oi luseis tis eksiswsis me Cholesky gia araious pinakes einai : \n");
    fflush(fp);
    
    cs_ipvec(S->pinv, B_new, x, size);
    cs_lsolve(N->L, x);
    cs_ltsolve(N->L, x);
    cs_pvec(S->pinv, x, B_new, size);
        
    fprintf(fp,"---------------------------------------------------------------------\n");
    
    for (i=0; i<size; i++) {
	  fprintf(fp, "x[%d] =  %g\n", i, x[ S->pinv[i] ]);
	  fflush(fp);
    }
    
    if (thesiDC != -1){
	  for (inc=dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	    B[thesiDC] = inc;
	    
	    for (i = 0; i < size; i++){
	      B_new[i] =  B[i];
	    }
	    
	    cs_ipvec(S->pinv, B_new, x, size);
	    cs_lsolve(N->L, x);
	    cs_ltsolve(N->L, x);
	    cs_pvec(S->pinv, x, B_new, size);
	    
	    fprintf(fp,"---------------------------------------------------------------------\n");
	      
	    fprintf(fp,"for %s:  %lf\t", dc_table->input_var, inc);
	    for (cur = head_plot; cur != NULL; cur = cur->next){
	      for (i = 0; i < size; i++){
		if ( strcmp(X[i].name,cur->node) == 0 ){
		  fprintf(fp," %c(%3s): %g \t",cur->type, X[i].name, x[ S->pinv[i] ] );
		}
	      }
	    }
	    fprintf(fp,"\n");
	  }
    }else if (dcInodes.node1 != -1){				
       for (inc = dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc += dc_table->increment){
	  
	  if (dcInodes.flag==0){
	    B[dcInodes.node1] = -inc;
	    if (dcInodes.node2 != -1)
	      B[dcInodes.node2] = inc;	
	  }else if (dcInodes.flag==1) {
	    B[dcInodes.node1] = inc;
	    if (dcInodes.node2 != -1)
	     B[dcInodes.node2] = -inc;
	  }
	  
	  for (i = 0; i < size; i++){
	    B_new[i] =  B[i];
	  }
			  
	  cs_ipvec(S->pinv, B_new, x, size);
	  cs_lsolve(N->L, x);
	  cs_ltsolve(N->L, x);
	  cs_pvec(S->pinv, x, B_new, size);
	  
	  fprintf(fp,"---------------------------------------------------------------------\n");
	    
	  fprintf(fp, "for %s:  %lf\t", dc_table->input_var, inc);
	  for(cur = head_plot; cur != NULL; cur = cur->next){
	    for(i = 0; i < size; i++){
	      if ( strcmp(X[i].name,cur->node) == 0 ){
		fprintf(fp," %c(%3s): %g \t",cur->type, X[i].name , x[ S->pinv[i]] );
	      }
	    }
	  }
	  fprintf(fp,"\n");
	}
    }
    
    fprintf(fp,"---------------------------------------------------------------------\n");

    fclose(fp);
    
    cs_sfree(S);
    cs_nfree(N);
    cs_free(x);
    free(B);
    
    return(0);
}
