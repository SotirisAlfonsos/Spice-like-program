/*		Circuit Simulation 2014 3 Solution of MNA System
 * 
 * Onoma programmatisti:	Alfonsos Swthrhs 959
 * 				Servou Katerina 852
 * 				Stauropoulos Antonis 1141
 * 				Vakerlis Christos 982
 * 
 * Hmeromhnia:		01/11/2014
 */

int decoLU (double **A,double *B, int size, dc_analysis *dc_table, tPlot *head_plot, int thesiDC, tarrayX *X, dcIt dcInodes){								//14/11
    double patates[size*size-1];
    int j;  
    int i;
    int thesh=0;
    int s;
    double inc;
    tPlot *cur;
    
    FILE *fp;
    fp=fopen("results.txt","w");
    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) {
	patates[thesh] = A[i][j];
	thesh++;
      }
    }

    gsl_matrix_view m = gsl_matrix_view_array (patates, size, size); 	//m.matrix= A
    
    gsl_vector_view b = gsl_vector_view_array (B, size);		//b.vector= B

    gsl_vector *x = gsl_vector_alloc (size);				//alloc(x)

    gsl_permutation * p = gsl_permutation_alloc (size);			//alloc(p)
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);				
    
    fprintf(fp,"Oi luseis tis eksiswsis me LU einai : \n");
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    fprintf(fp,"---------------------------------------------------------------------\n");
    fprintf (fp,"x = \n");
    gsl_vector_fprintf (fp, x, "%g");
    
    if (thesiDC != -1){
      for (inc=dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	b.vector.data[thesiDC] = inc;
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	fprintf(fp,"---------------------------------------------------------------------\n");
	  
	fprintf(fp,"for %s:  %lf\t",dc_table->input_var,inc);
	for(cur=head_plot; cur!=NULL; cur=cur->next){
	  for(i=0; i<size; i++){
	    if ( strcmp(X[i].name,cur->node) == 0 ){
	      fprintf(fp," %c(%3s): %g \t",cur->type,X[i].name,gsl_vector_get (x, i));
	    }
	  }
	}
	fprintf(fp,"\n");
      }
    }else if (dcInodes.node1 != -1){				//14/11
       for (inc=dc_table->start_value; inc< (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	  
	  if (dcInodes.flag==0){
	    b.vector.data[dcInodes.node1] = -inc;
	    if (dcInodes.node2 != -1)
	      b.vector.data[dcInodes.node2] = inc;	
	  }else if(dcInodes.flag==1) {
	    b.vector.data[dcInodes.node1] = inc;
	    if (dcInodes.node2 != -1)
	      b.vector.data[dcInodes.node2] = -inc;
	  }
	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	  fprintf(fp,"---------------------------------------------------------------------\n");
	    
	  fprintf(fp,"for %s:  %lf\t",dc_table->input_var,inc);
	  for(cur=head_plot; cur!=NULL; cur=cur->next){
	    for(i=0; i<size; i++){
	      if ( strcmp(X[i].name,cur->node) == 0 ){
		fprintf(fp," %c(%3s): %g \t",cur->type,X[i].name,gsl_vector_get (x, i));
	      }
	    }
	  }
	  fprintf(fp,"\n");
	}
    }
    
    fprintf(fp,"---------------------------------------------------------------------\n");

    fclose(fp);
    gsl_permutation_free (p);
    gsl_vector_free (x);
    return 0;
}


int decoCHOLESKY (double **A,double *B, int size, dc_analysis *dc_table, tPlot *head_plot, int thesiDC, tarrayX *X, dcIt dcInodes){								//14/11
    double patates[size*size-1];
    int j;  
    int i;
    int thesh=0;
    double inc;
    FILE *fp;
    tPlot *cur;
    fp=fopen("results.txt","w");
    
    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) {
	patates[thesh] = A[i][j];
	thesh++;
      }
    }
    gsl_set_error_handler_off();

    gsl_matrix_view m = gsl_matrix_view_array (patates, size, size);
    
    gsl_vector_view b = gsl_vector_view_array (B, size);

    gsl_vector *x = gsl_vector_alloc (size);

    gsl_linalg_cholesky_decomp (&m.matrix);

    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, x);
    fprintf(fp,"Oi luseis tis eksiswsis me CHOLESKY einai : \n");
    fprintf(fp,"---------------------------------------------------------------------\n");
    fprintf (fp,"x = \n");
    gsl_vector_fprintf (fp, x, "%g");
    
    if (thesiDC != -1){
      for (inc=dc_table->start_value; inc< (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	b.vector.data[thesiDC] = inc;
	gsl_linalg_cholesky_solve (&m.matrix, &b.vector, x);
	fprintf(fp,"---------------------------------------------------------------------\n");
	
	fprintf(fp,"for %s:  %lf\t",dc_table->input_var,inc);
	for(cur=head_plot; cur!=NULL; cur=cur->next){
	  for(i=0; i<size; i++){
	    if ( strcmp(X[i].name,cur->node) == 0 ){
	      fprintf(fp," %c(%3s): %g \t",cur->type,X[i].name,gsl_vector_get (x, i));
	    }
	  }
	}
	fprintf(fp,"\n");
      }
    }else if (dcInodes.node1 != -1){				//14/11
       for (inc=dc_table->start_value; inc< (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	
	if (dcInodes.flag==0){
	  b.vector.data[dcInodes.node1] = -inc;
	  if (dcInodes.node2 != -1)
	    b.vector.data[dcInodes.node2] = inc;	
	}else if(dcInodes.flag==1) {
	  b.vector.data[dcInodes.node1] = inc;
	  if (dcInodes.node2 != -1)
	    b.vector.data[dcInodes.node2] = -inc;
	}
	gsl_linalg_cholesky_solve (&m.matrix, &b.vector, x);
	fprintf(fp,"---------------------------------------------------------------------\n");
	  
	fprintf(fp,"for %s:  %lf\t",dc_table->input_var,inc);
	for(cur=head_plot; cur!=NULL; cur=cur->next){
	  for(i=0; i<size; i++){
	    if ( strcmp(X[i].name,cur->node) == 0 ){
	      fprintf(fp," %c(%3s): %g \t",cur->type,X[i].name,gsl_vector_get (x, i));
	    }
	  }
	}
	fprintf(fp,"\n");
      }
    }
    
    
    fprintf(fp,"---------------------------------------------------------------------\n");
    
    fclose(fp);
    gsl_vector_free (x);
    
    return 0;
}


  