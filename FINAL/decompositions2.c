/*		Circuit Simulation 2014 4 General Purpose Preconditioners
* 
* Onoma programmatisti:	Alfonsos Swthrhs 959
* 				Servou Katerina 852
* 				Stauropoulos Antonis 1141
* 				Vakerlis Christos 982
* 
* Hmeromhnia:		15/11/2014
*/


#define EPS 1e-18

int CG (double **A,double *B, tarrayX *X, double itol, int size, dc_analysis *dc_table, tPlot *head_plot, int thesiDC, dcIt dcInodes){
  
    double A_monodiastatos[size*size-1];
    int j;  
    int i;
    int thesh=0;
    int iter;
    double euclidean_B;
    double euclidean_r;
    double M[size];
    double rho;
    double p_inv_q;
    double rho1;
    double alpha;
    double beta;
    double inc;
    tPlot *cur;
	
    FILE *fp;
    fp=fopen("results.txt","w");

    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) {
	A_monodiastatos[thesh] = A[i][j];
	thesh++;
      }
    }
    
    gsl_set_error_handler_off();

    gsl_matrix_view a = gsl_matrix_view_array (A_monodiastatos, size, size);
    gsl_vector_view b = gsl_vector_view_array (B, size);
    gsl_vector *x = gsl_vector_calloc (size);				//x[]=0
    gsl_vector *y = gsl_vector_alloc (size);
    gsl_vector *r = gsl_vector_calloc (size);
    gsl_vector *temp_b = gsl_vector_calloc (size);
    gsl_vector *z = gsl_vector_calloc (size);
    gsl_vector *p = gsl_vector_calloc (size);
    gsl_vector *temp_z = gsl_vector_calloc (size);
    gsl_vector *q = gsl_vector_calloc (size);
    gsl_vector *q_temp = gsl_vector_calloc (size);
  
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, x, 0.0, y);		//y = A*x
    gsl_vector_memcpy(temp_b, &b.vector);				//temp_b = &b.vector
    gsl_vector_sub(temp_b,y);						//temp_b = b-Ax
    gsl_vector_memcpy(r, temp_b);					//r = b-A*x
    iter = 0;
    
    for (i=0; i<size; i++) {						//M=diag(A)
	  if(A[i][i]==0){
	      M[i] = 1;
	  }else{
	      M[i] = 1/A[i][i];
	  }
    }
    
    gsl_vector_view m = gsl_vector_view_array (M, size);		//m = M
    
    euclidean_r = gsl_blas_dnrm2(r);					//euclidean_r = ||r||
    euclidean_B = gsl_blas_dnrm2(&b.vector);				//     -||-
    if (euclidean_B == 0) euclidean_B = 1;
   
    gsl_vector_memcpy(z, r);					// z = r

    while (((euclidean_r/ euclidean_B) > itol) && (iter < size)) {
      iter++;
      gsl_vector_mul(z, &m.vector);				// z = r / M
      gsl_blas_ddot (r, z, &rho);				//rho=r*z
      
      if (iter == 1) 
	gsl_vector_memcpy(p, z);			//p=z
      else {
	beta = rho/rho1;
	gsl_vector_memcpy(temp_z, z);				// p= z+ beta*p
	gsl_blas_daxpy (beta, p, temp_z);
	gsl_vector_memcpy(p,temp_z);
      }
      
      rho1 = rho;
      gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, p, 0.0, q);	//q= Ap
      gsl_blas_ddot (p, q, &p_inv_q);		// alpha = rho / (pT * q)
      alpha = rho / p_inv_q; 
      gsl_blas_daxpy (alpha, p, x);		//x= x+ alpha * p	
      gsl_blas_daxpy (-alpha, q, r);		//r = r-alpha*q
      euclidean_r = gsl_blas_dnrm2(r);
      gsl_vector_memcpy(z, r);					// z = r
    }
    
    fprintf(fp,"Oi luseis tis eksiswsis me CG einai : \n");
    fprintf(fp,"---------------------------------------------------------------------\n");
    fprintf (fp,"x = \n");
    gsl_vector_fprintf (fp, x, "%g");
    printf("iter: %d",iter);
    
    if ( thesiDC != -1){
      for (inc=dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	  b.vector.data[thesiDC] = inc;
	  
	    gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, x, 0.0, y);	//y = A*x
	    gsl_vector_memcpy(temp_b, &b.vector);				//temp_b = &b.vector
	    gsl_vector_sub(temp_b,y);					//temp_b = b-Ax
	    gsl_vector_memcpy(r, temp_b);					//r = b-A*x
	    iter = 0;

	    euclidean_r = gsl_blas_dnrm2(r);				//euclidean_r = ||r||
	    euclidean_B = gsl_blas_dnrm2(&b.vector);			//     -||-
	    if (euclidean_B == 0) euclidean_B = 1;
	  
	    gsl_vector_memcpy(z, r);					// z = r

	    while (((euclidean_r/ euclidean_B) > itol)  && (iter < size)) {			//&& (iter < size)) {
	      iter++;
	      gsl_vector_mul(z, &m.vector);				// z = r / M
	      gsl_blas_ddot (r, z, &rho);					//rho=r*z
	      
	      if (iter == 1) 
		gsl_vector_memcpy(p, z);					//p=z
	      else {
		beta = rho/rho1;
		gsl_vector_memcpy(temp_z, z);				// p= z+ beta*p
		gsl_blas_daxpy (beta, p, temp_z);
		gsl_vector_memcpy(p,temp_z);
	      }
	      
	      rho1 = rho;
	      gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, p, 0.0, q);	//q= Ap
	      gsl_blas_ddot (p, q, &p_inv_q);				// alpha = rho / (pT * q)
	      alpha = rho / p_inv_q; 
	      gsl_blas_daxpy (alpha, p, x);				//x= x+ alpha * p	
	      gsl_blas_daxpy (-alpha, q, r);				//r = r-alpha*q
	      euclidean_r = gsl_blas_dnrm2(r);
	      gsl_vector_memcpy(z, r);					// z = r
	    }
	  
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
	  printf("iter: %d",iter);
      }
    }else if (dcInodes.node1 != -1){				//14/11
       for (inc=dc_table->start_value; inc< (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	  
	  if (dcInodes.flag==0){
	    b.vector.data[dcInodes.node1] = -inc;
	    if (dcInodes.node2 != -1)
	      b.vector.data[dcInodes.node2] = inc;	
	  }else if (dcInodes.flag==1) {
	    b.vector.data[dcInodes.node1] = inc;
	    if (dcInodes.node2 != -1)
	      b.vector.data[dcInodes.node2] = -inc;
	  }
	  
	      gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, x, 0.0, y);	//y = A*x
	      gsl_vector_memcpy(temp_b, &b.vector);				//temp_b = &b.vector
	      gsl_vector_sub(temp_b,y);					//temp_b = b-Ax
	      gsl_vector_memcpy(r, temp_b);					//r = b-A*x
	      iter = 0;

	      euclidean_r = gsl_blas_dnrm2(r);				//euclidean_r = ||r||
	      euclidean_B = gsl_blas_dnrm2(&b.vector);			//     -||-
	      if (euclidean_B == 0) euclidean_B = 1;
	    
	      gsl_vector_memcpy(z, r);					// z = r

	      while (((euclidean_r/ euclidean_B) > itol)  && (iter < size)) {			//&& (iter < size)) {
		iter++;
		gsl_vector_mul(z, &m.vector);				// z = r / M
		gsl_blas_ddot (r, z, &rho);					//rho=r*z
		
		if (iter == 1) 
		  gsl_vector_memcpy(p, z);					//p=z
		else {
		  beta = rho/rho1;
		  gsl_vector_memcpy(temp_z, z);				// p= z+ beta*p
		  gsl_blas_daxpy (beta, p, temp_z);
		  gsl_vector_memcpy(p,temp_z);
		}
		
		rho1 = rho;
		gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, p, 0.0, q);	//q= Ap
		gsl_blas_ddot (p, q, &p_inv_q);				// alpha = rho / (pT * q)
		alpha = rho / p_inv_q; 
		gsl_blas_daxpy (alpha, p, x);				//x= x+ alpha * p	
		gsl_blas_daxpy (-alpha, q, r);				//r = r-alpha*q
		euclidean_r = gsl_blas_dnrm2(r);
		gsl_vector_memcpy(z, r);					// z = r
	      }
	  
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
	  printf("iter: %d",iter);
	}
    }
    
    fprintf(fp,"---------------------------------------------------------------------\n");

    fclose(fp);
    
    gsl_vector_free (y);
    gsl_vector_free (r);
    gsl_vector_free (temp_b);
    gsl_vector_free (z);
    gsl_vector_free (p);
    gsl_vector_free (temp_z);
    gsl_vector_free (q);
    gsl_vector_free (q_temp);
    gsl_vector_free (x);
    
    return(0);
}

int Bi_CG (double **A,double *B, tarrayX *X, double itol, int size, dc_analysis *dc_table, tPlot *head_plot, int thesiDC, dcIt dcInodes){
  
    double A_monodiastatos[size*size-1];
    double A_t_monodiastatos[size*size-1];
    double M[size];
    int iter;
    int i;
    int j;
    int thesh = 0;
    double rho;
    double rho1;
    double alpha;
    double beta;
    double omega;
    double euclidean_r;
    double euclidean_B;
    double inc;
    tPlot *cur;
    
    FILE *fp;
    fp=fopen("results.txt","w");
    
    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) {
	A_monodiastatos[thesh] = A[i][j];
	A_t_monodiastatos[thesh]=A[j][i];
	thesh++;
      }
    }

    for (i=0; i<size; i++) {						//M=diag(A)
	  if (A[i][i]==0) {
	      M[i] = 1;
	  }else {
	      M[i] = 1/A[i][i];
	  }
    }
    
    gsl_set_error_handler_off();
    
    gsl_matrix_view a = gsl_matrix_view_array (A_monodiastatos, size, size);
    gsl_matrix_view a_T = gsl_matrix_view_array (A_t_monodiastatos, size, size);
    gsl_vector_view m = gsl_vector_view_array (M, size);
    gsl_vector_view b = gsl_vector_view_array (B, size);
    gsl_vector *x = gsl_vector_calloc (size);
    gsl_vector *r = gsl_vector_calloc (size);
    gsl_vector *r_new = gsl_vector_calloc (size);
    gsl_vector *z = gsl_vector_calloc (size);
    gsl_vector *z_temp = gsl_vector_calloc (size);
    gsl_vector *temp_b = gsl_vector_calloc (size);
    gsl_vector *z_new = gsl_vector_calloc (size);
    gsl_vector *p = gsl_vector_calloc (size);
    gsl_vector *p_new = gsl_vector_calloc (size);
    gsl_vector *q = gsl_vector_calloc (size);
    gsl_vector *q_new = gsl_vector_calloc (size);
    gsl_vector *y = gsl_vector_alloc (size);
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, x, 0.0, y);		//y = A*x
    gsl_vector_memcpy(temp_b, &b.vector);				//temp_b = &b.vector
  
    gsl_vector_sub(temp_b,y);						//temp_b = b-Ax
    gsl_vector_memcpy(r, temp_b);	
    gsl_vector_memcpy(r_new, temp_b);
    iter = 0;
    
    euclidean_r = gsl_blas_dnrm2(r);					//euclidean_r = ||r||
    euclidean_B = gsl_blas_dnrm2(&b.vector);				//     -||-
    
    if (euclidean_B == 0) euclidean_B = 1;
  
    
    gsl_vector_memcpy(z, r);					// z = r
    gsl_vector_memcpy(z_new, r_new);				//z_new = r_new
  
    
    while (((euclidean_r/ euclidean_B) > itol)   ){//   && (iter < size)) {
      iter++;
      gsl_vector_mul(z, &m.vector);				// z = r / M
      gsl_vector_mul(z_new, &m.vector);				// z = r / M
      gsl_blas_ddot (z, r_new, &rho);				//rho=r*z
    
      
      if (fabs(rho)<EPS) exit(1);
	
      if (iter==1) {
	gsl_vector_memcpy(p,z);
	gsl_vector_memcpy(p_new,z_new);
      }else {
	beta=rho/rho1;
	
	gsl_vector_memcpy(z_temp,z);
	gsl_blas_daxpy (beta, p, z_temp);		//z_temp= z_temp+ beta * p
	gsl_vector_memcpy(p,z_temp);
	
	gsl_vector_memcpy(z_temp,z_new);
	gsl_blas_daxpy (beta, p_new, z_temp);		//z_temp= z_temp+ beta * p
	gsl_vector_memcpy(p_new,z_temp);
      }
      
      rho1=rho;
      gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, p, 0.0, q);	//q= Ap
      gsl_blas_dgemv(CblasNoTrans, 1.0, &a_T.matrix, p_new, 0.0, q_new);	//q_new= A_mono*p_new
      gsl_blas_ddot (p_new, q, &omega);
      
      if (fabs(omega)<EPS){ printf("omega %.20lf\n",omega);   exit(1);}
      
      alpha = rho/omega;
      gsl_blas_daxpy (alpha, p, x);		//x= x+ alpha * p
      gsl_blas_daxpy (-alpha, q, r);		//r= r- alpha * q
      gsl_blas_daxpy (-alpha, q_new, r_new);		//r_new= r_new- alpha * q_new   
      
      euclidean_r = gsl_blas_dnrm2(r);
      
      gsl_vector_memcpy(z, r);					// z = r
      gsl_vector_memcpy(z_new, r_new);				//z_new = r_new
    }
    
    
    fprintf(fp,"Oi luseis tis eksiswsis me BI_CG einai : \n");
    fprintf(fp,"---------------------------------------------------------------------\n");
    fprintf (fp,"x = \n");
    gsl_vector_fprintf (fp, x, "%g");
    printf("iter: %d ",iter);
    
    if (thesiDC != -1) {

      for (inc=dc_table->start_value; inc < (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	b.vector.data[thesiDC] = inc;

	    gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, x, 0.0, y);		//y = A*x
	    gsl_vector_memcpy(temp_b, &b.vector);				//temp_b = &b.vector
	    gsl_vector_sub(temp_b,y);						//temp_b = b-Ax
	    gsl_vector_memcpy(r, temp_b);	
	    gsl_vector_memcpy(r_new, r);
	    iter = 0;
	    euclidean_r = gsl_blas_dnrm2(r);					//euclidean_r = ||r||
	    euclidean_B = gsl_blas_dnrm2(&b.vector);				//     -||-
	    if (euclidean_B == 0) euclidean_B = 1;
	    gsl_vector_memcpy(z, r);					// z = r
	    gsl_vector_memcpy(z_new, r_new);				//z_new = r_new
	    
	    while (((euclidean_r/ euclidean_B) > itol)) {//   && (iter < size)) {
	      iter++;
	      gsl_vector_mul(z, &m.vector);				// z = r / M
	      gsl_vector_mul(z_new, &m.vector);				// z_new = r / M
	      gsl_blas_ddot (z, r_new, &rho);				//rho=r*z

	      if (fabs(rho)<EPS) exit(1);
		
	      if (iter==1) {
		gsl_vector_memcpy(p,z);
		gsl_vector_memcpy(p_new,z_new);
	      }else {
		beta=rho/rho1;
		gsl_vector_memcpy(z_temp,z);
		gsl_blas_daxpy (beta, p, z_temp);		//z_temp= z_temp+ beta * p
		gsl_vector_memcpy(p,z_temp);
		gsl_vector_memcpy(z_temp,z_new);
		gsl_blas_daxpy (beta, p_new, z_temp);		//z_temp= z_temp+ beta * p
		gsl_vector_memcpy(p_new,z_temp);
	      }
	      
	      rho1=rho;
	      gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, p, 0.0, q);	//q= Ap
	      gsl_blas_dgemv(CblasNoTrans, 1.0, &a_T.matrix, p_new, 0.0, q_new);	//q_new= A_mono*p_new
	      gsl_blas_ddot (p_new, q, &omega);
	     if (fabs(omega)<EPS){printf("omega %.20lf\n",omega);  exit(1);}
	      alpha = rho/omega;
	      gsl_blas_daxpy (alpha, p, x);		//x= x+ alpha * p
	      gsl_blas_daxpy (-alpha, q, r);		//r= r- alpha * q
	      gsl_blas_daxpy (-alpha, q_new, r_new);		//r_new= r_new- alpha * q_new   
	      euclidean_r = gsl_blas_dnrm2(r);
	      gsl_vector_memcpy(z, r);					// z = r
	      gsl_vector_memcpy(z_new, r_new);				//z_new = r_new
	    }
      
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
	printf("iter: %d ",iter);
      }
    }else if (dcInodes.node1 != -1){				//14/11

       for (inc=dc_table->start_value; inc< (dc_table->end_value + dc_table->increment); inc+=dc_table->increment){
	  
	  if (dcInodes.flag==0){
	    b.vector.data[dcInodes.node1] = -inc;
	    if (dcInodes.node2 != -1)
	      b.vector.data[dcInodes.node2] = inc;	
	  }else if (dcInodes.flag==1) {
	    b.vector.data[dcInodes.node1] = inc;
	    if (dcInodes.node2 != -1)
	      b.vector.data[dcInodes.node2] = -inc;
	  }
	  	  
	    gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, x, 0.0, y);		//y = A*x
	    gsl_vector_memcpy(temp_b, &b.vector);				//temp_b = &b.vector
	    gsl_vector_sub(temp_b,y);						//temp_b = b-Ax
	    gsl_vector_memcpy(r, temp_b);	
	    gsl_vector_memcpy(r_new, temp_b);
	    iter = 0;
	    euclidean_r = gsl_blas_dnrm2(r);					//euclidean_r = ||r||
	    euclidean_B = gsl_blas_dnrm2(&b.vector);				//     -||-
	    if (euclidean_B == 0) euclidean_B = 1;
	    gsl_vector_memcpy(z, r);					// z = r
	    gsl_vector_memcpy(z_new, r_new);				//z_new = r_new
	    
	    while (((euclidean_r/ euclidean_B) > itol) ) {//   && (iter < size)) {
	      iter++;
	      gsl_vector_mul(z, &m.vector);				// z = r / M
	      gsl_vector_mul(z_new, &m.vector);				// z = r / M
	      gsl_blas_ddot (z, r_new, &rho);				//rho=r*z
	    
	      if (fabs(rho)<EPS) exit(1);
		
	      if (iter==1) {
		gsl_vector_memcpy(p,z);
		gsl_vector_memcpy(p_new,z_new);
	      }else {
		beta=rho/rho1;
		gsl_vector_memcpy(z_temp,z);
		gsl_blas_daxpy (beta, p, z_temp);		//z_temp= z_temp+ beta * p
		gsl_vector_memcpy(p,z_temp);
		gsl_vector_memcpy(z_temp,z_new);
		gsl_blas_daxpy (beta, p_new, z_temp);		//z_temp= z_temp+ beta * p
		gsl_vector_memcpy(p_new,z_temp);
	      }
	      
	      rho1=rho;
	      gsl_blas_dgemv(CblasNoTrans, 1.0, &a.matrix, p, 0.0, q);	//q= Ap
	      gsl_blas_dgemv(CblasNoTrans, 1.0, &a_T.matrix, p_new, 0.0, q_new);	//q_new= A_mono*p_new
	      gsl_blas_ddot (p_new, q, &omega);
	      
	     if (fabs(omega)<EPS){ printf("omega %.20lf\n",omega);   exit(1);}
	      
	      alpha = rho/omega;
	      gsl_blas_daxpy (alpha, p, x);		//x= x+ alpha * p
	      gsl_blas_daxpy (-alpha, q, r);		//r= r- alpha * q
	      gsl_blas_daxpy (-alpha, q_new, r_new);		//r_new= r_new- alpha * q_new   
	      euclidean_r = gsl_blas_dnrm2(r);
	      gsl_vector_memcpy(z, r);					// z = r
	      gsl_vector_memcpy(z_new, r_new);				//z_new = r_new
	    }
	  
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
	  printf("iter: %d ",iter);
	}
    }
    
    fprintf(fp,"---------------------------------------------------------------------\n");
 
      
      
     gsl_vector_free (z_new);
     gsl_vector_free (r);
     gsl_vector_free (r_new);
     gsl_vector_free (temp_b);
     gsl_vector_free (z);
     gsl_vector_free (p);
     gsl_vector_free (p_new);
     gsl_vector_free (z_temp);
     gsl_vector_free (q);
     gsl_vector_free (q_new);
     gsl_vector_free (x);
     gsl_vector_free (y);
      
  return(0); 
}

    