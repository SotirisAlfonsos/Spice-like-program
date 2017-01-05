/*		Circuit Simulation 2014 6 Transient Analysis
 * 
 * Onoma programmatisti:	Alfonsos Swthrhs 959
 * 				Servou Katerina 852
 * 				Stauropoulos Antonis 1141
 * 				Vakerlis Christos 982
 * 
 * Hmeromhnia:		18/12/2014
 */


#define PI 3.141593

/************************************************************************************************************************
 * This function computes the [transient_spec] EXP(i1 i2 td1 tc1 td2 tc2). The parameters time and final_time are required	*
 * for finding which expression to be used. Returns the computed value. If the parameter time is out of time boundaries	*
 * the program terminates.										*
 ************************************************************************************************************************/
double EXP_func(double i1, double i2, double td1, double tc1, double td2, double tc2, double tk, double fin) {

	double exp_power1, exp_power2;

	if ( (tk >= 0.0) && (tk <= td1) ) return i1;
	
	else if ( (tk > td1) && (tk <= td2) ) {
		exp_power1 = (-1.0)*(tk - td1)/tc1;
		return ( i1 + (i2-i1)*(1.0 - exp(exp_power1)) );
	}
	else if ( (tk > td2) && (tk <= fin) ) {
		exp_power1 = (-1.0)*(tk - td1)/tc1;
		exp_power2 = (-1.0)*(tk - td2)/tc2;
		return ( i1 + (i2-i1)*(exp(exp_power2) - exp(exp_power1)) );
	}    
	else {
		printf("Error in EXP calculation\n");
		exit(1);
	}

}


/************************************************************************************************************************
 * This function computes the [transient_spec] SIN(i1 ia fr td df ph). The parameters time and final_time are required	*
 * for finding which expression to be used. Returns the computed value. If the parameter time is out of time boundaries	*
 * the program terminates.										*
 ************************************************************************************************************************/
double SIN_func(double i1, double ia, double fr, double td, double df, double ph, double tk, double fin) {

	double exp_power;

	if ( (tk >= 0.0) && (tk <= td) ) return( i1+ia*sin(PI*ph/180) );
	else if ( (tk > td) && (tk <= fin) ) { 
		exp_power = (-1.0) * (tk-td) * df;  
		return( i1 + ia*sin(2*PI*fr*(tk-td) + PI*ph/180)*exp(exp_power) );
	}
	else {
 		printf("Error in SIN calculation\n");
		exit(1);
	}
}


/************************************************************************************************************************
 * This function computes the [transient_spec] PULSE(i1 i2 td tr tf pw per). The parameters time and final_time are	*
 * required for finding which expression to be used. It returns the computed value. If the parameter time is out of	*
 * time boundaries, the program terminates.								*
 ************************************************************************************************************************/
double PULSE_func(double i1, double i2, double td, double tr, double tf, double pw, double per, double tk, double fin) {

	int k;
	
	k = (int)(tk/per);
	
	//if ( (tk >= 0) && (tk <= (td)))//+k*per)) ) 
	if ( (tk >= 0) && (tk <= (td+k*per)) ) 
	  return(i1);
	else if ( (tk >= (td + k*per)) && (tk <= (td + tr + k*per)) ) 
	  return( ((i2-i1)/(td + tr + k*per - td + k*per)) * (tk - td + k*per) + i1 );
	else if ( (tk >= (td + tr + k*per)) && (tk <= (td + tr + pw + k*per)) ) 
	  return i2;
	else if ( (tk >= (td + tr + pw + k*per)) && (tk <= (td + tr + pw + tf + k*per)) ) 
	  return( (i1-i2)/( (td + tr + pw + tf + k*per) - (td + tr + pw + k*per)) * (tk - (td + tr + pw + k*per)) + i2 );
	else if ( (tk >= (td + tr + pw + tf + k*per)) && (tk <= (td + per + k*per)) ) 
	  return i1;
	else {
	  printf("Error in PULSE calculation\n");
	  exit(1);
	}		

}


/************************************************************************************************************************
 * This function computes the [transient_spec] PWL(t1 v1 t2 v2 ... tn vn). The parameters time and final_time are	*
 * required for finding which expression to be used. It returns the computed value. If the parameter time is out of	*
 * time boundaries, the program terminates.								*
 ************************************************************************************************************************/
double PWL_func(int n, double *t, double *i, double tk, double fin) {
	double temp=0;
	int k;
	int flag = 0;
	
	for(k = 0; k < n; k++) {
	    if ( (t[k] <= tk) && (tk <= t[k+1]) )   {
		temp = temp + ( (i[k+1] - i[k])/(t[k+1] - t[k]) * (tk - t[k]) + i[k] );
		flag = 1;
	    }
	    
	}
	if ( tk > t[n-1] ){
		temp = i[n-1];
	  
	}else if (!flag){
	  printf("Error in PWL calculation\n");
	  exit(1); 
	}
	return temp;
			
}

 /************************************************************************************************************************
 * Calculates the left hand side of the equation for the Backward Euler (BE) method					 *
 ************************************************************************************************************************/
double **BE_A(double h, double **G, double **C, int size) {
	double **A;
	int i, j;
	double temp_elem;
	A = (double **)calloc(size, sizeof(double *));
	
	if (A != NULL) {
	  for (i = 0; i < size; i++) {
	    A[i] = (double *)calloc(size, sizeof(double));
	    if (A[i] == NULL) {
	      printf("Memory Allocation Failed...\n");
	      exit(1);
	    }
	  }
	}
	else {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	
        for (i = 0; i < size; i++) {
          for (j = 0; j < size; j++) {  
	    temp_elem = C[i][j]*1.0/h;
            A[i][j] = G[i][j] + temp_elem;  
          }
        }
        return A;
  
}

 /*************************************************************************************
 *Uses the Backward Euler method to calculate the new x values. 		      *
 *Calculates the right hand side of the equation, changes the e() array according to  *
 *the transient values of the changing elements and uses LU decomposition to find     *
 *the new x.									      *
 *************************************************************************************/ 
void BE_b(double *B, double h, double fin, int size, tElements *head_elem, tSources *sources, double **G, double **C, tPlot *head_plot, tarrayX *X) {
	int i=0;
	double x_prev[size];
	tElements *curr;
	tnumbers *curr2;
	tSources *curr3;
	tPlot *curPlot;
	double patates[size*size-1];
        int k,j,s,thesiPlot=0;
	int thesh=0;
	int sum_numbers = 0;
	double *numbers = NULL;
	double temp_array[size];
	double trans_value;
	double *t_pwl;
	double *i_pwl;
	double *time_vector;	
	int sumPlots=0;	
	char *PlotFileName;
	gnuplot_ctrl **h_plot;
	double **A=BE_A(h, G, C, size);
	double **liseis;
	
	double *tempB = (double*)calloc(size,sizeof(double));

	for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	    sumPlots++;
	}
	
	h_plot = (gnuplot_ctrl **)malloc(sumPlots*sizeof(gnuplot_ctrl *));
	if (h_plot == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	
	liseis = (double **)malloc(sumPlots*sizeof(double *));
	if (h_plot == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	else{
	    for(i=0; i<sumPlots; i++){
		liseis[i] = (double *)calloc( (int)(fin/h) , sizeof(double ) );
		if (liseis[i] == NULL) {
		    printf("Memory Allocation Failed...\n");
		    exit(1);
		}	  
	    }
	}
	
	time_vector = (double *)calloc( (int)(fin/h), sizeof(double) );
	if (time_vector == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}		
	
	for (i=0; i<(int)fin/h; i++) time_vector[i] = h*i;
	
	for (i=0; i<size; i++) {
	  for (j=0; j<size; j++) {
	    patates[thesh] = G[i][j];
	    thesh++;
	  }
	}
	thesh = 0;
	
	gsl_matrix_view m = gsl_matrix_view_array (patates, size, size); 	//m.matrix= A
	gsl_vector_view b = gsl_vector_view_array (B, size);			//b.vector= B
	gsl_vector *x = gsl_vector_alloc (size);				//alloc(x)
	gsl_permutation *p = gsl_permutation_alloc (size);			//alloc(p)
	
	gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	
	for(i=0; i<size; i++){
	  x_prev[i]=gsl_vector_get (x, i);
	}

	for(k=1; k*h<=fin; k++) {
	  for (i = 0; i < size; i++) {
	    temp_array[i] = 0;
	      for (j = 0; j < size; j++) {  
		  temp_array[i] = temp_array[i] + C[i][j] * 1.0/h * x_prev[j];
		  //if (i==1) 
		  //printf("temparray[%d}: %g C[%d][%d]: %g x_prev[%d]: %g\n",i,temp_array[i],i,j,C[i][j],j,x_prev[j]);
		  //A[i][j] = G[i][j] + temp_elem;  
	      }
	  }
	  // calculate e(t) which changes at every step of transient analysis  
	  for (curr = head_elem; curr != NULL; curr = curr->next) {
	      if (curr->trans!=NULL) {
		  sum_numbers=0;
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    sum_numbers++;
		  }
		  numbers = (double *)malloc(sum_numbers*sizeof(double));
		  if (numbers == NULL) {
		      printf("Memory Allocation Failed...\n");
		      exit(1);
		  }

		  sum_numbers = 0;
		  
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    numbers[sum_numbers] = curr2->number;
		    sum_numbers++;
		  }
		    
		  if (strcmp(curr->trans->method,"EXP") == 0 || strcmp(curr->trans->method,"exp") == 0) {
		    
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = EXP_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    else{
			trans_value = EXP_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    
		  }
		  else if (strcmp(curr->trans->method,"SIN") == 0 || strcmp(curr->trans->method,"sin") == 0) {
		    
		     if( strcmp(curr->posNode,"0")==0){
			trans_value = SIN_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }
		     else{
			trans_value = SIN_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }		      
		  }
		  else if (strcmp(curr->trans->method,"PULSE") == 0 || strcmp(curr->trans->method,"pulse") == 0) {
		    
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = PULSE_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }
		    else{ 
			trans_value = PULSE_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }  
		  }
		  else if (strcmp(curr->trans->method,"PWL") == 0 || strcmp(curr->trans->method,"pwl") == 0) {
		      t_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (t_pwl == NULL) {
			  printf("Memory Allocation Failed...\n");
			  exit(1);
		      }
		      
		      i_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (i_pwl == NULL) {
			  printf("Memory Allocation Failed...\n");
			  exit(1);
		      }
		    
		    for (j = 0; j<sum_numbers; j = j + 2) {   		      
		      t_pwl[j>>1] = numbers[j];
		      if( strcmp(curr->posNode,"0")==0){
			i_pwl[j>>1] = -numbers[j+1];
		      }
		      else{
			i_pwl[j>>1] = numbers[j+1];
		      }		      
		    }
		    
		    trans_value = PWL_func((int)sum_numbers/2, t_pwl, i_pwl, k*h, fin);
		    
		  }
		  else {
		    printf("Error in Transient Section\n");
		    exit(1);
		  }
		  if (curr->type == 'v' || curr->type == 'V'){
		      for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      tempB[curr3->thesiB] =  trans_value;//B[curr3->thesiB] + trans_value;
		  }else if ( (curr->type == 'i' || curr->type == 'I') ) {
		     for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      
		      /****************************************************************/
		      tempB[curr3->thesiB] = - trans_value; // B[curr3->thesiB] - trans_value; 
		      /***************************************************************/
		      if (curr3->thesiB_pairNode!=-1){
			  tempB[curr3->thesiB_pairNode] =  trans_value;		//B[curr3->thesiB_pairNode] + trans_value;	
			  
		      }
		  }
	
		 // free(numbers);
		  //free(t_pwl);
		 // free(i_pwl);
	      }
	      //free(numbers);
	  }

	  for (i=0; i<size; i++) {
	    temp_array[i] = temp_array[i] + tempB[i];
	    for (j=0; j<size; j++) {
	      patates[thesh] = A[i][j];
	      thesh++;
	    }
	  }
	  thesh = 0;
	  
	  m = gsl_matrix_view_array (patates, size, size); 	//m.matrix= A
	  b = gsl_vector_view_array (temp_array, size);		//b.vector= B

	  gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	  
	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	  for(i=0; i<size; i++){
	    x_prev[i]=gsl_vector_get (x, i);
	  }

	  thesiPlot=0;
	  for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	      for(i=0; i<size; i++){
		  
		  if ( strcmp(X[i].name,curPlot->node) == 0 ){
		    liseis[ thesiPlot ][k-1] = x_prev[i];
		  }
	      }
	      thesiPlot++;
	  }
	  
	}
	printf("BE\n");
	thesiPlot=0;
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot = curPlot->next){
	  
	      /***dimiourgia onomatos arxiou eksodou***/
	    sumPlots = strlen(curPlot->node) + 25;
	    PlotFileName = (char  *)malloc(sumPlots*sizeof( char ));
	    strcpy(PlotFileName,"set output \"BE_ \0");
	    PlotFileName[ strlen(PlotFileName)- 1] = curPlot->type;
	    strcat(PlotFileName, "(\0");
	    strcat(PlotFileName, curPlot->node);
	    strcat(PlotFileName, ").png\"\0");
	    
	      /***anoigma .png output file kai ektiposi tou apotelesmatos se auto***/
	    h_plot[thesiPlot] = gnuplot_init();
	    gnuplot_setstyle(h_plot[thesiPlot], "lines");
	    gnuplot_cmd(h_plot[thesiPlot], "set terminal png");
	    gnuplot_cmd(h_plot[thesiPlot], PlotFileName);
	    gnuplot_set_xlabel(h_plot[thesiPlot], "Time");
	    
	    if(curPlot->type=='I' || curPlot->type=='i') 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "I");
	    else 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "V");
	    
	    PlotFileName[ strlen(PlotFileName) - 5 ] = '\0';
	    gnuplot_plot_xy(h_plot[thesiPlot], time_vector, liseis[ thesiPlot ], (int)fin/h, &PlotFileName[15]);

	    gnuplot_close(h_plot[thesiPlot]);
	    thesiPlot++;
	}
	for (i=0; i< size; i++) {
	  free(A[i]);
	  
	}
	//free(tempB);
	//free(time_vector);
}

 /************************************************************************************************************************
 * Calculates the left hand side of the equation for the Trapezoidal (TR) method					 *
 ************************************************************************************************************************/
double **TR_A(double h, double **G, double **C, int size) {
	double **A;
	int i, j;
	double temp_elem;
	A = (double **)calloc(size, sizeof(double *));
	
	if (A != NULL) {
	  for (i = 0; i < size; i++) {
	    A[i] = (double *)calloc(size, sizeof(double));
	    if (A[i] == NULL) {
	      printf("Memory Allocation Failed...\n");
	      exit(1);
	    }
	  }
	}
	else {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	
        for (i = 0; i < size; i++) {
          for (j = 0; j < size; j++) {  
	    temp_elem = C[i][j]*2.0/h;
            A[i][j] = G[i][j] + temp_elem;  
          }
        }
        
        return A;
}


void TR_b(double *B, double h, double fin, int size, tElements *head_elem, tSources *sources, double **G, double **C, tPlot *head_plot, tarrayX *X) {
	int i=0;
	double x_prev[size];
	tElements *curr;
	tnumbers *curr2;
	tSources *curr3;
	tPlot *curPlot;
	double patates[size*size-1];
        int k,j,s,thesiPlot=0;
	int thesh=0;
	int sum_numbers = 0;
	double *numbers = NULL;
	double temp_array[size];
	double trans_value;
	double *t_pwl;
	
	double *i_pwl;
	double *time_vector;	
	int sumPlots=0;	
	char *PlotFileName;
	gnuplot_ctrl **h_plot;
	double **A=TR_A(h, G, C, size);
	double **liseis;
	
	double *tempB = (double*)calloc(size,sizeof(double));
	double *e_prev = (double*)malloc(size*sizeof(double)); //			<--------------
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	    sumPlots++;
	}
	
	h_plot = (gnuplot_ctrl **)malloc(sumPlots*sizeof(gnuplot_ctrl *));
	if (h_plot == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	
	liseis = (double **)malloc(sumPlots*sizeof(double *));
	if (liseis == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	else{
	    for(i=0; i<sumPlots; i++){
		liseis[i] = (double *)calloc( (int)(fin/h) , sizeof(double ) );
		if (liseis[i] == NULL) {
		    printf("Memory Allocation Failed...\n");
		    exit(1);
		}	  
	    }
	}
	
	time_vector = (double *)calloc( (int)(fin/h), sizeof(double) );
	if (time_vector == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}		
	
	for (i=0; i<(int)fin/h; i++) time_vector[i] = h*i;
	
	for (i=0; i<size; i++) {
	  e_prev[i] = B[i];		//    					<-----------------
	  for (j=0; j<size; j++) {
	    patates[thesh] = G[i][j];
	    thesh++;
	  }
	}
	thesh = 0;
	
	gsl_matrix_view m = gsl_matrix_view_array (patates, size, size); 	//m.matrix= A
	gsl_vector_view b = gsl_vector_view_array (B, size);			//b.vector= B
	gsl_vector *x = gsl_vector_alloc (size);				//alloc(x)
	gsl_permutation *p = gsl_permutation_alloc (size);			//alloc(p)
	
	gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	for(i=0; i<size; i++) {
	  x_prev[i]=gsl_vector_get (x, i);
	  //printf("x[%d]: %g\n",i,x_prev[i]);
	}

	for(k=1; k*h<=fin; k++) {
	  for (i = 0; i < size; i++) {
	    temp_array[i] = 0;
	      for (j = 0; j < size; j++) {  
		  temp_array[i] = temp_array[i] - (G[i][j]-C[i][j] * 2.0/h) * x_prev[j];
		  //A[i][j] = G[i][j] + temp_elem;  
	      }
	  }
	  // calculate e(t) which changes at every step of transient analysis  
	  for (curr = head_elem; curr != NULL; curr = curr->next) {
	      if (curr->trans!=NULL) {
		  sum_numbers=0;
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    sum_numbers++;
		  }
		  
		  numbers = (double *)malloc(sum_numbers*sizeof(double));
		  if (numbers == NULL) {
		      printf("Memory Allocation Failed...\n");
		      exit(1);
		  }

		  sum_numbers = 0;
		  
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    numbers[sum_numbers] = curr2->number;
		    sum_numbers++;
		  }

		   if (strcmp(curr->trans->method,"EXP") == 0 || strcmp(curr->trans->method,"exp") == 0) {
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = EXP_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    else{
			trans_value = EXP_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    
		  }
		   else if (strcmp(curr->trans->method,"SIN") == 0 || strcmp(curr->trans->method,"sin") == 0) {
		     if( strcmp(curr->posNode,"0")==0){
			trans_value = SIN_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }
		     else{
			trans_value = SIN_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }      
		  }
		else if (strcmp(curr->trans->method,"PULSE") == 0 || strcmp(curr->trans->method,"pulse") == 0) {	    
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = PULSE_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }
		    else{ 
			trans_value = PULSE_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }	    
		  }
		   else if (strcmp(curr->trans->method,"PWL") == 0 || strcmp(curr->trans->method,"pwl") == 0) {
		      t_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (t_pwl == NULL) {
			  printf("Memory Allocation Failed...\n");
			  exit(1);
		      }
		      
		      i_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (i_pwl == NULL) {
			  printf("Memory Allocation Failed...\n");
			  exit(1);
		      }
		    
		    for (j = 0; j<sum_numbers; j = j + 2) {   
		      
		      t_pwl[j>>1] = numbers[j];
		      if( strcmp(curr->posNode,"0")==0){
			i_pwl[j>>1] = -numbers[j+1];
		      }
		      else{
			i_pwl[j>>1] = numbers[j+1];
		      }      
		    }
		    
		    trans_value = PWL_func((int)sum_numbers/2, t_pwl, i_pwl, k*h, fin);    
		  }
		  else {
		    printf("Error in Transient Section\n");
		    exit(1);
		  }
		  if (curr->type == 'v' || curr->type == 'V'){
		      for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      tempB[curr3->thesiB] =  trans_value;//B[curr3->thesiB] + trans_value;
		  }else if ( (curr->type == 'i' || curr->type == 'I') ) {
		     for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      
		      /****************************************************************/
		      tempB[curr3->thesiB] = - trans_value; // B[curr3->thesiB] - trans_value; 
		      /***************************************************************/
		      if (curr3->thesiB_pairNode!=-1){
			  tempB[curr3->thesiB_pairNode] =  trans_value;		//B[curr3->thesiB_pairNode] + trans_value;	
			  
		      }  
		    
		  }
		  //if (numbers!=NULL)
		  //free(numbers);
	      }
	      
	  }

	  for (i=0; i<size; i++) {
	    temp_array[i] = temp_array[i] + tempB[i] + e_prev[i];
	    e_prev[i]=tempB[i];				//			<-----------------------
	    for (j=0; j<size; j++) {
	      patates[thesh] = A[i][j];
	      thesh++;
	    }
	  }
	  thesh = 0;
	  
	  m = gsl_matrix_view_array (patates, size, size); 	//m.matrix= A
	  b = gsl_vector_view_array (temp_array, size);		//b.vector= B

	  gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	  
	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	  for(i=0; i<size; i++){
	    x_prev[i]=gsl_vector_get (x, i);
	  }
	  
	  thesiPlot=0;
	  for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	      for(i=0; i<size; i++){
		  if ( strcmp(X[i].name,curPlot->node) == 0 ){
		    liseis[ thesiPlot ][k-1] = x_prev[i];
		  }
	      }
	      thesiPlot++;
	  }
	  
	}
	printf("TR\n");
	thesiPlot=0;
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot = curPlot->next){
	  
	      /***dimiourgia onomatos arxiou eksodou***/
	    sumPlots = strlen(curPlot->node) + 25;
	    PlotFileName = (char  *)malloc(sumPlots*sizeof( char ));
	    strcpy(PlotFileName,"set output \"TR_ \0");
	    PlotFileName[ strlen(PlotFileName)- 1] = curPlot->type;
	    strcat(PlotFileName, "(\0");
	    strcat(PlotFileName, curPlot->node);
	    strcat(PlotFileName, ").png\"\0");
	    
	      /***anoigma .png output file kai ektiposi tou apotelesmatos se auto***/
	    h_plot[thesiPlot] = gnuplot_init();
	    gnuplot_setstyle(h_plot[thesiPlot], "lines");
	    gnuplot_cmd(h_plot[thesiPlot], "set terminal png");
	    gnuplot_cmd(h_plot[thesiPlot], PlotFileName);
	    gnuplot_set_xlabel(h_plot[thesiPlot], "Time");
	    
	    if(curPlot->type=='I' || curPlot->type=='i') 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "I");
	    else 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "V");
	    
	    PlotFileName[ strlen(PlotFileName) - 5 ] = '\0';
	    gnuplot_plot_xy(h_plot[thesiPlot], time_vector, liseis[ thesiPlot ], (int)fin/h, &PlotFileName[15]);

	    gnuplot_close(h_plot[thesiPlot]);
	    thesiPlot++;
	    
	}	
	
}


cs *BE_A_sparse(double h, cs *G, cs *C, int size) {
	cs *A;
	
	A = (cs *)calloc(size, sizeof(cs));
	
	if (A == NULL) {
	   printf("Memory Allocation Failed...\n");
	   exit(1);
	}
	 
	A = cs_add(G, C, 1.0, 1.0/h);
	  
        return A;
  
}

 /*************************************************************************************
 *Uses the Backward Euler method to calculate the new x values. 		      *
 *Calculates the right hand side of the equation, changes the e() array according to  *
 *the transient values of the changing elements and uses LU decomposition to find     *
 *the new x.									      *
 *************************************************************************************/ 
void BE_b_sparse(double *B, double h, double fin, int size, tElements *head_elem, tSources *sources, cs *G, cs *C, tPlot *head_plot, tarrayX *X) {
	int i=0;
	double *x_prev;
	tElements *curr;
	tnumbers *curr2;
	tSources *curr3;
	tPlot *curPlot;
        int k,j,thesiPlot=0;
	int sum_numbers = 0;
	double *numbers = NULL;
	double *temp_array, *temp;
	double trans_value;
	double *t_pwl;
	double *i_pwl;
	double *time_vector;	
	int sumPlots=0;	
	char *PlotFileName;
	gnuplot_ctrl **h_plot;
	cs *A = BE_A_sparse(h, G, C, size);
	double **liseis;
	double *tempB = (double*)calloc(size,sizeof(double));
	css *S;
	csn *N;
	double *B_new;
	double *x;
	
	FILE * fp = fopen("results_trans.txt", "w");
	
	
	x_prev = (double *)calloc(size,sizeof(double));
	if (x_prev == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	temp = (double *)calloc(size,sizeof(double));
	if (temp == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	temp_array = (double *)calloc(size,sizeof(double));
	if (temp_array == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	x = (double *)calloc(size,sizeof(double));
	if (x == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	    sumPlots++;
	}
	
	h_plot = (gnuplot_ctrl **)calloc(sumPlots,sizeof(gnuplot_ctrl *));
	if (h_plot == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	
	liseis = (double **)calloc(sumPlots,sizeof(double *));
	if (liseis == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	else{
	    for(i=0; i<sumPlots; i++){
		liseis[i] = (double *)calloc( (int)(fin/h) , sizeof(double ) );
		if (liseis[i] == NULL) {
		    printf("Memory Allocation Failed...\n");
		    exit(1);
		}	  
	    }
	}
	
	time_vector = (double *)calloc( (int)(fin/h)+1, sizeof(double) );
	if (time_vector == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}		
	
	//for (i=0; i<(int)(fin/h); i++) time_vector[i] = h*i;  //TODO (fin/h) to the other functions, delete sxolia & prints
	for (i=0; i*h<fin; i++) { time_vector[i] = h*i; }

	
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
	
	S = cs_sqr(2, G, 0);
	N = cs_lu(G, S, 1);
      
//	cs_spfree(C);    
	
	//fprintf(fp,"Oi luseis tis eksiswsis me LU gia araious pinakes einai : \n");
	    
	cs_ipvec(N->pinv, B_new, x, size);
	cs_lsolve(N->L, x);
	cs_usolve(N->U, x);
	cs_ipvec(S->q, x, B_new, size); 
	
	
	/*
	gsl_matrix_view m = gsl_matrix_view_array (G, size, size); 	        //m.matrix= A
	gsl_vector_view b = gsl_vector_view_array (B, size);			//b.vector= B
	gsl_vector *x = gsl_vector_alloc (size);				//alloc(x)
	gsl_permutation *p = gsl_permutation_alloc (size);			//alloc(p)
	
	gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	*/
	for (i=0; i<size; i++) {		
	  for (j = 0; j < size; j++) {	//search for correct indexing after column permutations
		if (S->q[j] == i){
		   x_prev[i] = x[j];
		   break;
		}
	  }
	  //printf("x[%d] = %g\n", i, x_prev[i]);
	  //fprintf(fp, "x[%d] = %g\n", i, x[j]);
	}

	for(k=1; k*h<(fin+h); k++) {
	  for (i = 0; i < size; i++) {
	    temp_array[i] = 0;
	    temp[i] = 1.0/h * x_prev[i];
	  }
	      //for (j = 0; j < size; j++) {  
		  //temp_array[i] = temp_array[i] + C[i][j] * 1.0/h * x_prev[j];
	  cs_gaxpy(C, temp, temp_array);
		  //A[i][j] = G[i][j] + temp_elem;  
	     // }
	  
	  // calculate e(t) which changes at every step of transient analysis  
	  for (curr = head_elem; curr != NULL; curr = curr->next) {
	      if (curr->trans!=NULL) {
		  sum_numbers=0;
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    sum_numbers++;
		  }
		  numbers = (double *)malloc(sum_numbers*sizeof(double));
		  if (numbers == NULL) {
		      printf("Memory Allocation Failed...numbers\n");
		      exit(1);
		  }

		  sum_numbers = 0;
		  
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    numbers[sum_numbers] = curr2->number;
		    sum_numbers++;
		  }
		    
		   if (strcmp(curr->trans->method,"EXP") == 0 || strcmp(curr->trans->method,"exp") == 0) {
		    
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = EXP_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    else{
			trans_value = EXP_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    
		  }
		  else if (strcmp(curr->trans->method,"SIN") == 0 || strcmp(curr->trans->method,"sin") == 0) {
		    
		     if( strcmp(curr->posNode,"0")==0){
			trans_value = SIN_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }
		     else{
			trans_value = SIN_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }		      
		  }
		  else if (strcmp(curr->trans->method,"PULSE") == 0 || strcmp(curr->trans->method,"pulse") == 0) {
		    
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = PULSE_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }
		    else{ 
			trans_value = PULSE_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }  
		  }
		   else if (strcmp(curr->trans->method,"PWL") == 0 || strcmp(curr->trans->method,"pwl") == 0) {
		      t_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (t_pwl == NULL) {
			  printf("Memory Allocation Failed...t_pwl\n");
			  exit(1);
		      }
		      
		      i_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (i_pwl == NULL) {
			  printf("Memory Allocation Failed...i_pwl\n");
			  exit(1);
		      }
		    
		    for (j = 0; j<sum_numbers; j = j + 2) {   		      
		      t_pwl[j>>1] = numbers[j];
		      if( strcmp(curr->posNode,"0")==0){
			i_pwl[j>>1] = -numbers[j+1];
		      }
		      else{
			i_pwl[j>>1] = numbers[j+1];
		      }		      
		    }
		    
		    trans_value = PWL_func((int)sum_numbers/2, t_pwl, i_pwl, k*h, fin);
		    
		  }
		  else {
		    printf("Error in Transient Section\n");
		    exit(1);
		  }
		  if (curr->type == 'v' || curr->type == 'V'){
		      for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      tempB[curr3->thesiB] =  trans_value;//B[curr3->thesiB] + trans_value;
		  }else if ( (curr->type == 'i' || curr->type == 'I') ) {
		     for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      
		      /****************************************************************/
		      tempB[curr3->thesiB] = - trans_value; // B[curr3->thesiB] - trans_value; 
		      /***************************************************************/
		      if (curr3->thesiB_pairNode!=-1){
			  tempB[curr3->thesiB_pairNode] =  trans_value;		//B[curr3->thesiB_pairNode] + trans_value;	
			  
		      }  
		    
		  }
		  //if (numbers!=NULL)
		  //free(numbers);
	      }
	      
	  }

	  for (i=0; i<size; i++) {
	    temp_array[i] = temp_array[i] + tempB[i];
	    //for (j=0; j<size; j++) {
	     // patates[thesh] = A[i][j];
	     // thesh++;
	    //}
	  }
	  //thesh = 0;
	S = (css *) cs_calloc(1, sizeof(css));
	N = (csn *) cs_calloc(1, sizeof(csn));
	
	for (i = 0; i < size; i++){
	    B_new[i] =  temp_array[i];
	}
	
	S = cs_sqr(2, A, 0);
	N = cs_lu(A, S, 1); 
	
	//fprintf(fp,"Oi luseis tis eksiswsis me LU gia araious pinakes einai : \n");
	    
	cs_ipvec(N->pinv, B_new, x, size);
	cs_lsolve(N->L, x);
	cs_usolve(N->U, x);
	cs_ipvec(S->q, x, B_new, size); 

	  
//	  m = gsl_matrix_view_array (A, size, size); 	//m.matrix= A
//	  b = gsl_vector_view_array (temp_array, size);		//b.vector= B

//	  gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	  
//	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	  for (i=0; i<size; i++) {		
	  for (j = 0; j < size; j++) {	//search for correct indexing after column permutations
		if (S->q[j] == i){
		   x_prev[i] = x[j];
		   break;
		}
	  }
	  //printf("x[%d] = %g\n", i, x_prev[i]);
	  //fprintf(fp, "x[%d] = %g\n", i, x[j]);
	}
	  
	  thesiPlot=0;
	  for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	      for(i=0; i<size; i++){
		  if ( strcmp(X[i].name,curPlot->node) == 0 ){
		    liseis[ thesiPlot ][k-1] = x_prev[i];
		  }
	      }
	      thesiPlot++;
	  }
	  
	}
	printf("BE\n");
	thesiPlot=0;
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot = curPlot->next){
	  
	      /***dimiourgia onomatos arxiou eksodou***/
	    sumPlots = strlen(curPlot->node) + 25;
	    PlotFileName = (char  *)malloc(sumPlots*sizeof( char ));
	    strcpy(PlotFileName,"set output \"BE_ \0");
	    PlotFileName[ strlen(PlotFileName)- 1] = curPlot->type;
	    strcat(PlotFileName, "(\0");
	    strcat(PlotFileName, curPlot->node);
	    strcat(PlotFileName, ").png\"\0");
	    
	      /***anoigma .png output file kai ektiposi tou apotelesmatos se auto***/
	    h_plot[thesiPlot] = gnuplot_init();
	    gnuplot_setstyle(h_plot[thesiPlot], "lines");
	    gnuplot_cmd(h_plot[thesiPlot], "set terminal png");
	    gnuplot_cmd(h_plot[thesiPlot], PlotFileName);
	    gnuplot_set_xlabel(h_plot[thesiPlot], "Time");
	    
	    if(curPlot->type=='I' || curPlot->type=='i') 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "I");
	    else 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "V");
	    
	    PlotFileName[ strlen(PlotFileName) - 5 ] = '\0';
	    gnuplot_plot_xy(h_plot[thesiPlot], time_vector, liseis[ thesiPlot ], (int)(fin/h), &PlotFileName[15]);

	    gnuplot_close(h_plot[thesiPlot]);
	    
	    fprintf(fp,"\n\n%c(%s)\n",curPlot->type, curPlot->node);
	    for(i=0; i*h<fin; i++) {
	      fprintf(fp, "time: %10g  value: %10g\n", time_vector[i], liseis[ thesiPlot ][i]);     
	      fflush(fp);
	    }
	    
	    thesiPlot++;
	}	
	fclose(fp);
}

 /************************************************************************************************************************
 * Calculates the left hand side of the equation for the Trapezoidal (TR) method					 *
 ************************************************************************************************************************/
cs *TR_A_sparse(double h, cs *G, cs *C, int size) {
	
	cs *A;
	
	A = (cs *)calloc(size, sizeof(cs));
	
	if (A == NULL) {
	   printf("Memory Allocation Failed...\n");
	   exit(1);
	}
	 
	A = cs_add(G, C, 1.0, 2.0/h);
        return A;
}


void TR_b_sparse(double *B, double h, double fin, int size, tElements *head_elem, tSources *sources, cs *G, cs *C, tPlot *head_plot, tarrayX *X) {
	int i=0;
	double *x_prev;
	tElements *curr;
	tnumbers *curr2;
	tSources *curr3;
	tPlot *curPlot;
        int k,j,thesiPlot=0;
	int sum_numbers = 0;
	double *numbers = NULL;
	double *temp_array;
	double trans_value;
	double *t_pwl;	
	double *i_pwl;
	double *time_vector;	
	int sumPlots=0;	
	char *PlotFileName;
	gnuplot_ctrl **h_plot;
	cs *A=TR_A_sparse(h, G, C, size);
	double **liseis;
	cs *C_temp;
	css *S;
	csn *N;
	double *B_new;
	double *x;
	char s[30];
	FILE *fp;
	double *tempB = (double*)calloc(size,sizeof(double));
	double *e_prev = (double*)malloc(size*sizeof(double)); //			<--------------
	


	fp = fopen("results_trans.txt", "w");
	
	x = (double *)calloc(size,sizeof(double));
	if (x == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	x_prev = (double *)calloc(size,sizeof(double));
	if (x_prev == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	temp_array = (double *)calloc(size,sizeof(double));
	if (temp_array == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	    sumPlots++;
	}
	
	C_temp = (cs *)calloc(size, sizeof(cs));
	if (C_temp == NULL) {
	   printf("Memory Allocation Failed...\n");
	   exit(1);
	}
	
	h_plot = (gnuplot_ctrl **)malloc(sumPlots*sizeof(gnuplot_ctrl *));
	if (h_plot == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	
	liseis = (double **)calloc(sumPlots,sizeof(double *));
	if (liseis == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}
	else{
	    for(i=0; i<sumPlots; i++){
		liseis[i] = (double *)calloc( (int)(fin/h) , sizeof(double ) );
		if (liseis[i] == NULL) {
		    printf("Memory Allocation Failed...\n");
		    exit(1);
		}	  
	    }
	}

	time_vector = (double *)calloc( ((int)(fin/h)+1), sizeof(double) );
	if (time_vector == NULL) {
	    printf("Memory Allocation Failed...\n");
	    exit(1);
	}		
	//printf("%d\n",(int)(fin/h));
	//for (i=0; i<(int)(fin/h); i++) time_vector[i] = h*i;  //TODO (fin/h) to the other functions, delete sxolia & prints
	
	for (i=0; i*h<fin; i++){ time_vector[i] = h*i; }//printf("time vector %g\n", time_vector[i]);}

	for (i=0; i<size; i++) {
	  e_prev[i] = B[i];		//    					<-----------------
	}
	
	S = (css *) cs_calloc(1, sizeof(css));
	N = (csn *) cs_calloc(1, sizeof(csn));
	
	B_new = (double *)calloc(size,sizeof(double));
	if (B_new == NULL) {
	  printf("Memory allocation failed. Exiting...\n");
	  exit(1);
	}

	for (i = 0; i < size; i++){
	    B_new[i] =  B[i];
	}
	
	S = cs_sqr(2, G, 0);
	N = cs_lu(G, S, 1);
      //printf("hi\n");
//	cs_spfree(C);    
	
	//fprintf(fp,"Oi luseis tis eksiswsis me LU gia araious pinakes einai : \n");
	    
	cs_ipvec(N->pinv, B_new, x, size);
	cs_lsolve(N->L, x);
	cs_usolve(N->U, x);
	cs_ipvec(S->q, x, B_new, size); 
	
//	gsl_matrix_view m = gsl_matrix_view_array (G, size, size); 	//m.matrix= A
//	gsl_vector_view b = gsl_vector_view_array (B, size);			//b.vector= B
//	gsl_vector *x = gsl_vector_alloc (size);				//alloc(x)
//	gsl_permutation *p = gsl_permutation_alloc (size);			//alloc(p)
	
//	gsl_linalg_LU_decomp (&m.matrix, p, &s);				
//	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	//for(i=0; i<size; i++) {
	  //x_prev[i] = x[i];
	  //printf("x[%d]: %g\n",i,x[i]);
	//}
	for (i=0; i<size; i++) {		
	  for (j = 0; j < size; j++) {	//search for correct indexing after column permutations
		if (S->q[j] == i){
		   x_prev[i] = x[j];
		   break;
		}
	  }
	 // printf("x[%d] = %g\n", i, x[j]);
	  //fprintf(fp, "x[%d] = %g\n", i, x[j]);
	}
	
	
//printf("hi2\n");
	for(k=1; k*h<(fin+h); k++) {
	 // sprintf(s, "%d", k);
	  //strcat(s,"C_temp\0");

	  for (i = 0; i < size; i++) {
	    temp_array[i] = 0.0;
	   // x_prev[i] = -x_prev[i];
	  }
	      
	  C_temp = cs_add(G, C, 1.0, (-2.0)/h);
	  //if(k<3)
	  //cs_print(C_temp, s, 0);
	  cs_gaxpy (C_temp, x_prev, temp_array);
	  //for (i=0; i<size; i++) {
	  //printf("x[%d] = %g\n", i, x_prev[i]);
	  //}
	  //for (i=0;i<size;i++){
	 //  if ( k<3)
	  //printf("temp_array[%d]: %g\n",i,temp_array[i]);
	 // }
	  //cs_print(C_temp, s, 0);
	  //for ( i=0; i<size
	  // calculate e(t) which changes at every step of transient analysis  
	  for (curr = head_elem; curr != NULL; curr = curr->next) {
	      if (curr->trans!=NULL) {
		  sum_numbers=0;
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    sum_numbers++;
		  }

		  numbers = (double *)calloc(sum_numbers,sizeof(double));
		  if (numbers == NULL) {
		      printf("Memory Allocation Failed...\n");
		      exit(1);
		  }

		  sum_numbers = 0;
		  
		  for (curr2 = curr->trans->num; curr2 != NULL; curr2 = curr2->next) {
		    numbers[sum_numbers] = curr2->number;
		    sum_numbers++;
		  }

		   if (strcmp(curr->trans->method,"EXP") == 0 || strcmp(curr->trans->method,"exp") == 0) {
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = EXP_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    else{
			trans_value = EXP_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		    }
		    
		  }
		   else if (strcmp(curr->trans->method,"SIN") == 0 || strcmp(curr->trans->method,"sin") == 0) {
		     if( strcmp(curr->posNode,"0")==0){
			trans_value = SIN_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }
		     else{
			trans_value = SIN_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], k*h, fin);
		     }      
		  }
		  else if (strcmp(curr->trans->method,"PULSE") == 0 || strcmp(curr->trans->method,"pulse") == 0) {	    
		    if( strcmp(curr->posNode,"0")==0){
			trans_value = PULSE_func(-numbers[0], -numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }
		    else{ 
			trans_value = PULSE_func(numbers[0], numbers[1], numbers[2], numbers[3], numbers[4], numbers[5], numbers[6], k*h, fin);
		    }	    
		  }
		   else if (strcmp(curr->trans->method,"PWL") == 0 || strcmp(curr->trans->method,"pwl") == 0) {
		      t_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (t_pwl == NULL) {
			  printf("Memory Allocation Failed...\n");
			  exit(1);
		      }
		      
		      i_pwl = (double *)calloc((int)sum_numbers/2,sizeof(double));
		      if (i_pwl == NULL) {
			  printf("Memory Allocation Failed...\n");
			  exit(1);
		      }
		    
		    for (j = 0; j<sum_numbers; j = j + 2) {   
		      
		      t_pwl[j>>1] = numbers[j];
		      if( strcmp(curr->posNode,"0")==0){
			i_pwl[j>>1] = -numbers[j+1];
		      }
		      else{
			i_pwl[j>>1] = numbers[j+1];
		      }      
		    }

		    trans_value = PWL_func((int)sum_numbers/2, t_pwl, i_pwl, k*h, fin);    
		  }
		  else {
		    printf("Error in Transient Section\n");
		    exit(1);
		  }
		  if (curr->type == 'v' || curr->type == 'V'){
		      for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      tempB[curr3->thesiB] =  trans_value;//B[curr3->thesiB] + trans_value;
		  }else if ( (curr->type == 'i' || curr->type == 'I') ) {
		     for(curr3=sources; strcmp(curr3->elemName,curr->name)!=0; curr3=curr3->next);
		      
		      /****************************************************************/
		      tempB[curr3->thesiB] = - trans_value; // B[curr3->thesiB] - trans_value; 
		      /***************************************************************/
		      if (curr3->thesiB_pairNode!=-1){
			  tempB[curr3->thesiB_pairNode] =  trans_value;		//B[curr3->thesiB_pairNode] + trans_value;	
			  
		      }  
		    
		  }
		  //if (numbers!=NULL)
		  //free(numbers);
	      }
	      
	  }

	  for (i=0; i<size; i++) {
	    temp_array[i] = - temp_array[i] + tempB[i] + e_prev[i];
	    e_prev[i]=tempB[i];	
	  }
	  
//	  S = (css *) cs_calloc(1, sizeof(css));
	  //N = cs_calloc (1, sizeof (csn)); 
//	  N = (csn *) cs_calloc(1, sizeof(csn));
	  
/*	  x = (double *)calloc(size,sizeof(double));
	  if (x == NULL) {
	    printf("Memory allocation failed. Exiting...\n");
	    exit(1);
	  }
	  
	  B_new = (double *)calloc(size,sizeof(double));
	  if (B_new == NULL) {
	    printf("Memory allocation failed. Exiting...\n");
	    exit(1);
	  }
*/
	  for (i = 0; i < size; i++){
	      B_new[i] =  temp_array[i];
	  }
	  
	  S = cs_sqr(2, A, 0);
	  N = cs_lu(A, S, 1);
	
//	  cs_spfree(C);    
	  
	  //fprintf(fp,"Oi luseis tis eksiswsis me LU gia araious pinakes einai : \n");
	      
	  cs_ipvec(N->pinv, B_new, x, size);
	  cs_lsolve(N->L, x);
	  cs_usolve(N->U, x);
	  cs_ipvec(S->q, x, B_new, size); 
	    
//	  m = gsl_matrix_view_array (A, size, size); 	//m.matrix= A
//	  b = gsl_vector_view_array (temp_array, size);		//b.vector= B

//	  gsl_linalg_LU_decomp (&m.matrix, p, &s);				
	  
//	  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	  for (i=0; i<size; i++) {		
	    for (j = 0; j < size; j++) {	//search for correct indexing after column permutations
		  if (S->q[j] == i){
		    x_prev[i] = x[j];
		    break;
		  }
	    }
	   // printf("x[%d] = %g\n", i, x_prev[i]);
	  }
	  
	  thesiPlot=0;
	  for(curPlot=head_plot; curPlot!=NULL; curPlot=curPlot->next){
	     //printf("%s\n",curPlot->node);
	      for(i=0; i<size; i++){
		  if ( strcmp(X[i].name,curPlot->node) == 0 ){
		    //printf("%s\n",X[i].name);
		    liseis[ thesiPlot ][k-1] = x_prev[i];
		    break;
		    //printf("liseis[%d][%d]: %g\n",thesiPlot,k-1,liseis[ thesiPlot ][k-1]);
		  }
	      }
	      thesiPlot++;
	  }
	  
	}
	printf("TR\n");
	thesiPlot=0;
	
	for(curPlot=head_plot; curPlot!=NULL; curPlot = curPlot->next){

	      /***dimiourgia onomatos arxiou eksodou***/
	    sumPlots = strlen(curPlot->node) + 25;
	    PlotFileName = (char  *)malloc(sumPlots*sizeof( char ));
	    strcpy(PlotFileName,"set output \"TR_ \0");
	    PlotFileName[ strlen(PlotFileName)- 1] = curPlot->type;
	    strcat(PlotFileName, "(\0");
	    strcat(PlotFileName, curPlot->node);
	    strcat(PlotFileName, ").png\"\0");
	    
	      /***anoigma .png output file kai ektiposi tou apotelesmatos se auto***/
	    h_plot[thesiPlot] = gnuplot_init();
	    gnuplot_setstyle(h_plot[thesiPlot], "lines");
	    gnuplot_cmd(h_plot[thesiPlot], "set terminal png");
	    gnuplot_cmd(h_plot[thesiPlot], PlotFileName);
	    gnuplot_set_xlabel(h_plot[thesiPlot], "Time");
	    
	    if(curPlot->type=='I' || curPlot->type=='i') 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "I");
	    else 
	      gnuplot_set_ylabel(h_plot[thesiPlot], "V");
	    
	    PlotFileName[ strlen(PlotFileName) - 5 ] = '\0';
	    gnuplot_plot_xy(h_plot[thesiPlot], time_vector, liseis[ thesiPlot ], (int)(fin/h), &PlotFileName[15]);

	    gnuplot_close(h_plot[thesiPlot]);

	    fprintf(fp,"\n\n%c(%s)\n",curPlot->type, curPlot->node);
	    for(i=0; i*h<fin; i++) {
	      fprintf(fp, "time: %10g  value: %10g\n", time_vector[i], liseis[ thesiPlot ][i]);     
	      fflush(fp);
	    }
	    thesiPlot++;
	}
	
	fclose(fp);
	
}

