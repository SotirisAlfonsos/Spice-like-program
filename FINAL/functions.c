/*		Circuit Simulation 2014 3 Solution of MNA System
 * 
 * Onoma programmatisti:	Alfonsos Swthrhs 959
 * 				Servou Katerina 852
 * 				Stauropoulos Antonis 1141
 * 				Vakerlis Christos 982
 * 
 * Hmeromhnia:		01/11/2014
 * 
 * ayto to arxeio xrhsimopoieitai apo to arxeio cs_parser.c kai periexei tis
 * synarthseis pou xrisimopoiei to programma
 */

/******************* 
 * Purpose :  Arxikopoihsh twn listwn 
 * Postconditions :  Epistrefetai i kefali tis listas arxikopoihmeni se NULL
    
 ******************/

#include "decompositions.c"
#include "decompositions2.c"
#include "cs_LU_Cholesky.c"
#include "cs_CG_BiCG.c"
#include "transient_analysis.c"
//#include "petros.c"

tElements *listElem_init() {
    tElements *head;
    head = NULL;
    return(head);
}

tDiodes *listDiodes_init() {
    tDiodes *head;
    head = NULL;
    return(head);
}

tMOS *listMOS_init() {
    tMOS *head;
    head = NULL;
    return(head);
}

tBJT *listBJT_init() {
    tBJT *head;
    head = NULL;
    return(head);
}

tPlot *listPlot_init() {
    tPlot *head;
    head = NULL;
    return(head);
}


/****************
 * Purpose :  Eisagwgi stoixeiwn stis listes kai dhmiourgia twn komvwn twn 4 listwn
 * Parameters :  head -> kefali tis listas
 *               counter -> metrhths pou dhlwnei to stoixeio pou eisagetai kathe fora ston komvo tis listas
 *               stoixeio -> to token pou diavazetai apo to netlist sti main
 * Preconditions :  Thewroume oti to netlist einai swsto kai ta stoixeia pou diavazontai antistoixoun se pragmatikes times
 * Postconditions : Epistrefetai kathe fora to root tis listas, wste na ginei i eisagwgh twn epomenwn stoixeiwn tou komvou h na dhmiourghthei o
 * 		     epomenos komvos kapoias listas.
 *****************/

tElements *insertElem(tElements *head, int counter, char *stoixeio) {		// V, I, L, R, C
    tElements *curr;
    tnumbers *cur_num; 
    double inserted;
    
    if (counter==1) {				//eisagwgi tou type kai tou name
      
      curr = (tElements *)malloc(sizeof(tElements));
      
      if (curr==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      curr->next = head;
      head = curr;
      head->type = stoixeio[0];
      head->trans=NULL;		//12/12

      strcpy(head->name,stoixeio);
      head->name[strlen(stoixeio)] = '\0';
    }
    else if (counter==2) {			//eisagwgi tou thetikou akrodekti
      strcpy(head->posNode,stoixeio);
      head->posNode[strlen(stoixeio)] = '\0';
    }
    else if (counter==3) {			//eisagwgi tou arnhtikou akrodekti
      strcpy(head->negNode,stoixeio);
      head->negNode[strlen(stoixeio)] = '\0';
    }
    else if (counter==4) {			//eisagwgi tou value
      inserted = atof(stoixeio);    
      head->value = inserted;
    }
    else if (counter==5) {			//eisagwgi tou source method 12/12
      head->trans = (tTrans_source *)malloc(sizeof(tTrans_source));
      if (head->trans==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      head->trans->method = strdup( stoixeio );
      head->trans->num = NULL;
    }
    else if (counter == 6 ) {			//eisagwgi twn arithmon tou kathe transient method 12/12
      cur_num = (tnumbers *)malloc(sizeof(tnumbers));
      if (cur_num==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      cur_num->number = atof( stoixeio );
      
      head->trans->num = cur_num;
      head->trans->num->next=NULL;
      
    }
    else if (counter > 6) {			//eisagwgi twn arithmon tou kathe transient method 12/12
      for (cur_num = head->trans->num; cur_num->next!=NULL; cur_num = cur_num->next);
      
      
      cur_num->next = (tnumbers *)malloc(sizeof(tnumbers));
      if (cur_num==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      cur_num->next->number = atof( stoixeio );
      cur_num->next->next=NULL;
     
    }
    
    return(head);
}   


tDiodes *insertDiodes(tDiodes *head, int counter, char *stoixeio) {		//DIODES
    tDiodes *curr;
    double inserted;
    
    if (counter==1) {
	curr = (tDiodes *)malloc(sizeof(tDiodes));
    
      if (curr==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      curr->next = head;
      head = curr;
      strcpy(head->name,stoixeio);
      head->name[strlen(stoixeio)] = '\0';
    }
    else if (counter==2) {
      strcpy(head->posNode,stoixeio);
      head->posNode[strlen(stoixeio)] = '\0';
    }
    else if (counter==3) {
      strcpy(head->negNode,stoixeio);
      head->negNode[strlen(stoixeio)] = '\0';
    }
    else if (counter==4) {
      strcpy(head->modName,stoixeio);
      head->modName[strlen(stoixeio)] = '\0';
    }
    else if (counter==5){
      inserted = atof(stoixeio);
      
      head->area = inserted;
    }
    return(head);
} 


tMOS *insertMOS(tMOS *head, int counter, char *stoixeio) {			//MOS
    tMOS *curr;
    double inserted;
    //char *end;
    
    if (counter==1) {
      
      curr = (tMOS *)malloc(sizeof(tMOS));
      if (curr==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      curr->next = head;
      head = curr;
      strcpy(head->name,stoixeio);
      head->name[strlen(stoixeio)] = '\0';
    }
    else if (counter==2) {
      strcpy(head->D,stoixeio);
      head->D[strlen(stoixeio)] = '\0';
    }
    else if (counter==3) {
      strcpy(head->G,stoixeio);
      head->G[strlen(stoixeio)] = '\0';
    }
    else if (counter==4) {
      strcpy(head->S,stoixeio);
      head->S[strlen(stoixeio)] = '\0';
    }
    else if (counter==5){
      strcpy(head->B,stoixeio);
      head->B[strlen(stoixeio)] = '\0';
    }
    else if (counter==6) {
      strcpy(head->modName,stoixeio);
      head->modName[strlen(stoixeio)] = '\0';
    }
    else if (counter==7) {
      inserted = atof(stoixeio);
      head->L = inserted;
    }
    else if (counter==8) {
      inserted = atof(stoixeio);
      head->W = inserted;
    }
    return(head);
} 


tBJT *insertBJT(tBJT *head, int counter, char *stoixeio) {			//BJT
    tBJT *curr;
    double inserted;
    
    if (counter==1) {
      curr = (tBJT *)malloc(sizeof(tBJT)); 
      if (curr==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
      }
      
      curr->next = head;
      head = curr;    
      strcpy(head->name,stoixeio);
      head->name[strlen(stoixeio)] = '\0';
    }
    else if (counter==2) {
      strcpy(head->C,stoixeio);
      head->C[strlen(stoixeio)] = '\0';
    }
    else if (counter==3) {
      strcpy(head->B,stoixeio);
      head->B[strlen(stoixeio)] = '\0';
    }
    else if (counter==4) {
      strcpy(head->E,stoixeio);
      head->E[strlen(stoixeio)] = '\0';
    }
    else if (counter==5) {
      strcpy(head->modName,stoixeio);
      head->modName[strlen(stoixeio)] = '\0';
    }
    else if (counter==6) {
      inserted = atof(stoixeio);
      head->area = inserted;
    }
    return(head);
} 


tPlot *insertPlot(tPlot *head, char *plot_elem, char initial_char) {
    tPlot *curr;
    
    curr = (tPlot *)malloc(sizeof(tPlot));
    if (curr==NULL) {
	printf("Memory Allocation Error\n");
	exit(1);
    }
    
    curr->next = head;
    head = curr;
     
     head->type = initial_char;
     head->node = strdup(plot_elem);
     
     return(head);
}
     
tSources *insertSourse(tSources *head, char *name, int node, int pairNode){
  tSources *curr;//,*curr2;
  
  curr = (tSources *)malloc(sizeof(tSources));
  if (curr == NULL) {
      printf("Memory Allocation Error\n");
      exit(1);
  }
  
  curr->elemName = strdup(name);
  curr->thesiB = node;
  curr->thesiB_pairNode = pairNode;
 // curr->next=NULL;
  
/*  if (head == NULL){
    head = curr;
    curr->next = NULL;
  }
  else{
    for(curr2=head; curr2->next!=NULL; curr2=curr2->next);
    curr2->next=curr;
    curr->next = NULL;
  }
  
  return head;
*/  
//allagh 2015
    curr->next = head;
    return curr;
}
    
     
/*tNodes *DCNodes(tElements *head_elem){		//vriskei posoi DC nodes iparxoun  kai toys antistoixei me diaforetiko id ton kathena
  	tElements *curElem;			
  	tNodes *head_tNodes = NULL;
	tNodes *cur_tNodes,*cur_tNodes2;
	long int id = 0;
	int flag;
	
  	for (curElem = head_elem; curElem != NULL; curElem = curElem->next) {  
	  
	    flag = 0;
	    for (cur_tNodes = head_tNodes; cur_tNodes != NULL; cur_tNodes = cur_tNodes->next) {  
		if (strcmp((curElem->posNode),"0")==0 || strcmp((curElem->posNode),cur_tNodes->name)==0) {
		    flag=1;
		    break;
		}
	    }
	    if (flag==0 && strcmp((curElem->posNode),"0")!=0) {
	      if (head_tNodes==NULL) {
		  if ((head_tNodes = (tNodes *)malloc(sizeof(tNodes))) == NULL) {
		      return NULL;
		  }
		  head_tNodes->next = NULL;
		  if ((head_tNodes->name = strdup(curElem->posNode)) == NULL) {
		      return NULL;
		  }
		  head_tNodes->id = id++;
	      }  
	      else {
		  if ((cur_tNodes = (tNodes *)malloc(sizeof(tNodes))) == NULL) {
		      return NULL;
		  }
		  cur_tNodes->next = NULL;
		  if ((cur_tNodes->name = strdup(curElem->posNode)) == NULL) {
		      return NULL;
		  }
		  cur_tNodes->id = id++; 
		  for (cur_tNodes2 = head_tNodes; cur_tNodes2->next!=NULL; cur_tNodes2=cur_tNodes2->next);
		  cur_tNodes2->next = cur_tNodes;
	      }
	    }  
	    
	    flag=0;
	    for (cur_tNodes = head_tNodes; cur_tNodes!=NULL; cur_tNodes=cur_tNodes->next) {  
	      if (strcmp((curElem->negNode),"0")==0 || strcmp((curElem->negNode),cur_tNodes->name)==0) {
		  flag = 1;
		  break;
	      }
	    }
	    if (flag==0 && strcmp((curElem->negNode),"0")!=0){
	      if (head_tNodes==NULL) {
		  if ((head_tNodes = (tNodes *)malloc(sizeof(tNodes))) == NULL) {
		      return NULL;
		  }
		  head_tNodes->next = NULL;
		  if ((head_tNodes->name = strdup(curElem->negNode)) == NULL) {
		      return NULL;
		  }
		  head_tNodes->id = id++;
	      }  
	      else {
		  if ((cur_tNodes = (tNodes *)malloc(sizeof(tNodes))) == NULL) {
		      return NULL;
		  }
		  cur_tNodes->next=NULL;
		  if ((cur_tNodes->name = strdup(curElem->negNode)) == NULL) {
		      return NULL;
		  }
		  cur_tNodes->id = id++;
		  for (cur_tNodes2 = head_tNodes; cur_tNodes2->next!=NULL; cur_tNodes2 = cur_tNodes2->next);
		  cur_tNodes2->next = cur_tNodes;
	      }
	    }
	}
	return head_tNodes;
}

*/

/* Create a new hashtable. */
hashtable_t *ht_create(long int size) {
 
	hashtable_t *hashtable = NULL;
	int i;
 
	if (size < 1) return NULL;
 
	/* Allocate the table itself. */
	if ((hashtable = malloc(sizeof(hashtable_t))) == NULL) {
		return NULL;
	}
 
	/* Allocate pointers to the head nodes. */
	if ((hashtable->table = malloc(sizeof(entry_t *) * size)) == NULL) {
		return NULL;
	}
	for (i = 0; i < size; i++) {
		hashtable->table[i] = NULL;
	}
 
	hashtable->size = size;
	return hashtable;	
}
 

/* Create a key-value pair. */
entry_t *ht_newpair(char *key, double value , char *name, char* pairNodeName) {
	entry_t *newpair;
 
	if ((newpair = malloc(sizeof(entry_t))) == NULL) {
		return NULL;
	}
 
	if ((newpair->key = strdup(key)) == NULL) {
		return NULL;
	}
	
	if ((newpair->name = strdup(name)) == NULL) {
		return NULL;
	}

	if (pairNodeName!=NULL){						//vale to onoma tou pairnode   OR   NULL ean einai to "0"
		if ( (newpair->pairNodeName = strdup(pairNodeName)) == NULL ) {
			return NULL;
		}
	}
	else{
		newpair->pairNodeName = NULL;
	}
	
	newpair->pairNode = -1;
	newpair->value =  value;
	newpair->next = NULL;
 
	return newpair;
}
 
int ht_hash( hashtable_t *hashtable, char *key ) {

	unsigned long int hashval=101;
	int c = 0;

	while ((c = *key++))
		hashval = c + (hashval << 8) + (hashval << 7) - hashval;  //((hashval << 8) + hashval) + c;
	
	return hashval % hashtable->size;
} 
 
/* Insert a key-value pair into a hash table. */
void ht_set( hashtable_t *hashtable, char *key, double value , char *name , long int id , char * pairNodeName, hashtable_t *remainHash) {
	entry_t *newpair = NULL;
	entry_t *next = NULL;
//	entry_t *last = NULL;
	int i=0;
//	printf("mpika\n");
	next = hashtable->table[ id ];
 
	if (next==NULL){		//ean pseydes simenei oti exei ginei lathos tou hashing function: 2 nodes same key
//		printf("mpika\n");  
//	    while( next != NULL && next->key != NULL && strcmp( key, next->key ) > 0 ) {
//		    last = next;
//		    next = next->next;
//	    }
    
	    newpair = ht_newpair( key, value , name , pairNodeName);

//	    if( next == hashtable->table[ id ] ) {		/* We're at the start of the linked list in this bin. */
	    newpair->next = next;
	    hashtable->table[ id ] = newpair;
//	    } 
//	    else if ( next == NULL ) {			/* We're at the end of the linked list in this bin. */
//		    last->next = newpair;
//	    } 
//	    else  {						/* We're in the middle of the list. */
//		    newpair->next = next;
//		    last->next = newpair;
//	    }
//	  	printf("mpika\n");  
	}
	else if( strcmp(next->key,key)==0 ) {
	    newpair = ht_newpair( key, value , name , pairNodeName);
	    newpair->next = next;
	    hashtable->table[ id ] = newpair;
	}
	else{							//dimiourgia enos 2ou hashtable gia ta lathi tou hash function
	    while (remainHash->table[i]!=NULL){
		if  ( strcmp( remainHash->table[i]->key, key )==0 )
		      break;
		i++;
	    }
	    next = remainHash->table[i];
	    
	    newpair = ht_newpair( key, value , name , pairNodeName);
	    
	    newpair->next = next;
	    remainHash->table[i] = newpair;	    
	  
	}
}

int getHashId( hashtable_t *hashtable, hashtable_t *remainHash, tNodes **negList, char *pairNodeName, int *k, hashtable_t* misHash){
	int i = 0,a;
	tNodes *cur_tNodes2 = *negList;
	tNodes *cur_tNodes, *temp = *negList;
	entry_t *newpair = NULL;
	entry_t *next = NULL;
		
	i = ht_hash( hashtable, pairNodeName );

	if ( hashtable->table[i] == NULL) goto L1;	//gia opoion dusmoiro programmatisti exei tin kaki tuxi na diavasei auton ton kwdika:
							//na kserete oti auto to teratourghma einai ergo XRHSTOU
	if ( strcmp( hashtable->table[i]->key, pairNodeName )==0){			//search in hashtable
		return hashtable->table[i]->id;
	}
	else {
		i=0;
		while (remainHash->table[i]!=NULL){					//search in remainHash
		  	if ( strcmp( remainHash->table[i]->key, pairNodeName )==0){
				return remainHash->table[i]->id;
			}
			i++;
		}
		
		L1:	  					//node has never found as posNode
//		printf("xs %s\n", pairNodeName);
		
		i = ht_hash( misHash, pairNodeName );
		next = misHash->table[i];
		if ( next == NULL) {
			newpair = ht_newpair( pairNodeName, 0 , "\0" ,  "\0" );
			newpair->id = *k;
			
			newpair->next = NULL;
			misHash->table[ i ] = newpair;
			    
			*k=*k+1;
			return *k - 1;

		}
		else if( strcmp(next->key, pairNodeName)==0 ) {
			return next->id;
		}
		
		

		while (cur_tNodes2 != NULL){						//search in remainHash 
//			printf("%c\t",cur_tNodes2->name[0]);//, pairNodeName);
			a=strcmp( cur_tNodes2->name, pairNodeName );
			if (a == 0){
				return cur_tNodes2->id;
			}
			else if (a>0)
				break;
			
			temp = cur_tNodes2;
			cur_tNodes2 = cur_tNodes2->next;
		}
		//printf("\n");
		
		if ((cur_tNodes = (tNodes *)malloc(sizeof(tNodes))) == NULL) {		//has not been found so insert and return new id k++
		      return NULL;
		}
		if ((cur_tNodes->name = strdup(pairNodeName)) == NULL) {
		      return NULL;
		}
		cur_tNodes->id = *k;
		
		cur_tNodes->next = cur_tNodes2;
		if (temp==NULL)
		      *negList = cur_tNodes;
		else if (strcmp( (*negList)->name, pairNodeName )>0)
		      *negList = cur_tNodes;
		else
		      temp->next = cur_tNodes;
		//*negList = cur_tNodes;

		*k=*k+1;
		return *k - 1;
	}
	
}  
/*	int i = 0;
	tNodes *cur_tNodes = *negList;
	i = ht_hash( hashtable, pairNodeName );

	if ( hashtable->table[i] == NULL) goto L1;
	
	if ( strcmp( hashtable->table[i]->key, pairNodeName )==0){			//search in hashtable
		return hashtable->table[i]->id;
	}
	else {
		i=0;
		while (remainHash->table[i]!=NULL){					//search in remainHash
		  	if ( strcmp( remainHash->table[i]->key, pairNodeName )==0){
				return remainHash->table[i]->id;
			}
			i++;
		}
		
		L1:	  
		while (cur_tNodes != NULL){						//search in remainHash
//			printf("mpika %s %s\n",cur_tNodes->name, pairNodeName);
			if ( strcmp( cur_tNodes->name, pairNodeName ) == 0){
				return cur_tNodes->id;
			}  
			cur_tNodes = cur_tNodes->next;
		}
		
		if ((cur_tNodes = (tNodes *)malloc(sizeof(tNodes))) == NULL) {		//has not been found so insert and return new id k++
		      return NULL;
		}
		cur_tNodes->next = *negList;
		if ((cur_tNodes->name = strdup(pairNodeName)) == NULL) {
		    return NULL;
		}
		cur_tNodes->id = *k;
		
		*negList = cur_tNodes;

		*k=*k+1;
		return *k - 1;
	}
	
}*/

hashtable_t *sumNodes(tElements *head_elem,int flag_cholesky, dc_analysis *dc_table, tPlot *head_plot,
			  int flag_biCG, int flag_CG ,double ITOL, int flag_sparse,double *trans_h_fin,int trans_method,int lines){
	int i,flag,j=0,l=0,k=0,b=0;
	flag=0;
	int posNode,negNode;
	int size=0;
	int *validHash;
	int valids = 0;
	int misHashSize = 0;
  	tElements *curElem;			//voithitikoi deiktes gia prospelasi twn listwn	
	tNodes *negList=NULL, *cur_tNodes;
	entry_t* ht;
	
	dcIt dcInodes = {-1,-1,0};								//14/11
	
	int m1=0;
	int m2=0;
	double **A,**C_trans;
	double *B;
	tarrayX *X;
	double apol;
	int thesiDC = -1;
	
	int nz = 0, nz2=0;				//21/11
	cs *Asp, *C;				//21/11
	cs *C_trans_sparce,*C_trans_sparce_compressed;				//21/11
	int spCoun = 0;				//21/11
	int spTrCoun = 0;
	int piknos = flag_sparse ? 0 : 1;				//21/11 
	double sum;
	int flag_trans = (trans_h_fin[0]==0 && trans_h_fin[1]==0) ? 0 : 1;
	tSources *sourcesRoot;
  
	sourcesRoot = NULL;
	
	
	/*** create&set hashtable, find m1&m2   ***/
	printf("/*** create & set hashtable  , calc m1 & m2   ***/\n");
	hashtable_t *hashtable = ht_create( lines<<6); //prosoxi tis olisthisis: analogi tis Ram tou pc!!!!!
	hashtable_t *remainHash = ht_create( lines );
	hashtable_t *misHash = ht_create( lines<<2 );
  	for (curElem = head_elem; curElem!=NULL; curElem=curElem->next) {  
	  
	    if(curElem->name[0]=='r' || curElem->name[0]=='R' || curElem->name[0]== 'i' || curElem->name[0]== 'I'){
	      m1++;
	    }else if(curElem->name[0]== 'l' || curElem->name[0]== 'L' || curElem->name[0]== 'v' || curElem->name[0]== 'V') {
	      m2++;
	    }
	    
	    posNode = ht_hash(hashtable, curElem->posNode );

	    if (strcmp(curElem->negNode,"0")==0 ){
		ht_set( hashtable, curElem->posNode, curElem->value, curElem->name , posNode, NULL, remainHash);		
	    }else if (strcmp(curElem->posNode,"0")==0 ){
		negNode = ht_hash(hashtable, curElem->negNode );
		ht_set( hashtable, curElem->negNode, -curElem->value, curElem->name , negNode, NULL, remainHash);
	    }else{
		ht_set( hashtable, curElem->posNode, curElem->value, curElem->name , posNode, curElem->negNode, remainHash);
	    }
	   
	}
	
	if(!flag_trans){    /***  free head_elem list  ***/
	  for (curElem = head_elem; curElem!=NULL; curElem=head_elem) {		
	    head_elem = head_elem->next;
	    free(curElem);
	  }
	}
	
	
	printf("\n/*** Create table with hashtable's valid places and calculation of n-1 ***/\n");
	
	for( i=0; i<(hashtable->size); i++ ){
	    if ( hashtable->table[i] != NULL ){
		size++;
	    }
	}
	
	validHash = (int *)malloc( size * sizeof(int) );
	
	for( i=0; i<(hashtable->size); i++ ){
	    if ( hashtable->table[i] != NULL ){
		hashtable->table[i]->id = valids;
		validHash[valids++] = i;	

	    }
	}
	
	i = 0;
	
	while (remainHash->table[i]!=NULL){
	    remainHash->table[i]->id = i+size;
	    i++;
	}
	size = size + i;   
	
	
	k=size;
	printf("/*** eisgogi twn negNodes IDs sto hashtable kai\n     dhmiourgia IDs gia negNodes pou pou den vrethikan\n     os posNode i den eixan posNode==0 ***/\n");  
	for(j=0; j<size; j++){
		if (j < valids){
		    ht = hashtable->table[ validHash[j]  ];
		}
		else{  
		    ht = remainHash->table[ j - valids];
		}		
		//if(j%10000==0) printf("%d\t",size-j);
		while (ht != NULL){    
		    if ( ht->pairNodeName != NULL){
			ht->pairNode  = getHashId( hashtable, remainHash, &negList, ht->pairNodeName, &k, misHash);	 
		    }		  
		    ht= ht->next;
		}

	}
	
	for (j=0; j<misHash->size; j++){
	  if( misHash->table[j] != NULL ) misHashSize++;
	}
	printf("\nm1: %d\tm2: %d\tn-1: %d \t= { valids: %d\tremainHash places: %d\tHashed negnodes: %d\tunregisterd negNodes: %d }\tlines: %d\n\n",m1,m2,k,valids,i,misHashSize,k-size-misHashSize,lines);
	
	
	
	/** piknoi + araioi pinakes **/				//21/11 
	/*** create & init A , B and X tables  ***/
	printf("/*** create & init A , B and X tables  ***/\n");
	B=(double *)malloc((m2+k)*sizeof(double));
	if (B==NULL){
		printf("No memory available. Exiting...");
		exit (1);
	}
	for(i=0;i<(m2+k);i++){
	    B[i]=0;
	}
	
			/******** aloc mem for A ********/
	if (piknos){				//21/11 
	    A = (double **)malloc((m2+k)*sizeof(double *));
	    if (A!=NULL) {
		    for (i=0; i<=(m2+k); i++) {
			    A[i] = (double *)malloc((m2+k)*sizeof(double));
			    if (A[i]==NULL) {
				    printf("No memory available. Exiting...");
				    exit (1);
			    }
		    }
	    }
	    else  {
		    printf("No memory available. Exiting...");
		    exit (1);
	    }
	    for(i=0;i<(m2+k);i++){
	      for(j=0;j<(m2+k);j++){
		A[i][j]=0;
	      }
	    }
	    
	    if(flag_trans){
	    	  C_trans = (double **)malloc((m2+k)*sizeof(double *));
		  if (C_trans!=NULL) {
			  for (i=0; i<=(m2+k); i++) {
				  C_trans[i] = (double *)malloc((m2+k)*sizeof(double));
				  if (C_trans[i]==NULL) {
					  printf("No memory available. Exiting...");
					  exit (1);
				  }
			  }
		  }
		  else  {
			  printf("No memory available. Exiting...");
			  exit (1);
		  }
		  for(i=0;i<(m2+k);i++){
		    for(j=0;j<(m2+k);j++){
		      C_trans[i][j]=0;
		    }
		  }
	      
	    }
	    
	}
	else{		/*****    cacl nz    *****/	//sparse table		//21/11 
	    printf("/***    cacl nz    ***/\n");
	    for( i=0; i<size; i++ ){
	      
		if (i < valids){
		    ht = hashtable->table[ validHash[i]  ];
		}
		else{  
		    ht = remainHash->table[ i - valids];
		}
			
		for( ht = ht; ht!=NULL; ht=ht->next ){
		  if( ht->name[0] == 'r' || ht->name[0] == 'R' ){			//dhmiourgia tou n-1 x n-1 aristerou pano merous tou A
			if( ht->pairNode != -1 ){
			  nz+=4;
			}else{
			  nz++;
			}
		  }
		  else if ( ht->name[0] == 'v' || ht->name[0] == 'V' || ht->name[0] == 'l' || ht->name[0] == 'L') {
			nz+=2;
			if (ht->pairNode!=-1) nz+=2;
			if ( (ht->name[0] == 'l' || ht->name[0] == 'L') && flag_trans ) nz2++;
		  }
		  else if ( (ht->name[0] == 'C'|| ht->name[0] == 'c') && flag_trans ) {				//transient if C
			if( ht->pairNode != -1 ){
			      nz2 += 4;
			}else{
			      nz2++;
			}     
		  }
		}
	    }

	    Asp = cs_spalloc(m2+k, m2+k, nz, 1, 1);			//21/11 
	    Asp->nz = nz;
	    if (flag_trans){
		  C_trans_sparce = cs_spalloc(m2+k, m2+k, nz2, 1, 1);			//21/11 
		  C_trans_sparce->nz = nz2;
	    }
	}
	
	X=(tarrayX *)malloc((m2+k)*sizeof(tarrayX));
	if (X==NULL){
	    printf("No memory available. Exiting...");
	    exit (1);
	}
	
	
	/*** set A ,B and X tables  ***/
	printf("/*** set A ,B and X tables  ***/\n");

	cur_tNodes = negList;					//set X from negList
	while(cur_tNodes != NULL){ //me ==
	  if ( (X[cur_tNodes->id].name = strdup(cur_tNodes->name )) == NULL) {
	      return NULL;
	  } 
	  cur_tNodes=cur_tNodes->next;
	}
	for (j=0; j<misHash->size; j++){
	  if( misHash->table[j] != NULL ){
	    if ( (X[ misHash->table[j]->id  ].name = strdup(misHash->table[j]->key )) == NULL) {
	      return NULL;
	    } 
	    free( misHash->table[j] );
	  }
	}
	 free( misHash );
	
	for( i=0; i<size; i++ ){
	  
	  if (i < valids){
	      ht = hashtable->table[ validHash[i]  ];
	  }
	  else{  
	      ht = remainHash->table[ i - valids];
	  }


	  if (ht==NULL)printf("SOS SOS SOS************************mpika vlaka************************SOS SOS SOS\n");//{					
	  if ((X[i].name = strdup(ht->key)) == NULL) {		//X[i].Name=...
	      return NULL;
	  }

	  
	  for( ht = ht; ht!=NULL; ht=ht->next ){	 
	    if( ht->name[0] == 'r' || ht->name[0] == 'R' ){			//dhmiourgia tou n-1 x n-1 aristerou pano merous tou A
	      
	      if((ht->value)>=0)
		apol = ht->value;
	      else 
		apol = -ht->value;
	      	  	
	      if (piknos){
		  if( ht->pairNode != -1 ){
		    A[i][i] += 1.0/apol;
		    A[ht->pairNode][ht->pairNode] += 1.0/apol;
		    A[ht->pairNode][i] -= 1.0/apol;
		    A[i][ht->pairNode] -= 1.0/apol;
		  }else{
		    A[i][i] += 1.0/apol; 
		  }
	      }else{			//sparse table			//21/11 
		  if( ht->pairNode != -1 ){

		    Asp->i[spCoun] = i;
		    Asp->p[spCoun] = i;
		    Asp->x[spCoun] = 1.0/apol;
		    spCoun++;
		    	      
		    Asp->i[spCoun] = ht->pairNode;
		    Asp->p[spCoun] = ht->pairNode;
		    Asp->x[spCoun] = 1.0/apol;
		    spCoun++; 
		    
		    Asp->i[spCoun] = ht->pairNode;
		    Asp->p[spCoun] = i;
		    Asp->x[spCoun] = -1.0/apol;
		    spCoun++;
		    
		    Asp->i[spCoun] = i;
		    Asp->p[spCoun] = ht->pairNode;
		    Asp->x[spCoun] = -1.0/apol;
		    spCoun++;
		    
		  }else{
		    
		    Asp->i[spCoun] = i;
		    Asp->p[spCoun] = i;
		    Asp->x[spCoun] = 1.0/apol;
		    spCoun++;
		    
		  }
		
	      }
	    }
	    else if ( ht->name[0] == 'v' || ht->name[0] == 'V' || ht->name[0] == 'l' || ht->name[0] == 'L') {
	      if (piknos){
		  A[i][k+b] += 1;
		  if (ht->pairNode!=-1) A[ht->pairNode][k+b] -= 1;
		  
		  A[k+b][i] += 1;
		  if (ht->pairNode!=-1) A[k+b][ht->pairNode] -= 1;
	      }
	      else{			//sparse table			//21/11 
		  Asp->i[spCoun] = i;
		  Asp->p[spCoun] = k+b;
		  Asp->x[spCoun] = 1;
		  spCoun++;
		  if (ht->pairNode!=-1){
		    Asp->i[spCoun] = ht->pairNode;
		    Asp->p[spCoun] = k+b;
		    Asp->x[spCoun] = -1;
		    spCoun++; 
		  }
		  Asp->i[spCoun] = k+b;
		  Asp->p[spCoun] = i;
		  Asp->x[spCoun] = 1;
		  spCoun++;
		  if (ht->pairNode!=-1){
		    Asp->i[spCoun] = k+b;
		    Asp->p[spCoun] = ht->pairNode;
		    Asp->x[spCoun] = -1;
		    spCoun++;		  
		  }
	      }
	      	  	
	      if (ht->name[0] == 'v' || ht->name[0] == 'V'){
		B[k+b]=ht->value;
		if(flag_trans){
		    sourcesRoot = insertSourse(sourcesRoot, ht->name, b+k, -1);
		    //printf("source: %s %d 
		}
	      }
	      else if(flag_trans){				//transient if L
		if (piknos){
		    C_trans[ b+k][k + b] -= ht->value; 
		}
		else{
		    C_trans_sparce->i[spTrCoun] = k + b;
		    C_trans_sparce->p[spTrCoun] = k + b;
		    C_trans_sparce->x[spTrCoun++] = -ht->value; 
		}
	      }
	      
	      if ((X[b+k].name = strdup( ht->name)) == NULL) {
		return NULL;
	      }
	      b++;
	    }
	    else if ( (ht->name[0] == 'i' || ht->name[0] == 'I') ){// && !flag_trans ) {
	      B[i] -= ht->value;
	     
	      if(flag_trans) {
		  sourcesRoot = insertSourse(sourcesRoot, ht->name,  i, ht->pairNode);
		  //printf("source: %s pos:%d  neg:%d B[i]:%g\n", ht->name,  i, ht->pairNode,B[i]);
	      }
	      
	      if (ht->value < 0) dcInodes.flag=1;
	      if (strcmp(dc_table->input_var,ht->name)==0)   dcInodes.node1=i;					//14/11
	      if (ht->pairNode!=-1){
		  B[ht->pairNode] += ht->value;	
		  if (strcmp(dc_table->input_var,ht->name)==0)  dcInodes.node2 = ht->pairNode;			//14/11
	      }
	    }
	    else if ( (ht->name[0] == 'C'|| ht->name[0] == 'c') && flag_trans) {				//transient if C
	      
	      if((ht->value)>=0)
		  apol = ht->value;
	      else 
		  apol = -ht->value;
	      
	      if (piknos){
		  if( ht->pairNode != -1 ){
			C_trans[i][i] += apol;
			C_trans[ht->pairNode][ht->pairNode] += apol;
			C_trans[ht->pairNode][i] -= apol;
			C_trans[i][ht->pairNode] -= apol;
		  }
		  else{
			C_trans[i][i] += apol; 
		  }     
	      }
	      else{			
		  if( ht->pairNode != -1 ){
			C_trans_sparce->i[spTrCoun] =i;
			C_trans_sparce->p[spTrCoun] =i;
			C_trans_sparce->x[spTrCoun++] = apol;
					    
			C_trans_sparce->i[spTrCoun] = ht->pairNode;
			C_trans_sparce->p[spTrCoun] = ht->pairNode;
			C_trans_sparce->x[spTrCoun++] = apol;
			
			C_trans_sparce->i[spTrCoun] = ht->pairNode;
			C_trans_sparce->p[spTrCoun] = i;
			C_trans_sparce->x[spTrCoun++] = -apol;
			
			C_trans_sparce->i[spTrCoun] = i;
			C_trans_sparce->p[spTrCoun] = ht->pairNode;
			C_trans_sparce->x[spTrCoun++] = -apol;
		  }
		  else{
			C_trans_sparce->i[spTrCoun] = i;
			C_trans_sparce->p[spTrCoun] = i;
			C_trans_sparce->x[spTrCoun++] = apol;		    
		  } 
	      }
	      
	    }	 
	    
	  }
	}
	
		/*** free head_DCNodes list ***/
		printf("/*** free some memory ***/\n");
 /* 	for (cur_tNodes = negList; cur_tNodes!=NULL; cur_tNodes=negList){
	    printf("%s -> %d\n", cur_tNodes->name, cur_tNodes->id);
	    negList = negList->next;
	    free(cur_tNodes);
	}
*/	
		/********** print hashtable *********/
	if(lines<100){
	    printf("\nid: 0\n");
	    for(i=0; i<size; i++){
	      
		if (i < valids){
		    ht = hashtable->table[ validHash[i]  ];
		}
		else{  
		    ht = remainHash->table[ i - valids];
		    if (i == valids) printf("\nremainHash:\n");
		}		

		while (ht != NULL){     
		    if (flag==1){
			printf("\nid: %d\n",i);
			flag=0;
			printf("node: %s\t nodeId: %d\t pairName %s\t pairId: %ld\t value: %10lf\t edge: %s\n", ht->key,ht->id,ht->pairNodeName, ht->pairNode, ht->value , ht->name ) ;
		    }
		    else
			printf("node: %s\t pairName %s\t pairId: %ld\t value: %10lf\t edge: %s\n", ht->key,ht->pairNodeName, ht->pairNode, ht->value , ht->name ) ;
		    ht= ht->next;
		}
		flag=1;	
	    }
	    printf("\n");
	}		

		/***  free hashtables  ***/
	printf("free hashtables\n");
	for(i=0; i<size; i++){
	    if (i < valids){
		ht = hashtable->table[ validHash[i]  ];
		for (ht = ht; ht!=NULL; ht=hashtable->table[ validHash[i] ]) {		
		    hashtable->table[ validHash[i] ] = hashtable->table[ validHash[i] ]->next;
		    free(ht);
		}
	    }
	    else{  


		ht = remainHash->table[ i - valids];

		for (ht = ht; ht!=NULL; ht=remainHash->table[i - valids]) {		
		    remainHash->table[i - valids] = remainHash->table[i - valids]->next;
		    free(ht);
		}		
		
	    }
	}
			printf("free end\n");

	free(hashtable);
	free(remainHash);
	free(validHash);

		/********** print A B X*********/
	if (nz<100) {
	  for(i=0; i<(12*(k+m2)); i++)
	    printf("-");
	}
	
	if(piknos){
	    for(i=0; i<(k+m2); i++){
	      printf("\nA[%3d]: ",i);
	      for (j=0; j<(k+m2); j++){
		if (A[i][j]>=0)
		  printf("  %5lf  ",A[i][j]);
		else
		  printf(" %5lf  ",A[i][j]);
	      }
	    }
	    if(flag_trans){
	      printf("\n");
	      for(i=0; i<(k+m2); i++){
		printf("\nC_trans[%3d]: ",i);
		for (j=0; j<(k+m2); j++){
		  if (C_trans[i][j]>=0)
		    printf("  %5lf  ",C_trans[i][j]);
		  else
		    printf(" %5lf  ",C_trans[i][j]);
	      }
	    }
	    }
	}
	else{   			//sparse table			//21/11 
	   
	  if (nz<100){
	      for(i=0; i<(k+m2); i++){
		  printf("\nA[%3d]: ",i);
		  for (j=0; j<(k+m2); j++){
		      sum=0;
		      
		      for (l=0; l < Asp->nz; l++){
			if (Asp->i[l]==i && Asp->p[l]==j)
			    sum = sum + Asp->x[l];
		      }
		      
		      if (sum>0)
			printf("  %5lf  ",sum);
		      else if ( sum<0)
			printf(" %5lf  ",sum);
		      else
			printf("    --    ");
		      
		  }
	      }
	      if (flag_trans){
		  printf("paok\n");
		    for(i=0; i<(k+m2); i++){
			printf("\nC[%3d]: ",i);
			for (j=0; j<(k+m2); j++){
			    sum=0;
			    
			    for (l=0; l < C_trans_sparce->nz; l++){
			      if (C_trans_sparce->i[l]==i && C_trans_sparce->p[l]==j)
				  sum = sum + C_trans_sparce->x[l];
			    }
			    
			    if (sum>0)
			      printf("  %5lf  ",sum);
			    else if ( sum<0)
			      printf(" %5lf  ",sum);
			    else
			      printf("    --    ");
			    
			}
		    }
	      }	      
	      
	      
	  }
	  
	}
	
	if (nz<100) {
	    printf("\n");
	    for(i=0; i<(12*(k+m2)); i++)
		printf("-");
	    printf("\n");
	}
	
	for(i=0;i<(k+m2);i++){	  
	  if (nz<100) {
	      if(i<k){
		if (i==(k+m2)/2)
		  printf("X  V(%5s)  =  B[%3d]: %lf\n",X[i].name,i,B[i]);
		else
		  printf("   V(%5s)     B[%3d]: %lf\n",X[i].name,i,B[i]);
	      }
	      else 
		printf("   I(%5s)     B[%3d]: %lf\n",X[i].name,i,B[i]);
	  }
	  if (strcmp(dc_table->input_var,X[i].name)==0)
	    thesiDC = i;
	}
	
	if (nz<100) {
	  for(i=0; i<(12*(k+m2)); i++)
	    printf("-");
	  printf("\n");
	}  
	
	if (!piknos) {
	    C = cs_compress(Asp);
	    cs_spfree(Asp);
	    cs_dupl(C);
	    if(flag_trans){
		C_trans_sparce_compressed = cs_compress(C_trans_sparce);
		cs_spfree(C_trans_sparce);
		cs_dupl(C_trans_sparce_compressed);
	    }
	} 
	
	printf("thesiDC-SwapNode %d, flag_biCG %d, flag_CG %d, flag_cholesky %d, itol %lf\n", thesiDC, flag_biCG, flag_CG, flag_cholesky, ITOL);
	
	
		/*** call solving functions ***/
	printf("/*** call solving function ***/\n");	
	if(flag_trans){
	  if(piknos){
	    printf(" h prin %g\n",trans_h_fin[0]);
	      if (trans_method == 0) {
		TR_b( B, trans_h_fin[0], trans_h_fin[1], k+m2, head_elem, sourcesRoot, A, C_trans, head_plot, X);
	      }else {
		BE_b( B, trans_h_fin[0], trans_h_fin[1], k+m2, head_elem, sourcesRoot, A, C_trans, head_plot, X);
	      }
	  }
	  else{
	      if (trans_method == 0) {
		printf("tr sparse\n");
		TR_b_sparse( B, trans_h_fin[0], trans_h_fin[1], k+m2, head_elem, sourcesRoot, C, C_trans_sparce_compressed, head_plot, X);
	      }else {
		printf("be sparse\n");
		BE_b_sparse( B, trans_h_fin[0], trans_h_fin[1], k+m2, head_elem, sourcesRoot, C, C_trans_sparce_compressed, head_plot, X);
	      }	    
	  }
	}
  	else if (piknos){
	  if (flag_biCG)
	      Bi_CG(A, B, X, ITOL, k+m2, dc_table, head_plot, thesiDC, dcInodes);
	  else if (flag_CG) 
	      CG(A, B, X, ITOL, k+m2, dc_table, head_plot, thesiDC, dcInodes);
	  else if (flag_cholesky)
	      decoCHOLESKY(A, B,k+m2, dc_table, head_plot, thesiDC, X, dcInodes);								//14/11
	  else 
	      decoLU(A,B,k+m2, dc_table, head_plot, thesiDC, X, dcInodes);								//14/11
	}
	else{
	  if (flag_biCG){
	      printf("Bi_CG\n");
	      cs_Bi_CG(C, B, X, ITOL, k+m2, dc_table, head_plot, thesiDC, dcInodes);
	  }else if (flag_CG){
	      printf("CG\n");
	      cs_CG(C, B, X, ITOL, k+m2, dc_table, head_plot, thesiDC, dcInodes);
	  }else if (flag_cholesky)
	      cs_Cholesky(C, B, k+m2, dc_table, head_plot, thesiDC, X, dcInodes);								//14/11
	  else 
	      cs_LU(C, B, k+m2, dc_table, head_plot, thesiDC, X, dcInodes);		  	  
	}
	 
	if(flag_trans){    /***  free head_elem list  ***/
	  for (curElem = head_elem; curElem!=NULL; curElem=head_elem) {		
	    head_elem = head_elem->next;
	    free(curElem);
	  }
	}
	printf("THE END\n");
	return hashtable;
}

  
  
