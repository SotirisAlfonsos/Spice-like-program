/*		Circuit Simulation 2014 3 Solution of MNA System
 * 
 * Onoma programmatisti:	Alfonsos Swthrhs 959
 * 				Servou Katerina 852
 * 				Stavropoulos Antonis 1141
 * 				Vakerlis Christos 982
 * 
 * Hmeromhnia:		01/11/2014

 *
 * ayto to arxeio xrhsimopoieitai ws vivliothiki apo to arxeio cs_parser.c 
 * kai periexei tis domes struct kai ta tis katholikes metavlites pou xrisimopoiei to programma
 */

  #define WORDSIZE 30
  #define MAX_SIZE_GRAMMIS 500

  
  
   typedef struct numbers{									//12/12 transient
	double number;
	struct numbers *next;
  }tnumbers;
  
  
   typedef struct trans_source{								//12/12 transient
	char *method;
	tnumbers *num;
  }tTrans_source;
  
  
  
  //domes pou tha sxhmatisoun aples listes apo elements 
  
  typedef struct elements {	//includes L,C,V,I,R  circuit's elements 
	 char type;			//<type>
	 char name[WORDSIZE];		//<name>
	 char posNode[WORDSIZE];	//<+>
	 char negNode[WORDSIZE];	//<->
	 double value;			//<value>
	 tTrans_source *trans;
	 struct elements * next;	//deiktis pou deixnei ston epomeno komvo tis listas
  }tElements;

  typedef struct diodes {	//includes circuit's diodes 
	 char name[WORDSIZE];		//<name>
	 char posNode[WORDSIZE];	//<+>
	 char negNode[WORDSIZE];	//<->
	 char modName[WORDSIZE];	//<model-name>	
	 double area;			//<area>
	 struct diodes * next;		//deiktis pou deixnei ston epomeno komvo tis listas
  }tDiodes;
    
  typedef struct MOS {		//includes circuit's Mos Transistors
	 char name[WORDSIZE];		//<name>
	 char D[WORDSIZE];		//<D>
	 char G[WORDSIZE];		//<G>
	 char S[WORDSIZE];		//<S>
	 char B[WORDSIZE];		//<B>
	 char modName[WORDSIZE];	//<model-name>	
	 double L;			//<L>	
	 double W;			//<W>	
	 struct MOS * next;		//deiktis pou deixnei ston epomeno komvo tis listas
  }tMOS;
  
  typedef struct BJT {		//includes circuit's Bjt Transistors
	 char name[WORDSIZE];		//<name>
	 char C[WORDSIZE];		//<C>
	 char B[WORDSIZE];		//<B>
	 char E[WORDSIZE];		//<E>
	 char modName[WORDSIZE];	//<model-name>	
	 double area;			//<area>
	 struct BJT * next;		//deiktis pou deixnei ston epomeno komvo tis listas
  }tBJT;
  
  typedef struct entry_s {
	char *key;
	long int pairNode;
	char *name;
	int id;		//mporei kai me *id => malloc stin dilosi timis
	double value;
	char *pairNodeName;
	struct entry_s *next;
  }entry_t;
  
  typedef struct hashtable_s {
	int size;
	struct entry_s **table;	
  }hashtable_t;
  
  typedef struct nodes{
	 char *name;
	 long int id;
	 struct nodes *next;
	 struct  nodes **table;	
  }tNodes;
 
   typedef struct arrayX{
	 char *name;
//	 double value;
  }tarrayX;
  
  typedef struct plot{
	char *node;
	char type;
	struct plot *next;
  }tPlot;
	
  typedef struct dc{
	char *input_var;
	double start_value;
	double end_value;
	double increment;
  }dc_analysis;
  
  typedef struct dcI{								//14/11
	int node1;
	int node2;
	int flag;
  }dcIt;
  
  typedef struct Sources{								//13/12
	char *elemName;
	int thesiB;
	int thesiB_pairNode;
	struct Sources *next;
  }tSources;
  
  


  
	
	