#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

/* wrapper for malloc */
void *cs_malloc (int n, size_t size)
{
    return (malloc (CS_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *cs_calloc (int n, size_t size)
{
    return (calloc (CS_MAX (n,1), size)) ;
}

/* wrapper for free */
void *cs_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *cs_realloc (void *p, int n, size_t size, int *ok)
{
    void *pnew ;
    pnew = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

/* free a sparse matrix */
cs *cs_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    cs_free (A->p) ;
    cs_free (A->i) ;
    cs_free (A->x) ;
    return (cs_free (A)) ;      /* free the cs struct and return NULL */
}

/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, int ok)
{
    cs_free (w) ;                       /* free workspace */
    cs_free (x) ;
    return (ok ? C : cs_spfree (C)) ;   /* return result if OK, else free it */
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (int m, int n, int nzmax, int values, int triplet)
{
    cs *A = cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = cs_malloc (triplet ? nzmax : n+1, sizeof (int)) ;
    A->i = cs_malloc (nzmax, sizeof (int)) ;
    A->x = values ? cs_malloc (nzmax, sizeof (double)) : NULL ;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
int cs_sprealloc (cs *A, int nzmax)
{
    int ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (CS_CSC (A)) ? (A->p [A->n]) : A->nz ;
    A->i = cs_realloc (A->i, nzmax, sizeof (int), &oki) ;
    if (CS_TRIPLET (A)) A->p = cs_realloc (A->p, nzmax, sizeof (int), &okj) ;
    if (A->x) A->x = cs_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
double cs_cumsum (int *p, int *c, int n)
{
    int i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid int overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

/* C = compressed-column form of a triplet matrix T */
cs *cs_compress (const cs *T)
{
    int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    double *Cx, *Tx ;
    cs *C ;
    if (!CS_TRIPLET (T)) return (NULL) ;                /* check inputs */
    m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
    C = cs_spalloc (m, n, nz, Tx != NULL, 0) ;          /* allocate result */
    w = cs_calloc (n, sizeof (int)) ;                   /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;           /* column counts */
    cs_cumsum (Cp, w, n) ;                              /* column pointers */
    for (k = 0 ; k < nz ; k++)
    {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }
    return (cs_done (C, w, NULL, 1)) ;      /* success; free w and return C */
}

/* remove duplicate entries from A */
int cs_dupl (cs *A)
{
    int i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
    double *Ax ;
    if (!CS_CSC (A)) return (0) ;               /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    w = cs_malloc (m, sizeof (int)) ;           /* get workspace */
    if (!w) return (0) ;                        /* out of memory */
    for (i = 0 ; i < m ; i++) w [i] = -1 ;      /* row i not yet seen */
    for (j = 0 ; j < n ; j++)
    {
        q = nz ;                                /* column j will start at q */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;                        /* A(i,j) is nonzero */
            if (w [i] >= q)
            {
                Ax [w [i]] += Ax [p] ;          /* A(i,j) is a duplicate */
            }
            else
            {
                w [i] = nz ;                    /* record where row i occurs */
                Ai [nz] = i ;                   /* keep A(i,j) */
                Ax [nz++] = Ax [p] ;
            }
        }
        Ap [j] = q ;                            /* record start of column j */
    }
    Ap [n] = nz ;                               /* finalize A */
    cs_free (w) ;                               /* free workspace */
    return (cs_sprealloc (A, 0)) ;              /* remove extra space from A */
}

/* print a sparse matrix */
int cs_print (const cs *A, int brief)
{
    int p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
	FILE *fptr;
	fptr = fopen("matrix.txt","w");
    if (!A) { fprintf (fptr,"(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    if (nz < 0)
    {
        fprintf (fptr,"%d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, Ap [n]) ;
        for (j = 0 ; j < n ; j++)
        {
            fprintf (fptr,"    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                fprintf (fptr,"      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { fprintf (fptr,"  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        fprintf (fptr,"triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            fprintf (fptr,"    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { fprintf (fptr,"  ...\n") ; return (1) ; }
        }
    }
	fclose(fptr);
    return (1) ;
}

main()
{
	int n,b,bv,j;
	cs *Qt,*Qc;
	int *A,*Qti,*Qtj;
	double idt;
	double *c,*l,*r,*g,*u,*Qtx;

	n = 6; /* number of nodes */
	b = 7; /* number of branches connected to 2 nodes */
	bv = 1; /* number of branches connected to supply (or gnd) in one end */
	idt = 1; /* inverse time step (in THz or ps^-1) */

	Qt = cs_spalloc(n,n,3*b+bv+n,1,1);
	
	A = cs_malloc(2*b+bv,sizeof(int));
	c = cs_malloc(n,sizeof(double));
	l = cs_malloc(b+bv,sizeof(double));
	r = cs_malloc(b+bv,sizeof(double));
	g = cs_malloc(b+bv,sizeof(double));
	u = cs_malloc(b+bv,sizeof(double));

	if (!Qt||!A||!c||!l||!r||!g||!u)
	{
		fprintf(stderr,"Memory allocation failure\n");
		exit(1);
	}

	/* Node-to-Branch incidence matrix in sparse form */
	/* For each branch 1...b only the 2 nodes where it is connected are indicated (SMALLEST index FIRST) */
	/* The branches b+1...b+bv are signified by 1 node and are enumerated LAST */
	A[0]=0;A[1]=1;A[2]=1;A[3]=2;A[4]=3;A[5]=4;A[6]=4;A[7]=5;A[8]=0;A[9]=3;A[10]=1;A[11]=4;A[12]=2;A[13]=5;A[14]=5;
	
	c[0]=0.0525;c[1]=0.07;c[2]=0.0525;c[3]=0.0525;c[4]=0.07;c[5]=0.0625; /* node capacitances (in pF) */
	l[0]=17.5;l[1]=17.5;l[2]=17.5;l[3]=17.5;l[4]=35.0;l[5]=35.0;l[6]=35.0;l[7]=10.0; /* branch self-inductances (in pH) */
	r[0]=17.5;r[1]=17.5;r[2]=17.5;r[3]=17.5;r[4]=35.0;r[5]=35.0;r[6]=35.0;r[7]=50.0; /* branch resistances (in Ohms) */
	g[0]=0.0;g[1]=0.0;g[2]=0.0;g[3]=0.0;g[4]=0.0;g[5]=0.0;g[6]=0.0;g[7]=0.0; /* branch conductances in parallel (in Ohms^-1) */

	/* System matrix creation */
	Qti = Qt->i; Qtj = Qt->p; Qtx = Qt->x;
	for(j=0;j<b;j++)
	{
		u[j] = 1/(r[j]+l[j]*idt)+g[j];
		Qti[3*j] = A[2*j];
		Qtj[3*j] = A[2*j];
		Qtx[3*j] = u[j];
		Qti[3*j+1] = A[2*j+1];
		Qtj[3*j+1] = A[2*j+1];
		Qtx[3*j+1] = u[j];
		Qti[3*j+2] = A[2*j];
		Qtj[3*j+2] = A[2*j+1];
		Qtx[3*j+2] = -u[j];		
	}
	for(j=0;j<bv;j++)
	{
		u[b+j] = 1/(r[b+j]+l[b+j]*idt)+g[b+j];
		Qti[3*b+j] = A[2*b+j];
		Qtj[3*b+j] = A[2*b+j];
		Qtx[3*b+j] = u[b+j];
	}
	for(j=0;j<n;j++)
	{
		Qti[3*b+bv+j] = j;
		Qtj[3*b+bv+j] = j;
		Qtx[3*b+bv+j] = c[j]*idt;
	}
	Qt->nz = 3*b+bv+n;

	Qc = cs_compress(Qt);
	cs_dupl(Qc);

	cs_print(Qc,0);
	
	cs_spfree(Qt);
	cs_spfree(Qc);
	
	cs_free(A);
	cs_free(c);
	cs_free(l);
	cs_free(r);
	cs_free(g);
	cs_free(u);
}