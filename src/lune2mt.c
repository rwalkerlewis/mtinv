/*********************************************************************************************************/
/* 
lune2mt.c - G. Ichinose 
Sat Feb 11 15:24:11 PST 2023

Creates a grid of lune_lon, lune_lat, piso, pdc, pclvd, mxx, myy, mzz, mxy, mxz, myz
	evenly spaced 1 by 1 degree intervals 
	random even distribution (mteig.c) style
	output in csv format 

This does not forward compute synthetics and variance reduction (see mteig.c) 
only a grid for testing in python classifiers

*/
/********************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>

#include "../include/mteig.h" /** global datatype and structure declarations **/

char progname[128];

int main( int ac, char **av )
{
/************************/
/*** local variables ***/
/************************/

	Results *r; /*** size is nsim_evec; this holds all results for single ***/
                    /*** createlambda() iteration randomEigVec() ***/

	Results *rbest; /*** size is nsim_eig; from all r[] this is the best ***/
			/*** fitting (highest var_red) for each createlambda() run ***/
	
	void print_results( Results *r, FILE *fp ); /*** just prints ptr to 1-row ***/
	void print_results_all( Results *r, int nrows, char *filename ); /*** prints all ***/

	Eigenvalue e[3] = { { 0, { 0, 0, 0 }, 0, 0 }, 
			    { 0, { 0, 0, 0 }, 0, 0 },
			    { 0, { 0, 0, 0 }, 0, 0 } };
/*** sorted eig ***/
	Eigenvalue e_sort[3] = { { 0, { 0, 0, 0 }, 0, 0 }, 
                                 { 0, { 0, 0, 0 }, 0, 0 },
                                 { 0, { 0, 0, 0 }, 0, 0 } };

	float Mo = 1.0E+24;
	int i, j, k;
	int nsim_eig = 2000;
	int nsim_evec= 4000;  /*** obsolete due to eigvec_fac, still used as a place holder ***/
	int eigvec_fac = 17000; /*** von Mises-Fisher PDF default nsim_evec range 1000-8000 kappa=1.3 ***/
	int sum = 0;

/*** variable nsim_evec as a function of lunelat/lon vMF3d PDF distrib ***/
	int nevec;
	int iteration = 0;
	Tensor M;
	int seed = 1;   /*** default use random seed on each run ***/
	int verbose = 0;

	char pathname[128], currentdir[128];

/*** this adds the DC (0,0) and +iso(0,90) and -iso (0,-90) DEFAULT-off ***/
	int Add_DC_iso = 1;   /*** default is to add these lune points to search ***/

	int itag, Ntags=5;   /*5: 0=DC 1=+iso 2=-iso 3=+CLVD 4=-CLVD */
	void Add_DC_iso_2eigenvals( int itag, Results *rbest, int eigvec_fac, float Mo );

/*** adds a user supplied moment tensor ***/
	int Add_user_MT = 0;
	Tensor Musr;

/*** adds a user supplied eigenvalues ***/
	int Add_user_Eig = 0;
	float e_user[3];

	void abs_sort( float *e, int verbose );

	void scalevec( float *x, float sc );

	void eig2eig( float *lamb1, float *lamb2, float *lamb3,
                float *percent_dc, float *percent_clvd, float *percent_iso, 
                float *lune_lat, float *lune_lon, float *hudson_k, float *hudson_T, 
		float *gcarc_dc, float *gcarc_exp, int *nevec, int eigvec_fac, int verbose );

/**********************************/
/*** Internal Local subroutines ***/
/**********************************/

        float pdc = 0, pclvd = 0, piso = 0;
        float e0 = 0, e1 = 0, e2 = 0;
        float lune_lat, lune_lon, hudson_k, hudson_T;
	float gcarc_dc = 0, gcarc_exp = 0;

        void createlamb(
                float *pdc,
		float *pclvd,
		float *piso,
                float *lamb1, float *lamb2, float *lamb3,
                float *lune_lat, float *lune_lon,
		float *hudson_k, float *hudson_T,
		float *gcarc_dc, float *gcarc_exp, 
                int seed,
		int *nevec,
		int eigvec_fac,
		int verbose );

/*** interpolate along the deviatoric disc ***/
	float dMo = 0.05;  /*** preset the spacing on the deviatoric disc ***/
	Results *create_clvd_eigs( Results *rbest, int *nsim_eig, int eigvec_fac, float Mo, float dMo, int verbose );

	void create_lune2mt_grid(
		int nsim_eig,
		int nsim_evec, 
		Results *rbest, 
		float Mo,
		int seed, int verbose );

/********************************/
/*** Old External subroutines ***/
/********************************/

	int mtdegfree = 6;
	
/*** shorten_path.c ***/
	char *shorten_path( char *pathname, char *filename );

/************************/
/*** misc             ***/
/************************/
        int setpar(int,char **), mstpar(), getpar();
        void endpar();

/******************************/
/*** start program ***/
/******************************/
	strcpy( progname, av[0] );
	strcpy( pathname, progname );
	shorten_path( pathname, progname );
	getcwd( currentdir, 128 );

	if( verbose )
        {
          fprintf( stdout, "%s: %s: %s: STDOUT: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname,
                __FILE__,
                __func__,
                Version_Label, Version_Date, pathname );
        }

        fprintf( stderr, "%s: %s: %s: STDERR: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname,
                __FILE__,
                __func__,
                Version_Label, Version_Date, pathname );

/*** get command line parameters ***/

	setpar(ac,av);
	getpar( "verbose",  "b", &verbose);
	getpar( "eigvec_fac", "d", &eigvec_fac );
/*** fractional moment for search points on deviatoric disc ***/
/*** override the spacing of search points on deviatoric disc ***/
	getpar( "dMo", "f", &dMo );
	if( dMo <= 0 ) dMo = 0.001;
	if( dMo >= 1.0 ) dMo = 0.5;

	getpar( "Mo", "f", &Mo );  
	mstpar( "nsim_eig", "d", &nsim_eig );
        getpar( "nsim_evec", "d", &nsim_evec );
        getpar( "seed", "d", &seed );
	getpar( "Add_DC_iso", "b", &Add_DC_iso );
	getpar( "Add_user_Eig", "b", &Add_user_Eig );
	if( Add_user_Eig )
	{ 	/*** required parameters if option selected ***/
		mstpar("e0", "f", &(e_user[0]) );
		mstpar("e1", "f", &(e_user[1]) );
		mstpar("e2", "f", &(e_user[2]) );
	}
	endpar();


/**** start processing ***/

	fprintf( stderr, "%s: %s: %s: STDERR:\n", progname, __FILE__, __func__ );

	fprintf( stderr, "%s: %s: %s: STDERR: START! nsim_eig=%d nsim_evec=%d total-iter=%d\n",
		progname, __FILE__, __func__, nsim_eig, nsim_evec, (nsim_eig*nsim_evec) );

	fprintf( stderr, "%s: %s: %s: STDERR:          currentdir=%s\n", progname, 
		__FILE__, __func__, currentdir );

	fprintf( stderr, "%s: %s: %s: STDERR:\n", progname, __FILE__, __func__ );

/* Mw   = log10( Mo )/1.5 - 10.73; */

/***********************************************************************************/
/*** set inital random seed                                                      ***/
/***********************************************************************************/

	if(seed > 0)
		srandom(seed);
	else
		srand(abs(seed));

/***********************************************************************************/
/*** creat eigenvalues first and save lon,lat,k,T,e1,e2,e3 in results structure ****/
/***********************************************************************************/

	rbest = (Results *) calloc( nsim_eig, sizeof(Results) );

	for( j = 0; j < nsim_eig; j++ )
	{
		createlamb(
			&pdc, &pclvd, &piso, &e0, &e1, &e2, 
			&lune_lat, &lune_lon, &hudson_k, &hudson_T, 
			&gcarc_dc, &gcarc_exp, seed, &nevec, eigvec_fac, verbose );
		
		rbest[j].pdc      = pdc;
		rbest[j].pclvd    = pclvd;
		rbest[j].piso     = piso;
		rbest[j].e0       = e0;
		rbest[j].e1       = e1;
		rbest[j].e2       = e2;
		rbest[j].lune_lat = lune_lat;
                rbest[j].lune_lon = lune_lon;
		rbest[j].k        = hudson_k;
		rbest[j].T        = hudson_T;
		rbest[j].Mtotal   = Mo;
		rbest[j].FullMT.eig[1].val = e0;
		rbest[j].FullMT.eig[2].val = e1;
		rbest[j].FullMT.eig[3].val = e2;
		rbest[j].gcarc_dc  = gcarc_dc;
                rbest[j].gcarc_exp = gcarc_exp;
		rbest[j].nevec     = nevec;
	}
	
/*** add the deviatoric disk ***/
/**** non-returning function create_clvd_eigs() was not reallocating memory ****/
          /* nclvd = (((Mclvd1 - (Mclvd0)) / dMo ) + 1); */
          /* rbest = realloc( rbest, ( nsim_eig + (((1 - (-1)) / dMo ) + 1) )*sizeof(Results) ); */

/***
	  fprintf( stderr, "%s: %s: %s: STDERR calling create_clvd_eigs() nsim_eig = %d\n",
		progname, __FILE__, __func__, nsim_eig );
	  rbest = create_clvd_eigs( rbest, &nsim_eig, eigvec_fac, Mo, dMo, verbose );
	  fprintf( stderr, "%s: %s: %s: STDERR done with create_clvd_eigs() nsim_eig = %d increased\n",
		progname, __FILE__, __func__, nsim_eig );
***/

/*** add the special lune verticies in source-type space ***/

	if( Add_DC_iso )
	{
		fprintf( stderr,
"%s: %s: %s: STDERR Add_DC_iso=%d allocating memory for rbest(nsim_eig=%d Ntags=%d) and calling Add_DC_iso_2eigenvals()\n",
			progname, __FILE__, __func__, Add_DC_iso, nsim_eig, Ntags );

		rbest = realloc( rbest, (nsim_eig+Ntags)*sizeof(Results) );

		itag=0;
		for( j = nsim_eig; j < nsim_eig+Ntags; j++ )
		{
			Add_DC_iso_2eigenvals( itag, &(rbest[j]), eigvec_fac, Mo );
			itag++;
		}
		nsim_eig += Ntags;

		fprintf( stderr, "%s: %s: %s: STDERR Add_DC_iso=%d done with Add_DC_iso_2eigenvals()\n",
			progname, __FILE__, __func__, Add_DC_iso );
	}

/*** add user supplied special lat/lon lune point from inversion to test ***/

	if( Add_user_Eig )
	{
		fprintf( stderr, "%s: %s: %s: STDERR Add_user_Eig=%d allocating memory for rbest and calling eig2eig()\n",
                        progname, __FILE__, __func__, Add_user_Eig );

		rbest = realloc( rbest, (nsim_eig+2)*sizeof(Results) );
		
		eig2eig( &(e_user[0]), &(e_user[1]), &(e_user[2]), &pdc, &pclvd, &piso, 
                        &lune_lat, &lune_lon, &hudson_k, &hudson_T, 
			&gcarc_dc, &gcarc_exp, &nevec, eigvec_fac, verbose );

		j = nsim_eig;
		rbest[j].pdc      = pdc;
                rbest[j].pclvd    = pclvd;
                rbest[j].piso     = piso;
                rbest[j].e0       = e_user[0];
                rbest[j].e1       = e_user[1];
                rbest[j].e2       = e_user[2];
                rbest[j].lune_lat = lune_lat;
                rbest[j].lune_lon = lune_lon;
                rbest[j].k        = hudson_k;
                rbest[j].T        = hudson_T;
                rbest[j].Mtotal   = Mo;
                rbest[j].FullMT.eig[1].val = e_user[0];
                rbest[j].FullMT.eig[2].val = e_user[1];
                rbest[j].FullMT.eig[3].val = e_user[2];
		rbest[j].gcarc_dc  = gcarc_dc;
		rbest[j].gcarc_exp = gcarc_exp;
		rbest[j].nevec     = nevec;

		nsim_eig++;
	}

/**********************************************************************************/
/*** recount the number of eigenvectors ***/
/**********************************************************************************/

	for( sum = 0, j = 0; j < nsim_eig; j++ )
		sum += rbest[j].nevec;

	fprintf( stderr,
	  "%s: %s: %s: STDERR: Resized nsim_eig= %d total eigvecs=%d writting file eigvals.out\n", 
		progname, __FILE__, __func__, nsim_eig, sum );
	
/* this only has eigenvalues */
	/* print_results_all( rbest, nsim_eig, "eigvals.out" ); */

/**********************************************************************************/
/*** loop over all eigenvalues stored in rbest[] and search random eigenvectors ***/
/**********************************************************************************/

	fprintf( stderr, "%s: %s: %s: STDERR: calling create_lune2mt_grid():\n", progname, __FILE__, __func__ );

	create_lune2mt_grid( nsim_eig, nsim_evec, rbest, Mo, seed, verbose );

/****************************/
/*** free allocated memory ***/
/****************************/

        free(rbest);
        fprintf( stderr, "%s: %s: %s: STDERR: END OF PROGRAM\n\n", progname, __FILE__, __func__ );
	exit(0);
  /****************************/
} /*** end of main() ***/
  /****************************/


