#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/nrutil.h"
#include "../include/mt.h"

char progname[128];

typedef struct {
	int ista;
	int nz;
	float *z;
} DepthVector;

int main( int ac, char **av )
{

/************************/
/*** depth vector     ***/
/************************/

	int iz;
	DepthVector *z;

/************************/
/*** event info stuff ***/
/************************/

	EventInfo *ev;
	EventInfo *glib2inv_get_input_parameters( char *, EventInfo *, int *, int );
	int ista, nsta;
	char evinfo_filename[256];

/********************/
/*** Greens stuff ***/
/********************/
	
	Greens **grn;
	FILE *fpin;

/*******************/
/*** local stuff ***/
/*******************/

	int iparallel = 0;
	int debug = 0;
	int verbose = 0;
	int idumpgrn = 0;
	int idumpsac = 0;
	int test_special = 0;
	int DIFFoperator = 3;
	char pathname[128];

/***************************/
/*** function prototypes ***/
/***************************/

	/*** glib2inv_serial.c ***/
	void glib2inv_serial( Greens **grn, EventInfo *ev, DepthVector *z, int nsta, int idumpsac, int idumpgrn, int verbose );

	/*** glib2inv_parallel.c ***/
	void glib2inv_parallel( Greens **grn, EventInfo *ev, DepthVector *z, int nsta, int idumpsac, int idumpgrn, int verbose );

	/*** wrtgrn2sac.c ***/
	/* void wrtgrn2sac( Greens *g, int ista ); */
	void wrtgrn2sac( Greens *g, int ista, char *wavetype );

	/*** shorten_path.c ***/
	char *shorten_path( char *pathname, char *filename );

	/*** load_special_grns.c ***/
	float *load_special_grns( EventInfo *ev, Greens **grn, int nsta, int *nz_tmp, int verbose );

	/*** glib2inv_serial.c ***/
	void glib2inv_special( Greens **grn, EventInfo *ev, DepthVector *z, int nsta, int idumpsac, int idumpgrn, int verbose );

	int setpar(int ac, char **av);
	int mstpar(), getpar();
	void endpar(void);
	void Usage_Print();

/***************************************************************************************************/
/**** begin program                                                                              ***/
/***************************************************************************************************/

/*****************************************************/
/*** get the input parameters foreach each station ***/
/*****************************************************/

	strcpy( pathname, av[0] );
	shorten_path( pathname, progname );

	fprintf( stderr, "%s: STDERR: Version=%s ReleaseDate=%s exec full path=%s\n",
		progname, Version_Label, Version_Date, pathname );

/*****************************************************/
/*** usage                                         ***/
/*****************************************************/

	if( ac <= 1 ) Usage_Print();

/******************************************************************************/
/*** command line arguments                                                 ***/
/******************************************************************************/

	setpar( ac, av );
	mstpar( "par",      "s", &evinfo_filename );
	getpar( "verbose",  "b", &verbose );
	getpar( "dumpgrn",  "b", &idumpgrn );
	getpar( "dumpsac",  "b", &idumpsac );
	getpar( "parallel", "b", &iparallel );
	getpar( "test_special", "b", &test_special );
	endpar();

	if( test_special ) iparallel = 0;

	if(verbose)
	{
	  fprintf( stdout, "%s: STDOUT: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname, Version_Label, Version_Date, pathname );
	}

/***********************************/
/*** load in par file parameters ***/
/***********************************/

	ev  = (EventInfo *)malloc(sizeof(EventInfo));
	ev  = (EventInfo *)
		glib2inv_get_input_parameters( evinfo_filename, ev, &nsta, verbose );

	if( verbose )
	{
		fprintf( stdout, "%s: %s: %s: STDOUT: nsta=%d\n", 
			progname, __FILE__, __func__, nsta );
	
		for( ista = 0; ista < nsta; ista++ )
		{
		  fprintf( stdout, "%s: %s: %s: ista=%03d data=%s glib=%s ginv=%s npole=%d npass=%d lf=%g ",
			progname, __FILE__, __func__, ista, ev[ista].data_filename, ev[ista].glib_filename,
			ev[ista].ginv_filename, ev[ista].npole, ev[ista].npass, ev[ista].lf );

		  fprintf( stdout, "hf=%g nt=%d dt=%g tr=%g tt=%g velordisp=%d mulfac=%g iused=%d\n",
			ev[ista].hf, ev[ista].nt, ev[ista].dt, ev[ista].tr, ev[ista].tt,
			ev[ista].grd_mo_type, ev[ista].mul_factor, ev[ista].iused );
		}

		for( ista = 0; ista < nsta; ista++ )
		{
			fprintf( stdout, "%s: \t", progname );
			WriteMyTime2STDOUT( &(ev[ista].ot) );
		}
	}

/******************************************************************************/
/*** loop over stations and just read in the Green functions                ***/
/******************************************************************************/

	if(test_special)
	{
		float *z_tmp;
		int nz = 1;
		grn = (Greens **) malloc( nsta * sizeof(Greens *) );

		z_tmp = (float *) load_special_grns( ev, grn, nsta, &nz, verbose );

		z = (DepthVector *) malloc( nsta * sizeof(DepthVector) );

		for( ista = 0; ista < nsta; ista++ )
		{
			z[ista].nz = 1;
			z[ista].ista = ista;
			z[ista].z = (float *) calloc(z[ista].nz, sizeof(float));
			z[ista].z = z_tmp;
		}	

		/**************************************/
		/*** process the GFs                ***/
		/**************************************/
		glib2inv_special( grn, ev, z, nsta, idumpsac, idumpgrn, verbose );

	} /*** end test_special option true ***/

/******************************************************************************/
/*** loop over stations and just read in the Green functions                ***/
/******************************************************************************/

	if(!test_special)
	{
	  grn = (Greens **) malloc( nsta * sizeof(Greens *) );

	  z = (DepthVector *) malloc( nsta * sizeof(DepthVector) );

	  for( ista = 0; ista < nsta; ista++ )
	  {
		if( (fpin = fopen( ev[ista].glib_filename, "rb" ) ) == NULL )
		{	
			fprintf(stderr, "%s: glib2inv.c: glib2inv(): STDERR: Fatal Error, cannot open file %s\n",
				progname, ev[ista].glib_filename );
			exit(-1);
		}

		fprintf( stderr, "%s: %s: %s: STDERR: reading file %s\n", 
			progname, __FILE__, __func__, ev[ista].glib_filename );
		fprintf( stdout, "%s: %s: %s: STDOUT: reading file %s\n",                
                        progname, __FILE__, __func__, ev[ista].glib_filename );

/******************************************************************************/
/*** get the depth range info from glib file and write it into the ginv file ***/
/******************************************************************************/
		
		fread( &(z[ista].nz), sizeof(int), 1, fpin );

		z[ista].z = (float *)calloc( z[ista].nz, sizeof(float) );
		
		fread( &(z[ista].z[0]), z[ista].nz * sizeof(float), 1, fpin );
	
/**********************************************************/
/*** loop over depth and read in the Green's functions ***/
/**********************************************************/

		grn[ista] = (Greens *)malloc( z[ista].nz * sizeof(Greens) );

		for( iz = 0; iz < z[ista].nz; iz++ )
		{
			fread( &(grn[ista][iz]), sizeof(Greens), 1, fpin );

			if(debug)
			{
				fprintf( stdout,
"ista=%d %s sta=%s net=%s stla=%g stlo=%g stel=%g evla=%g evlo=%g evdp=%g rdist=%g az=%g baz=%g t0=%g dt=%g twin=%g fmax=%g damp=%g eps=%g smin=%g rigidity=%g redv=%g ts0=%g tstart=%g tend=%g Ptakeoff=%g Pray=%g Ptime=%g Pray=%g kmax=%d nt=%d %s %s nlay=%d maxlay=%d rss=%g\n", 
				ista, 
				grn[ista][iz].filename,
				grn[ista][iz].stnm,
				grn[ista][iz].net,
				grn[ista][iz].stla,
				grn[ista][iz].stlo,
				grn[ista][iz].stel,
				grn[ista][iz].evla,
				grn[ista][iz].evlo,
				grn[ista][iz].evdp,
				grn[ista][iz].rdist,	
				grn[ista][iz].az,
				grn[ista][iz].baz,
				grn[ista][iz].t0,
				grn[ista][iz].dt,
				grn[ista][iz].twin,
				grn[ista][iz].fmax,
				grn[ista][iz].damp,
				grn[ista][iz].eps,
				grn[ista][iz].smin,
				grn[ista][iz].rigidity,
				grn[ista][iz].redv,
				grn[ista][iz].ts0,
				grn[ista][iz].tstart,
				grn[ista][iz].tend,
				grn[ista][iz].Ptakeoff,
				grn[ista][iz].Prayparameter,
				grn[ista][iz].Pttime,
				grn[ista][iz].Praybottom,
				grn[ista][iz].kmax,
				grn[ista][iz].nt,
				grn[ista][iz].v.modfile,
				grn[ista][iz].v.modpath,
				grn[ista][iz].v.nlay,
				grn[ista][iz].v.maxlay,
				grn[ista][iz].g.rss[0] );
			}

			if( verbose )
			{
			  fprintf( stdout, "%s: %s: %s: STDOUT: iz=%02d z=%3g %-8.8s rdist=%g az=%g ",
				progname, __FILE__, __func__,
				iz, 
				z[ista].z[iz], 
				grn[ista][iz].stnm, 
				grn[ista][iz].rdist, 
				grn[ista][iz].az );

			  fprintf( stdout, " evdp=%g t0=%g dt=%g nt=%d %s\n",
				grn[ista][iz].evdp,
				grn[ista][iz].t0,
				grn[ista][iz].dt,
				grn[ista][iz].nt,
				grn[ista][iz].filename );
			}

			if( ev[ista].nt > grn[ista][iz].nt )
			{
			  fprintf( stderr, "%s: %s: %s: STDERR: ERROR nt=%d of othe data greater than nt=%d ",
				progname, __FILE__, __func__, ev[ista].nt, grn[ista][iz].nt );
			  fprintf( stderr, "of the Green's functions for ista=%d sta=%s.%s\n",
				ista, ev[ista].stnm, ev[ista].net );
			  exit(-1);
			}

			if( ev[ista].dt < grn[ista][iz].dt )
			{
			  fprintf( stderr, "%s: %s: %s: STDERR: ERROR dt=%g of the data is less than dt=%g ",
				progname, __FILE__, __func__, ev[ista].dt, grn[ista][iz].dt );

			  fprintf( stderr, "of the Green's function for ista=%d sta=%s.%s\n",
				ista, ev[ista].stnm, ev[ista].net );

			  exit(-1);
			}

		} /*** loop over depth - iz ***/
	
		fclose(fpin);

	  } /*** loop over stations - ista  ***/

	/**************************************************/
	/*** do the processing                          ***/
	/**************************************************/

	  if( iparallel )
	  {
		glib2inv_parallel( grn, ev, z, nsta, idumpsac, idumpgrn, verbose );
	  }
	  else
	  {
		glib2inv_serial( grn, ev, z, nsta, idumpsac, idumpgrn, verbose );
	  }

	  fprintf( stdout, "%s: %s: %s: STDOUT: freeing memory\n", progname, __FILE__, __func__ );
	  fprintf( stderr, "%s: %s: %s: STDERR: freeing memory\n", progname, __FILE__, __func__ );

	  free(z);

	} /*** if not test_special ***/

	free(ev);
	free(grn);

	fprintf( stderr, "%s: %s: %s: STDERR: Finished Program. Bye-Bye! \n\n\n", progname, __FILE__, __func__ );
	fprintf( stdout, "%s: %s: %s: STDOUT: Finished Program. Bye-Bye! \n\n\n", progname, __FILE__, __func__ );

	exit(0);
}

void Usage_Print()
{
        fprintf(stderr,
          "\n USAGE: %s par=mtinv.par [no]verbose [no]dumpsac [no]dumpgrn [no]parallel [no]test_special\n",
          progname );

        fprintf(stderr, "\n" );
        fprintf(stderr, "\t REQUIRED PARAMETERS: \n" );
        fprintf(stderr, "\t par=mtinv.par    station parameter file\n" );
        fprintf(stderr, "\n" );

        fprintf(stderr, "\t OPTIONAL PARAMETERS: \n" );
        fprintf(stderr, "\t [no]verbose           be verbosy [boolean DEFAULT off]\n" );
        fprintf(stderr, "\t [no]dumpsac           compute 3-C synthetics from Green functions using str/dip/rak,Mo,depth \n" );
        fprintf(stderr, "\t                       from par file [boolean DEFAULT off]\n" );
        fprintf(stderr, "\t [no]dumpgrn           write out the Green functions as SAC formatted files [boolean DEFAULT off]\n" );
        fprintf(stderr, "\t [no]parallel          process Greens functions using pThreads/multithreading [boolean DEFAULT on] \n" );
        fprintf(stderr, "\t [no]test_special      read in Greens functions from SAC files in Mij(r,t,z) format [boolean DEFAULT off]\n" );
        fprintf(stderr, "\t                       test_special option requires serial option (noparallel) automatically unset\n" );
        fprintf(stderr, "\n" );
	fprintf(stderr, "\t DESCRIPTION:\n");
	fprintf(stderr, "\t Reads station, filters, models, from PAR file mtinv.par and processes Green's function libraries.\n" );         
	fprintf(stderr, "\t Outputs processed Green's function libraries for input into mtinv. For processing data see sacdata2inv.\n" );
	fprintf(stderr, "\n" );
}
