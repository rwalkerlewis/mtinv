#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <sys/wait.h>

typedef struct {
	float dt;
	char modfilename[64];
	char staname[8];
	char network[8];
	char location[8];
	pid_t pid;
} PARFILE;

#define MAXTHREADS 60  /*** this is the maximum number of threads (stations/rows) allowed to be launched at once ***/

char progname[128];

int main( int ac, char *av[] )
{
	int NumPars = 5; /* par= stnm= net= loc= dt= */
	char args[][64] = { "/Users/ichinose/Work/mtinv.v3.0.3/bin/mkgrnlib", "par=", "stnm=", "net=", "loc=", "dt=", "\0" };
	char executable_pathname[256];

	/* char args[128][128]; */

	PARFILE pf[MAXTHREADS];
	int i, j;
	int NumSta;
	int retv, wait_status, wait_pid, last_pid;
	char parfilename[128];
	char executable_path[256];
	
	PARFILE *read_parfile( char *, PARFILE *, int * );

	int setpar(int, char **), mstpar(), getpar();
	void endpar();
	int verbose = 0;

	void  Print_Usage(void);

/**** main ****/

	strcpy( progname, av[0] );

	if( ac <= 1 ) Print_Usage();

	setpar(ac,av);
	mstpar("parfile", "s", parfilename);
	mstpar("executable_pathname", "s", executable_pathname );
	getpar("verbose", "b", &verbose );
	endpar();

	fprintf(stderr, "%s: STDERR: start: parfile=%s executable_pathname=%s\n",
		progname, parfilename, executable_pathname );
	if(verbose)
	{
	  fprintf( stdout, "%s: STDOUT executable_pathname=%s parfilename=%s\n", 
		progname, executable_pathname, parfilename );
	  fflush(stdout);
	}

/*******************************************************/
/*** do a safety check for file existance            ***/
/*** and avoid some problems with runaway fork bombs ***/
/*******************************************************/

	if( access( executable_pathname, F_OK|X_OK ) != 0 ) 
	{
		fprintf( stderr, "%s: STDERR: cannot find executable_pathname=%s. Ensure it's set correctly\n",
			 progname, executable_pathname );
		fprintf( stdout, "%s: STDOUT: cannot find executable_pathname=%s. Ensure it's set correctly \n", 
                         progname, executable_pathname );
		exit(-1);
	}

	if( access( parfilename, F_OK ) != 0 ) 
	{
                fprintf( stderr, "%s: STDERR: cannot find parfile=%s. Ensure it's set correctly\n",
                         progname, parfilename );
                fprintf( stdout, "%s: STDOUT: cannot find parfile=%s. Ensure it's set correctly \n", 
                         progname, parfilename );
                exit(-1);
        }


/*** do a safety check for file existance ***/

	read_parfile( parfilename, pf, &NumSta );

	if( NumSta > MAXTHREADS )
	{
	  fprintf( stderr, "%s: STDERR: NumSta = %d MAXTHREADS %d exceeded. EXIT!\n", progname, NumSta, MAXTHREADS );
	  fprintf( stdout, "%s: STDOUT: NumSta = %d MAXTHREADS %d exceeded. EXIT!\n", progname, NumSta, MAXTHREADS );
	  fflush(stderr);
	  fflush(stdout);
	  exit(-1);
	}

	if(verbose)
	{
		fprintf( stdout, "nrows read NumSta=%d\n", NumSta );
		fflush(stdout);
	}

	if(verbose)
	{
	  for( i = 0; i < NumSta; i++ )
	  {
		strcpy(  args[0], executable_pathname );
                sprintf( args[1], "par=%s",  pf[i].modfilename );
                sprintf( args[2], "stnm=%s", pf[i].staname );
                sprintf( args[3], "net=%s",  pf[i].network );
		sprintf( args[4], "loc=%s",  pf[i].location );
                sprintf( args[5], "dt=%g",   pf[i].dt );

		for( j = 0; j <= NumPars; j++ )
			fprintf( stdout, "(%d)%s ", j, args[j] );

		fprintf( stdout, "\n" );
                fflush(stdout);
	  }
	}
	sleep(3);

/***********************************************************************************/
	for( i = 0; i < NumSta; i++ )
	{
		pf[i].pid = fork();

		if( pf[i].pid < 0 )
		{
			fprintf(stdout, "fork error");
		}
		else if( pf[i].pid == 0 )
		{
			strcpy( args[0], executable_pathname );
			sprintf( args[1], "par=%s",  pf[i].modfilename );
			sprintf( args[2], "stnm=%s", pf[i].staname );
			sprintf( args[3], "net=%s",  pf[i].network );
			sprintf( args[4], "loc=%s",  pf[i].location );
			sprintf( args[5], "dt=%g",   pf[i].dt );

			retv = execlp( args[0], args[0], args[1], args[2], args[3], args[4], args[5], (char *)0 );

		} /* if fork pid */

	} /* loop over NumSta */

	while( (wait_pid = waitpid( -1, NULL, 0 ))  ) 
	{
		if( errno == ECHILD )
		{
			break;
		}
	}
	if(verbose)fprintf( stdout, "%s: STDOUT finished.\n", progname );
	fprintf(stderr, "%s: STDERR: finished.\n\n", progname );
}

PARFILE *read_parfile( char *parfilename, PARFILE *pf, int *nrows )
{
	FILE *fp;
	int i = 0;
	char rec[256];
	extern char progname[128];

	if( (fp = fopen(parfilename, "r")) == NULL )
	{
		fprintf( stderr, "%s: cannot open parfilename=%s\n", 
			progname, parfilename );
		exit(-1);
	}
	
	while( fgets(rec,256,fp) != NULL )
	{
		/**** skip comment ****/
		if( rec[0] == '#' ) continue;

		sscanf(rec,"%s %s %s %s %f",
			pf[i].staname,
			pf[i].network,
			pf[i].location, 
			pf[i].modfilename,
			&(pf[i].dt) 
		);
		i++;
	}
	fclose(fp);
	*nrows = i;
	return (PARFILE *)pf;
}

void Print_Usage( void )
{
	fprintf( stderr, "\n" );
	fprintf( stderr, "Usage: %s parfile=(string) executable_pathname=(string) [no]verbose\n", progname );
	fprintf( stderr, "\n" );
	fprintf( stderr, "REQUIRED:\n" );
	fprintf( stderr, "\t parfile=(string) this is the mkgrnlib.par file\n" );
	fprintf( stderr, "\t executable_pathname=(string) this is the full path to mkgrnlib\n" );
	fprintf( stderr, "\n" );
        fprintf( stderr, "OPTIONAL:\n" );
	fprintf( stderr, "\t [no]verbose\n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "FORMAT mkgrnlib.par: \n" );
	fprintf( stderr, "### station-code network-code location-code mkgrnlib.parfile dt(sec/sample) ### comments\n" );
	fprintf( stderr, "ANMO             IU          \"00\"          wus.par          0.1            ### comments\n" );
	fprintf( stderr, "XAN              IC          \"\"           asia.par          0.1            ### comments\n" );
	fprintf( stderr, "\n" );
        fprintf( stderr, "see mkgrnlib man page for FORMAT mkgrnlib.parfile \n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "DESCRIPTION: wrapper routine for multithreaded mkgrnlib using fork process, limits to 50 stations\n" );
	
}


