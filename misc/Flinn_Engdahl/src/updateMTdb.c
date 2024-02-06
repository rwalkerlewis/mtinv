#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

char TableName[] = { "MT_ORIGIN_STAGE" };
char DataBaseFile[256];  /*** getenv MT_DATABASE_FILE ***/
char FLINN_ENGDAHL_PROGRAM[256]; /*** getenv MTINV_PATH + "/bin/FlinnEngdahl" ***/

int main( int ac, char **av )
{
	FILE *fp;
	char *database;
	char table[256];
	char program[256];
	char commandline[1024];
	char *mtinv_path;

	int imax = 1;  /*** true just do the last load, lookup the last orid from MT_ORIGIN_STAGE table ***/
	int cleanup = 1; /** delete all query and tmp output files ***/

	int setpar(int,char **);
	int mstpar( char *, char *, void * );
	int getpar( char *, char *, void * );
	void endpar();
	int verbose = 1;

/*** start program ***/

	fprintf( stderr,
"%s: %s: This program automatically adds georegion grn and grname to sqlite3 database MT solution tables. Requires SQLITE3 application https://www.sqlite.org.\n",
		__FILE__, __func__ );

/*** get default from enviromental variable, else default is from static var ***/
		/** MT_DATABASE_FILE not MTINV_DATABASE_FILE **/

	if( (database = getenv( "MT_DATABASE_FILE" )) == NULL )
	{
	  if(verbose) 
	  {
	   fprintf( stderr, "%s: %s: ERROR: MT_DATABASE_FILE not set! ", __FILE__, __func__ );
	   fprintf( stderr, "Use setenv MT_DATABASE_FILE /home/user/mtinv.version/data/mt.db\n" );
	   fflush(stderr);
	  }
	  exit(-1);
	}
	else
	{
	  sprintf( DataBaseFile, "%s", database );
	  if(verbose)
	  {
	    fprintf( stderr, "%s: %s: MTINV_DATABASE_FILE set %s\n", __FILE__, __func__, database );
	    fflush(stderr);
	  }
	}

	strcpy( table, TableName );

	if( ( mtinv_path = getenv( "MTINV_PATH" ) ) == NULL )
	{
	  fprintf( stderr, "%s: %s: ERROR: enviromental varibale MTINV_PATH not set!\n",
		__FILE__, __func__ );
	  fprintf( stderr, "Use setenv MTINV_PATH /home/user/mtinv.version\n" );
	  exit(-1);
	}
	else
	{
	  sprintf( FLINN_ENGDAHL_PROGRAM, "%s/bin/FlinnEngdahl", mtinv_path );
          if(verbose)
   	  {
	    fprintf( stderr, "%s: %s: enviromental varibale MTINV_PATH set=%s FLINN_ENGDAHL_PROGRAM=%s\n",
		__FILE__, __func__, mtinv_path, FLINN_ENGDAHL_PROGRAM );
	  }
	}
	strcpy( program, FLINN_ENGDAHL_PROGRAM );

	setpar(ac,av);
        getpar("db","s", database );
        getpar("tb","s", table );
        getpar("prog","s", program );
        getpar("max","b", &imax );
        getpar("clean","b", &cleanup );
        endpar();

        fprintf( stderr, "%s: %s: database=%s table=%s program=%s imax=%d(1=just do last load) clean=%d(1=delete tmp files)\n",
                __FILE__, __func__, database, table, program, imax, cleanup );

/***************************************************************/
/*** begin query                                             ***/
/*** for the lat and lon from the last MT_ORIGIN_STAGE.orid  ***/
/***************************************************************/
	fp = fopen( "query.csh","w");
	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "sqlite3 %s << EOF\n", database );
	fprintf( fp, ".output dump.out\n" );
	fprintf( fp, ".mode column\n" );
	fprintf( fp, ".header off\n" );

	if(imax)
	  fprintf( fp, "select lon, lat, max(orid) from %s;\n", table );
	else
	  fprintf( fp, "select lon, lat, orid from %s;\n", table );

	fprintf( fp, ".quit\n" );
	fprintf( fp, "EOF\n" );
	fclose(fp);
	system( "/bin/csh ./query.csh;" );

	sprintf( commandline, "%s < dump.out > update.sql", program );
	system( commandline );
	/* system( "/Users/ichinose/bin/FlinnEngdahl < dump.out > update.sql" ); */


/*****************************************************************************/
/*** for the MT_ORIGIN_STAGE.orid update the MT_ORIGIN_STAGE.grn GregionID ***/
/*****************************************************************************/
	fp = fopen( "update.csh", "w" );
	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "sqlite3 %s << EOF\n", database );
	fprintf( fp, ".read update.sql\n" );
	fprintf( fp, ".quit\n" );
	fprintf( fp, "EOF\n" );
	fclose(fp);

	system( "/bin/csh ./update.csh;" );
	if(cleanup)
	{
	  system( "/bin/rm -f dump.out update.csh update.sql query.csh" );
	}
}
