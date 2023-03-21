#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

char TableName[] = { "MT_ORIGIN_STAGE" };

/***
char DataBaseFile[] = { "/Users/ichinose/Work/mtinv.v3.0.6/data/mt.db" };
char FLINN_ENGDAHL_PROGRAM[] = { "/Users/ichinose/Work/mtinv.v3.0.6/bin/FlinnEngdahl" };
***/
char DataBaseFile[256];
char FLINN_ENGDAHL_PROGRAM[256];

int main( int ac, char **av )
{
	FILE *fp;
	char *database;
	char table[256];
	char program[256];
	char commandline[1024];
	char *mtinv_path;

	int imax = 1;  /*** true just do the last load ***/
	int cleanup = 1; /** delete tmp files ***/

	int setpar(int,char **),getpar(),mstpar();
	void endpar();
	int verbose = 1;

/*** start program ***/

	fprintf( stderr,
"%s: %s: This program automatically adds georegion grn and grname to sqlite3 database MT solution tables. Requires SQLITE3 application https://www.sqlite.org.\n",
		__FILE__, __func__ );

/*** get default from enviromental variable, else default is from static var ***/

	if( (database = getenv( "MTINV_DATABASE_FILE" )) == NULL )
	{
	  if(verbose) 
	  {
		sprintf( DataBaseFile, "/Users/ichinose/Work/mtinv.v3.0.6/data/mt.db");
	   fprintf( stderr, "%s: %s: no MTINV_DATABASE_FILE set, using hardwired default %s\n",
		__FILE__, __func__, DataBaseFile);
	   fflush(stderr);
	  }

	  sprintf( DataBaseFile, "%s", database );

	  if(verbose)
	  {
	    fprintf( stderr, "%s: %s: MTINV_DATABASE_FILE set %s\n", __FILE__, __func__, database );
	    fflush(stderr);
	  }
	}
	else
	{
          if(verbose)
          {
	    fprintf( stderr, "%s: %s: enviromental varibale MTINV_DATABASE_FILE set=%s\n",
		__FILE__, __func__, database );
	  }
	}

	strcpy( table, TableName );

	if( ( mtinv_path = getenv( "MTINV_PATH" ) ) == NULL )
	{
	  sprintf( FLINN_ENGDAHL_PROGRAM, "/Users/ichinose/Work/mtinv.v3.0.6/bin/FlinnEngdahl" );	
		fprintf( stderr, "%s: %s: enviromental varibale MTINV_PATH not set using hardwired default %s\n",
			 __FILE__, __func__, FLINN_ENGDAHL_PROGRAM );
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

/*******************/
/*** begin query ***/
/*******************/
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
