#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char progname[128];

void main( int ac, char **av )
{
	float dmo = 0;
	void create_clvd_eigs( float dmo );

	int setpar(), mstpar();
	void endpar();

	strcpy( progname, av[0] );

	setpar(ac,av);
	mstpar( "dmo", "f", &dmo );
	endpar();

	create_clvd_eigs( dmo );
}
