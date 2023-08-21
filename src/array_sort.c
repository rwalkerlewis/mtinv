/***
Tue Aug  4 15:25:05 PDT 2020
	example C program to do unqiue sorting of float array 
***/

#include <stdio.h>
#include <stdlib.h>

int main()
{
	float array[] = { 4, 6, 5, 3, 4, 5, 2, 8, 7, 0 };
	int i, nbytes, memory, npts;
	void sort_unique_array( float *a, int *n );

	nbytes = sizeof(array);
	memory = sizeof(array[0]);

	npts = nbytes/memory;

	fprintf( stderr, "%s: %s: nbytes=%d memory=%d npts=%d\n",
		__FILE__, __func__, nbytes, memory, npts );

	for( i = 0; i < npts; i++ )
		fprintf( stdout, "%s: %s: i=%d array=%g\n",
                        __FILE__, __func__, i, array[i] );

	sort_unique_array( array, &npts );

	for( i = 0; i < npts; i++ )
	{
		fprintf( stdout, "%s: %s: i=%d array=%g\n",
			__FILE__, __func__, i, array[i] );
	}
}
	
void sort_unique_array( float *array, int *size )
{
	int i, j, n, count = 0;
	float tmp, array1[*size];

	n = *size;

	for( i = 0; i < n; i++ )
	{
		for( j = i+1; j < n; j++ )
		{
			if( array[i] == array[j] )
			{
				/***duplicate array element found***/
				break; 
			}
		}

		/*** if j is equal to n it means we traversed whole
			array and didnt find a duplicate in array ***/
		if( j == n )
		{
			array1[count++] = array[i];
		}
	}

	/*** sort ***/
	for( i = 0; i < count-1; i++ )
	{
		for( j = i+1; j < count; j++ )
		{
			if( array1[i] > array1[j] )
			{
				tmp = array1[i];
				array1[i] = array1[j];
				array1[j] = tmp;
			}
		}
	}

	*size = count;
	fprintf( stdout, "%s: %s: new size count = %d n = %d\n", 
		__FILE__, __func__, count, n );

	for( i = 0; i < n; i++ ) array[i] = 0; /***reinitialize input array for output ***/

	/*** transfer array1 to array ***/
	for( i = 0; i < count; ++i )
	{
		array[i] = array1[i];
		fprintf( stdout, "%s: %s: i=%d array=%g\n", 
			__FILE__, __func__, i, array[i] );
	}
	return;
}
