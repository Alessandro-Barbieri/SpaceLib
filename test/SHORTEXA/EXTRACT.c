/*
----------------------------------------------------------------------------
EXTRACT.C:
  Example of function   extract().  For detail see User's Manual.

		Università degli Studi di Brescia
		Dipartimento di Ingegneria Meccanica ed Industriale
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include "SPACELIB.H"

int main (void)
{
	MAT4 Q= { {0., 1., 0., 0.},
		  {-1., 0., 0., 2.},
		  {0., 0., 1., 0.},
		  {0., 0., 0., 1.} };
	AXIS u;
	real phi;

	extract(M Q,u,&phi,4);

	printf("\nThe angle phi is : %f\n",phi);
	printv("The axis u is :",u,3);
	return(0);
}
