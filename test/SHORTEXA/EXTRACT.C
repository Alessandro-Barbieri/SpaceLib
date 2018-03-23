/*
----------------------------------------------------------------------------
EXTRACT.C:
  Example of function   extract().  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ ing.unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <conio.h>
#include "SPACELIB.H"

void main (void)
{
	MAT4 Q= { {0., 1., 0., 0.},
		  {-1., 0., 0., 2.},
		  {0., 0., 1., 0.},
		  {0., 0., 0., 1.} };
	AXIS u;
	real fi;

	extract(M Q,u,&fi,4);

	printf("\nThe angle fi is : %f\n",fi);
	printv("The axis u is :",u,3);
}
