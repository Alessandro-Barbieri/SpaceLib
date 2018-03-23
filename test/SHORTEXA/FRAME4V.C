/*
----------------------------------------------------------------------------
FRAME4V.C:
  Example of function frame4V().  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ ing.unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#include "SPACELIB.H"

void main(void)
{
	MAT4 m01;
	POINT P1={5.,4.,3.,1.};
	VECTOR v1={0.,2.,1.};
	VECTOR v2={0.,1.,2.};

	frame4V(P1, v1,v2,X,Y,M m01);

	printm4("\nThe position matrix m01 is:",m01);
}
