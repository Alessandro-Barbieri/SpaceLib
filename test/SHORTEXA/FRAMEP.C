/*
----------------------------------------------------------------------------
FRAMEP.C:
  Example of function frameP().  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#include "SPACELIB.H"

void main(void)
{
	MAT3 R01;
	MAT4 m01;
	POINT P1={5.,4.,3.,1.};
	POINT P2={5.,6.,4.,1.};
	POINT P3={5.,5.,5.,1.};

	clear4(m01);
	frameP(P1,P2,P3,X,Y,M m01,4);
	frameP(P1,P2,P3,X,Y,M R01,3);

	printm3("\nThe rotation matrix R01 is:",R01);
	printm4("\nThe position matrix m01 is:",m01);
}
