/*
----------------------------------------------------------------------------
MAKE_P.C:
  Example of function makeL.  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ ing.unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include "spacelib.h"

void main (void)
{
	POINT O=ORIGIN;
	AXIS u=Xaxis;
	MAT4 L0;
	real pitch;

	makeL(Pri,u,pitch,O,L0);

	printm4("\nThe L matrix for the prismatic joint is:",L0);
}
