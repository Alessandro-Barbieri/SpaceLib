/*
----------------------------------------------------------------------------
WTOL_P.C:
  Example of function WtoL().  For detail see User's Manual.

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

void main (void)
{
	MAT4 W = { {0.,0.,0.,0.},     {0.,0.,0.,1.4142},
		   {0.,0.,0.,1.4142}, {0.,0.,0.,0.} };
	MAT4 L;

	WtoL(W,L);

	printm4("\nThe velocity matrix is:", W);
	printm4("\nThe L matrix for a prismatic joint is:",L);

}
