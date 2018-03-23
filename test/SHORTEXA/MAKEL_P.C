/*
----------------------------------------------------------------------------
MAKE_P.C:
  Example of function makeL.  For detail see User's Manual.

		Università degli Studi di Brescia
		Dipartimento di Ingegneria Meccanica ed Industriale
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include "SPACELIB.H"

void main (void)
{
	POINT O=ORIGIN;
	AXIS u=Xaxis;
	MAT4 L0;
	real pitch;

	makeL(Pri,u,pitch,O,L0);

	printm4("\nThe L matrix for the prismatic joint is:",L0);
}
