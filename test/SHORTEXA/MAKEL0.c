/*
----------------------------------------------------------------------------
MAKEL0.C:
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

int main (void)
{
	POINT O=ORIGIN;
	real pitch=0.;
	AXIS u=Zaxis;
	MAT4 L0;

	makeL(Rev,u,pitch,O,L0);

	printm4("The L matrix in the frame 0 is:",L0);
	return(0);
}
