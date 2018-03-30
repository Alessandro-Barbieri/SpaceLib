/*
----------------------------------------------------------------------------
CARDOPTO.C:
  Example of function   cardanto_OMEGAPTO().  For detail see User's Manual.

		Università degli Studi di Brescia
		Dipartimento di Ingegneria Meccanica ed Industriale
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#include "SPACELIB.H"

int main(void)
{
	VECTOR q={1.,1.2,3.};
	VECTOR qp={0.,1.,0.};
	VECTOR qpp={3.,2.5,4.01};
	VECTOR omegapto;

	cardanto_OMEGAPTO(q,qp,qpp,Y,X,Z,omegapto);

	printv("\nThe angular acceleration omegapto is:",omegapto,3);
	return(0);
}
