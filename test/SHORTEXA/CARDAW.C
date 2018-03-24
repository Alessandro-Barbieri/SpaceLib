/*
----------------------------------------------------------------------------
CARDAW.C:
  Example of function  cardantoW().  For detail see User's Manual.

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
	MAT4 W;
	POINT O={50.,10.,100.,1.};
	VECTOR q={0.1,0.5,0.9};
	VECTOR qp={0.2,0.4,0.1};

	cardantoW(q,qp,X,Z,Y,O,W);

	printm4("\nThe velocity matrix W is:",W);
}