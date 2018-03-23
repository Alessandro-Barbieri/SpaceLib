/*
----------------------------------------------------------------------------
CARDAH.C:
  Example of function  cardantoH().  For detail see User's Manual.

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

void main(void)
{
	MAT4 H;
	POINT O={50.,10.,100.,1.};
	VECTOR q={0.1,0.5,0.9};
	VECTOR qp={0.2,0.4,0.1};
	VECTOR qpp={0.5,1.2,0.3};

	cardantoH(q,qp,qpp,X,Z,Y,O,H);

	printm4("\nThe acceleration matrix H is:",H);
}
