/*
----------------------------------------------------------------------------
CARDA_OM.C:
  Example of function  cardanto_omega().  For detail see User's Manual.

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
	MAT4 A;
	POINT O={200.,40.,300.,1.};
	VECTOR q={10.,5.,12.};
	VECTOR qp={0.,2.,1.};
	clear4(A);
	cardanto_omega(q,qp,Z,Y,X,M A,4);

	printm4("\nThe angular velocity matrix A is:",A);
}
