/*
----------------------------------------------------------------------------
CARDAM.C:
  Example of function  cardantoM().  For detail see User's Manual.

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
	POINT O={100.,200.,300.,1.};
	MAT4 m;
	VECTOR q={1.,2.,1.5};

	cardantoM(q,X,Z,Y,O,m);

	printm4("\nThe position matrix m is:",m);
	return(0);
}
