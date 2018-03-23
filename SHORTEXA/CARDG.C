/*
----------------------------------------------------------------------------
CARDG.C:
  Example of function  cardanto_G().  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ ing.unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#include "spacelib.h"

void main(void)
{
	MAT4 A;
	VECTOR q={0.1,0.5,0.9};
	VECTOR qp={0.2,0.4,0.1};
	VECTOR qpp={0.5,1.2,0.3};

	clear4(A);
	cardanto_G(q,qp,qpp,Y,X,Z,M A,4);

	printm4("\nThe matrix A is:",A);
}
