/*
----------------------------------------------------------------------------
ROTAT.C:
  Example of function  rotat().  For detail see User's Manual.

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
	MAT4 A;
	AXIS u=Zaxis_n;
	real phi=PIG_2;

	clear4(A);
	rotat(u,phi,M A,4);

	printm4("The A matrix is :",A);
	return(0);
}
