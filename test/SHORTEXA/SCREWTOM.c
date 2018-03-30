/*
----------------------------------------------------------------------------
SCREWTOM.C:
  Example of function  screwtom().  For detail see User's Manual.

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
	MAT4 Q;
	real fi=PIG_2;
	real h=0;
	AXIS u=Zaxis_n;
	POINT P={1., 1., 0., 1.};

	screwtom(u,fi,h,P,Q);

	printm4("The Q matrix is :",Q);
	return(0);
}
