/*
----------------------------------------------------------------------------
FRAMEV.C:
  Example of function frameV().  For detail see User's Manual.

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
	MAT4 m01;
	VECTOR v1={0.,2.,1.};
	VECTOR v2={0.,1.,2.};

	clear4(m01);
	frameV(v1,v2,X,Y,M m01,4);

	printm4("\nThe position matrix m01 is:",m01);
}
