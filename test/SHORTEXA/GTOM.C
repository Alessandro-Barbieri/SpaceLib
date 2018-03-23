/*
----------------------------------------------------------------------------
GTOM.C:
  Example of function gtom().  For detail see User's Manual.

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
	MAT4 Hg;
	real gx=0.;
	real gy=0.;
	real gz=-9.81;

	gtom(gx,gy,gz,Hg);

	printm4("\nThe acceleration matrix Hg is:",Hg);
}
