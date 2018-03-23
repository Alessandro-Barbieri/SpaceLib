/*
----------------------------------------------------------------------------
VELWH2.C:
  Example of function velacctoWH2.  For detail see User's Manual.

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
	MAT4 W1,H1;
	real qp=1.5;
	real qpp=0.9;

	velacctoWH2(Rev,Z,qp,qpp,W1,H1);

	printm4("\nThe velocity matrix W in frame (1) is:",W1);
	printm4("\nThe acceleration matrix H in frame (1) is:",H1);
}
