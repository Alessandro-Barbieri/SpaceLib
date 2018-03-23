/*
----------------------------------------------------------------------------
MAKEL.C:
  Example of function makeL.  For detail see User's Manual.

		Università degli Studi di Brescia
		Dipartimento di Ingegneria Meccanica ed Industriale
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include "SPACELIB.H"

void main (void)
{
   POINT P={0.,1.2,1.2,1.};
   real pitch=0.;
   AXIS u=Xaxis_n;
   MAT4 Lk;

	makeL(Rev,u,pitch,P,Lk);

	printm4("The L matrix in the frame k is:",Lk);
}
