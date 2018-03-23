/*
----------------------------------------------------------------------------
PROJPONL.C:
  Example of function  projpon().  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ ing.unibs.it
-----------------------------------------------------------------------------
*/   
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include "SPACELIB.H"

void main(void)
{
   LINE l={ {3.,7.2,2.05,1.}, {0.7,4.,9.} };
   LINE m;
   POINT P={5.,1.,3.,1.};
   POINT P1;
	projponl(P,l,P1);
	line2p(P,P1,&m);
	printv("2nd point [projection of P] is",P1,4);
	printv("The origin of the new line is:",m.P,4);
	printv("The direction of the new line is:",m.dir,3);
}
