/*
----------------------------------------------------------------------------
TR_MAMT.C:
  Example of function  trasf_mamt().  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
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
	MAT4 GAMMA1={ {0.,19.,-2.5,9. }, {-19.,0.,-38.5,0.5},
		      {2.5,38.5,0.,13.5},{-9.,-0.5,-13.5,0.}};
	MAT4 m={ {0.,1.,0.,1.2}, {-1.,0.,0.,-0.5},
		 {0.,0.,1.,4. }, { 0.,0.,0.,1.  }};
	MAT4 GAMMA0;

	trasf_mamt(M GAMMA1, M m, M GAMMA0,4);

	printm4("The angular momentum matrix in frame (1) is:",GAMMA1);
	printm4("The position matrix M01 is:", m);
	printm4("The angular momentum matrix in frame (0) is:",GAMMA0);
}
