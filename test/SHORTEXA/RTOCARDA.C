/*
----------------------------------------------------------------------------
RTOCARDA.C:
  Example of function rtocardan().  For detail see User's Manual.

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
	MAT3 R={ { 0.840, -0.395, -0.371},
		 {-0.415, -0.029, -0.909},
		 { 0.348,  0.918, -0.189}  };
	VECTOR q1,q2;

	rtocardan(M R,3,Y,X,Z,q1,q2);

	printv("\nThe first Cardan angule q1 is:",q1,3);
	printv("\nThe second Cardan angule q2 is:",q2,3);
}
