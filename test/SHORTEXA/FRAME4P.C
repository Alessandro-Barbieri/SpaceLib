/*
----------------------------------------------------------------------------
FRAME4P.C:
  Example of function frame4P().  For detail see User's Manual.

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
	MAT4 m01,ma1, m1a, m0a;
	POINT Ptmp;

	POINT P1={ 100.0,40.0,25.0,1.};
	POINT P2={170.0,270.,105.,1.};
	POINT P3={220.0,0.0,98.0,1.};

	POINT Q1={496.70941281,-72.36019493,72.52357206,1.};
	POINT Q2={580.82911538,152.34803995,155.02773018,1.};
	POINT Q3={407.44531117,111.38489437,144.51470414,1.};

	frame4P(P1,P2,P3,X,Y,ma1);
	frame4P(Q1,Q2,Q3,X,Y,m01);
	invers( ma1, m1a );
	molt4(m01, m1a , m0a);

	printm4("\nThe position matrix m01 is:",m01);
	printm4("\nThe position matrix ma1 is:",ma1);
	printm4("\nThe position matrix m0a is:",m0a);
	getch();
	moltp( m0a, P1, Ptmp );
	printv("temp", Ptmp,4);
	printv("Q1", Q1,4);

	moltp( m0a, P2, Ptmp );
	printv("temp", Ptmp,4);
	printv("Q1", Q2,4);

	moltp( m0a, P2, Ptmp );
	printv("temp", Ptmp,4);
	printv("Q1", Q2,4);

	getch();
}
