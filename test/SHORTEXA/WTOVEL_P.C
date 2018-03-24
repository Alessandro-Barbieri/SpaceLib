/*
----------------------------------------------------------------------------
WTOVEL_P.C:
  Example of function Wtovel().  For detail see User's Manual.

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

int main(void)
{
	MAT4 W={ {0.,0.,0.,2.5},{0.,0.,0.,1.7},
		 {0.,0.,0.,3.2},{0.,0.,0.,0.}};
	AXIS u;
	real omega;
	real vel;
	POINT P;

	Wtovel(W,u,&omega,&vel,P);

	printm4("The Velocity matrix is :", W);
	printv("The axis u is:",u,3);
	printf("The angular speed omega is: %f",omega);
	printf("The linear velocity vel is: %f",vel);
	printv("The point P of the axis is:",P,4);
}
