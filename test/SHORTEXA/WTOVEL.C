/*
----------------------------------------------------------------------------
WTOVEL.C:
  Example of function Wtovel().  For detail see User's Manual.

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
	MAT4 W={{0.,-2.,2.5,2.5},  {2.,0.,-4.5,1.7},
		{-2.5,4.5,0.,3.2}, {0.,0.,0.,0.}  };
	AXIS u;
	real omega;
	real vel;
	POINT P;

	Wtovel(W,u,&omega,&vel,P);

	printm4("The Velocity matrix is :", W);
	printv("The axis of rotation <u> is:",u,3);
	printf("\nThe angular speed omega is: %f",omega);
	printf("\nThe linear velocity vel is: %f",vel);
	printv("\nThe point P of the axis is:",P,4);
}
