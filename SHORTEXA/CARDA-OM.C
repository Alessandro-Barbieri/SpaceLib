/*
----------------------------------------------------------------------------
CARD_OM.C:
  Example of function  cardanto_OMEGA().  For detail see User's Manual.

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
	VECTOR q={10.,5.,12.};
	VECTOR qp={0.,2.,1.};
	VECTOR omega;

	cardanto_OMEGA(q,qp,Z,Y,X,omega);

	printv("\nThe angular velocity omega is:",omega,3);
}
