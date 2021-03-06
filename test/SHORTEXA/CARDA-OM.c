/*
----------------------------------------------------------------------------
CARD_OM.C:
  Example of function  cardanto_OMEGA().  For detail see User's Manual.

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
	VECTOR q={10.,5.,12.};
	VECTOR qp={0.,2.,1.};
	VECTOR omega;

	cardanto_OMEGA(q,qp,Z,Y,X,omega);

	printv("\nThe angular velocity omega is:",omega,3);
	return(0);
}
