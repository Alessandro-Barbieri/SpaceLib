/*
----------------------------------------------------------------------------
TRASF_MH.C:
  Example of function  trasf_mami().  For detail see User's Manual.

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
  MAT4 H1={ {-2.25,-0.9,  0., 0.},
	    { 0.9, -2.25, 0., 0.},
	    { 0.,   0.,   0., 0.},
	    { 0.,   0.,   0., 0.}  };

  MAT4 m01={{ 1., 0., 0., 0.4},
	    { 0., 1., 0., 0.1},
	    { 0., 0., 1., 0. },
	    { 0., 0., 0., 1. }  };

  MAT4 H0;

	trasf_mami( H1, m01, H0);

	printm4("The acceleration matrix H in frame (1) is:",H1);
	printm4("The position matrix M01 is:", m01);
	printm4("The acceleration matrix H in frame (0) is:",H0);
}
