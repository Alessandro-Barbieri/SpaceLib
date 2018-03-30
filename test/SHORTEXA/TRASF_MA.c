/*
----------------------------------------------------------------------------
TRASF_MA.C:
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

int main(void)
{
  MAT4 W1={ { 0. , -1.5, 0., 0.},
	    { 1.5,  0. , 0., 0.},
	    { 0.,   0.,  0., 0.},
	    { 0.,   0.,  0., 0.} };

  MAT4 m01={{ 1., 0., 0., 0.4},
	    { 0., 1., 0., 0.1},
	    { 0., 0., 1., 0. },
	    { 0., 0., 0., 1. }  };
  MAT4 W0;

	trasf_mami( W1, m01, W0);

	printm4("The velocity matrix W in frame (1) is:",W1);
	printm4("The position matrix M01 is:", m01);
	printm4("The velocity matrix W in frame (0) is:",W0);
	return(0);
}
