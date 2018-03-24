/*
----------------------------------------------------------------------------
INTERSLP.C:
  Example of function interslpl().  For detail see User's Manual.

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
PLANE pl={0.,0.,1.,-5.};
VECTOR dir={0.,0.,1.},
       v;
POINT P={0.,6.,10.,1.},
      P1;
POINT Ps;
LINE l;
int type;
real d;

	pcopy (P,l.P);
	vcopy(dir,l.dir);
	interslpl(l,pl,P1,&type);
	printv("Intersection Point P1",P1,4);
	d=dist(P,P1);
	vector(dir,d,v);
	subv(P1,v,Ps);
	Ps[U]=1.;
	printv("The symmetric point Ps is:",Ps,4);
}
