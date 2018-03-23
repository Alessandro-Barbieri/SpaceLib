/*
----------------------------------------------------------------------------
JTOJ.C:
  Example of function jtoJ.  For detail see User's Manual.

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

                giovanni.legnani @ ing.unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#include "SPACELIB.H"

void main(void)

{
    MAT4 M10;

    MAT4 J0, J1;
    real m=15.71;
    real jxx=36.633;
    real jyy=36.633;
    real jzz=7.850;
    real jxy=0., jyz=0., jxz=0.;
    real xg=0.,  yg=0.,  zg=0.;

    jtoJ(m,jxx,jyy,jzz,jxy,jyz,jxz,xg,yg,zg, J0);
    printm4("Inertia Matrix in frame (0) J0:",J0);

    xg=0.;  yg=3.;  zg=4.;     /* coordinates of C. of M. in frame 1 */
    jtoJ(m,jxx,jyy,jzz,jxy,jyz,jxz,xg,yg,zg, J1);
    printm4("Inertia Matrix in frame (1) J1:",J1);
    idmat4(M10);
    M10[X][U]=xg; M10[Y][U]=yg;  M10[Z][U]=zg;
    trasf_mamt4(J0,M10,J1);
    printm4("Inertia Matrix in frame (1) J1:",J1);
}
