#include <stdio.h>
#include "SPACELIB.H"
#include "LINEAR.H"


MAT4 J={{21., -2., -2., 20.0},
{-2., 1., 0., 0.},
{-2., 0., 1., 0.},
{20., 0., 0.,20.}};

var[2][6] ={{1,1,1,1,1,1},{0,0,0,0,0,0}};
var2[2][6]={{0,0,0,0,0,0},{1,1,1,1,1,1}};
var3[2][6]={{0,1,0,1,0,1},{1,0,1,0,1,0}};
var4[2][6]={{1,0,1,1,0,1},{0,1,0,0,1,0}};
var5[2][6]={{0,1,0,0,1,0},{1,0,1,1,0,1}};

int main()
{
	MAT4 Wp,F,F2;

	real fx=0., fy=8.334, fz=0, cx=0.8334, cy=0., cz=8.35;
	real jxx,jyy,jzz, jxy,jxz,jyz, xg, yg, zg, m;
	
	printf("----------------\n");


	actom(fx,fy,fz,cx,cy,cz,F);
	printm4("J",J);
	printm4("F",F);
		
	clear4(Wp);

	dyn_eq(J,Wp,F,var);
	char a;
	scanf(" %c",&a);
	skew4(Wp,J,F2);
	printm4("Wp",Wp);
	printm4("F",F2);
	printm4("F2",F2);

	scanf(" %c",&a); printf("-------\n");

	jxx=1.01; jyy=1.02; jzz=1.03;
	jxy=0.1;  jxz=0.2;  jyz=0.3;

	xg=0.11; yg=0.22; zg=0.33;

	m=10;

	jtoJ(m,jxx,jyy,jzz,jxy,jyz,jxz,xg,yg,zg,J);

	fx=1; fy=2; fz=3;
	cx=1.1; cy=1.2; cz=1.3;
	actom(fx,fy,fz,cx,cy,cz,F);

	clear4(Wp);
	dyn_eq(J,Wp,F,var);
	skew4(Wp,J,F2);
	printm4("Wp",Wp);
	printm4("F",F);
	printm4("F2",F2);

	scanf(" %c",&a); printf("0-------\n");
	clear4(F);
	dyn_eq(J,Wp,F,var2);
	skew4(Wp,J,F2);
	printm4("Wp",Wp);
	printm4("F",F);
	printm4("F2",F2);

	scanf(" %c",&a); printf("1-------\n");
	Wp[Z][X]=999.; Wp[X][Z]=999.;
	Wp[X][U]=999.;
	Wp[Z][U]=999.;

	F[Y][Z]=999.; F[Z][Y]=999.;
	F[Y][X]=999.; F[X][Y]=999.;
	F[Y][U]=999.;


	dyn_eq(J,Wp,F,var3);
	skew4(Wp,J,F2);
	printm4("Wp",Wp);
	printm4("F",F);
	printm4("F2",F2);

	scanf(" %c",&a); printf("2-------\n");
	Wp[Y][Z]=999.; Wp[Z][Y]=999.;
	Wp[Y][X]=999.; Wp[X][Y]=999.;
	Wp[X][U]=999.;
	Wp[Z][U]=999.;

	F[X][Z]=999.; F[Z][X]=999.;
	F[Y][U]=999.;


	dyn_eq(J,Wp,F,var4);
	skew4(Wp,J,F2);
	printm4("Wp",Wp);
	printm4("F",F);
	printm4("F2",F2);

	scanf(" %c",&a); printf("3-------\n");
	Wp[X][Z]=999.; Wp[Z][X]=999.;
	Wp[Y][U]=999.;

	F[Y][Z]=999.; F[Z][Y]=999.;
	F[X][Y]=999.; F[Y][X]=999.;
	F[X][U]=999.;
	F[Z][U]=999.;


	dyn_eq(J,Wp,F,var5);
	skew4(Wp,J,F2);
	printm4("Wp",Wp);
	printm4("F",F);
	printm4("F2",F2);

}
