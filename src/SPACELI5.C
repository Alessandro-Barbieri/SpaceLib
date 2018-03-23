/*
------------------------------------------------------------------------------
			SPACELI5.C
			(C) copyright by

			G. Legnani, D. Manara

			version v2.2 - November 2005

			Developed for :

			Prof. G. LEGNANI
			University of Brescia
			Mechanical Eng. Department
			Via Branze 38
			25123 BRESCIA - ITALY

			giovanni.legnani @ bsing.unibs.it

------------------------------------------------------------------------------
			Don't remove this copyright notice
------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPACELIB.H"


/* == crossVtoM =========================================================== */
/* function that evaluates the cross product between vector a and vector b and
   put the result in the 3*3 upper-left part of matrix C*/

void crossVtoM(real *a, real *b, MAT C, int dim)

#define c(l,m)    (*((C)+(dim)*(l)+(m)))
{
      c(X,X) = 0.;
	  c(Y,Y) = 0.;
	  c(Z,Z) = 0.;
      c(Y,X) = a[X]*b[Y]-a[Y]*b[X]; /* cz */
	  c(X,Y) = -c(Y,X);
	  c(Z,X) = a[X]*b[Z]-a[Z]*b[X]; /* -cy */
	  c(X,Z) = -c(Z,X);
	  c(Z,Y) = a[Y]*b[Z]-a[Z]*b[Y]; /* cx */
	  c(Y,Z) = -c(Z,Y);
      #undef c
}

/* == crossMtoM =========================================================== */
/* function that evaluates the cross product between vector a and vector b stored in the
   matrices A and B and put the result in the 3*3 upper-left part of matrix C*/

void crossMtoM(MAT A, MAT B, MAT C, int dima, int dimb, int dimc)

#define a(l,m)    (*((A)+(dima)*(l)+(m)))
#define b(l,m)    (*((B)+(dimb)*(l)+(m)))
#define c(l,m)    (*((C)+(dimc)*(l)+(m)))
{
      c(X,X) = 0.;
	  c(Y,Y) = 0.;
	  c(Z,Z) = 0.;
      c(Y,X) = a(Z,Y)*b(X,Z)-a(X,Z)*b(Z,Y); /* cz */
	  c(X,Y) = -c(Y,X);
	  c(Z,X) = a(Z,Y)*b(Y,X)-a(Y,X)*b(Z,Y); /* -cy */
	  c(X,Z) = -c(Z,X);
	  c(Z,Y) = a(X,Z)*b(Y,X)-a(Y,X)*b(X,Z); /* cx */
	  c(Y,Z) = -c(Z,Y);

#undef a
#undef b
#undef c
}

/* == axis =========================================================== */
/* function that store the value of axis a into the vector A*/

void axis(int a, AXIS A)

{
	switch(a)
	  {
      case X:
		   A[X]=1.;
		   A[Y]=0.;
		   A[Z]=0.;
		  break;
      case Y:
		   A[X]=0.;
		   A[Y]=1.;
		   A[Z]=0.;
		  break;
	  case Z:
		   A[X]=0.;
		   A[Y]=0.;
		   A[Z]=1.;
		  break;
	  default:
	       A[X]=0.;
		   A[Y]=0.;
		   A[Z]=0.;
		   fprintf(stderr,"Error (axis): (%d) illegal axis number.\n",a);
		  break;

	  }

}

/* == traslat =========================================================== */
/* function that Builds the matrix m of a traslation h along a unit vector u.*/

void traslat(AXIS u, real h, MAT4 m)

{
	  unitv(u,u);

      m[X][X]=1.;m[Y][Y]=1.;m[Z][Z]=1.;m[U][U]=1.;
	  m[X][Y]=0.;m[X][Z]=0.;m[Y][X]=0.;m[Y][Z]=0.;
	  m[Z][X]=0.;m[Z][Y]=0.;m[U][X]=0.;m[U][Y]=0.;
	  m[U][Z]=0.;
      m[X][U]=u[X]*h;
	  m[Y][U]=u[Y]*h;
	  m[Z][U]=u[Z]*h;

}

/* == traslat2 =========================================================== */
/* function that Builds the matrix m of a traslation h along a frame's axes*/

void traslat2(int a, real h, MAT4 m)

{
	  AXIS u;
	  axis(a,u);
	  traslat(u,h,m);

}

/* == traslat24 =========================================================== */
/* function that Builds the matrix m of a traslation h along a frame's axes*/

void traslat24(int a, real h, POINT p, MAT4 m)

{
	  int i;
	  traslat2(a,h,m);
      for(i=X;i<=Z;i++)
		  {
			  m[i][U]=m[i][U]+p[i];
		  }
}


/* == psedot =========================================================== */
/* function that performs the Pseudo scalar product between two 4*4 matrices.*/

real psedot(MAT4 A, MAT4 B)

{
	  real ps;
	  ps=A[X][U]*B[X][U]+A[Y][U]*B[Y][U]+A[Z][U]*B[Z][U]+
	     A[Z][Y]*B[Z][Y]+A[X][Z]*B[X][Z]+A[Y][X]*B[Y][X];
	  return(ps);
}


