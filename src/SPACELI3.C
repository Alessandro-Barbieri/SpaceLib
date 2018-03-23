/*
------------------------------------------------------------------------------
		SPACELI3.C addendum to SPACELIB.C (c)

		       copyright by

		G. Legnani and R. Adamini

		last revision May 1997

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

		giovanni.legnani @ bsing.unibs.it

------------------------------------------------------------------------------
		Do not remove this copyright notice.
------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "SPACELIB.H"


/* --- Rotation and Position Matrices --- */

/* == rotat2 ============================================================== */
/* this function builds a rotation matrix describing a rotation of angle q
   around axis a. It is
      a : rotation axis (it must be either the constant X=0,Y=1, Z=2, U=3)
   The rotation matrix is stored in the 3*3 upper-left part of a dim*dim
   matrix. If a=U the rotation is assumed null */

int rotat2(int a, real q, MAT R, int dim)

#define r(i,j) (*((R)+(i)*dim+(j)))
{
	int x,y,z;
	real s,c;

	if (a<X || a>U)
	{
		fprintf(stderr,"Error (rotat2) : illegal axis\n");
		return(NOTOK);
	}
	if (a == U) 		  /* a==U traslation only */
	{
		#ifdef _BORLAND_
			*R      =  1; *(R+1)  =  0; *(R+2)  = 0;
			r(Y, X) =  0; r(Y, Y) =  1; r(Y, Z) = 0;
			r(Z, X) =  0; r(Z, Y) =  0; r(Z, Z) = 1;

		#else
			r(X, X) =  1; r(X, Y) =  0; r(X, Z) = 0;
			r(Y, X) =  0; r(Y, Y) =  1; r(Y, Z) = 0;
			r(Z, X) =  0; r(Z, Y) =  0; r(Z, Z) = 1;
		#endif
	}
	else			  /* a==X,Y or Z roto-traslation */
	{
		x=(a+1) % 3; y=(a+2) % 3; z=a;
		c = (real) cos(q);  s = (real) sin(q);

		r(x, x) =  c; r(x, y) = -s; r(x, z) = 0;
		r(y, x) =  s; r(y, y) =  c; r(y, z) = 0;
		r(z, x) =  0; r(z, y) =  0; r(z, z) = 1;
	}

#undef r /* Bug fixed 03.04.97 AmGhiPu */

return (OK);
}

/* == rotat24 ============================================================= */
/* this function builds a position matrix of a frame whose origin is stored
   in point O and rotated of angle q around axis a. It is
      O : homogeneous coordinates
      a : rotation axis (it must be either the constant X=0,Y=1, Z=2, U=3)
   If a=U the rotation is assumed null */

int rotat24(int a, real q, POINT O, MAT4 m)
{
	if (a<X || a>U)
	{
		fprintf(stderr,"Error (rotat24) : illegal axis\n");
		return(NOTOK);
	}
	if (O[U] != 1)
	{
		fprintf(stderr,"Error (rotat24) : illegal point");
		return(NOTOK);
	}
	rotat2(a, q, M m, 4);

	/* traslation and constant parts */
	m[X][U] = O[X]; m[Y][U] = O[Y]; m[Z][U] = O[Z];
	m[U][X] = m[U][Y] = m[U][Z] = 0 ; m[U][U] =  1 ;

	return (OK);
}


/* --- Speed and Acceleration Matrices --- */

/* == makeL =============================================================== */
/* this function builds a L matrix describing a rotation or a traslation
   about axis whose unit vector is u and which passes through point P.
   pitch is the pitch of the axis. It is jtype==Rev (rotation) or jtype==Pri
   (traslation). If jtype==Pri, pitch is ignored */

int makeL(int jtype, AXIS u, real pitch, POINT P, MAT4 L)
{
	VECTOR v;

	clear4(L);

	if (jtype==Pri)      /* traslation only */
	{
		L[X][U]=u[X]; L[Y][U]=u[Y]; /*L[Y][U]=u[X]; bug fixed */
		L[Z][U]=u[Z];
	}
	else if (jtype==Rev) /* revolute or screw */
	{
		L[Y][Z]=-(L[Z][Y]=u[X]);
		L[Z][X]=-(L[X][Z]=u[Y]);
		L[X][Y]=-(L[Y][X]=u[Z]);
		cross(u,P,v);
		L[X][U]=-v[X]+pitch*u[X];
		L[Y][U]=-v[Y]+pitch*u[Y];;
		L[Z][U]=-v[Z]+pitch*u[Z];;
	}
	else
	{
		fprintf(stderr,"Error (makeL) : illegal joint type\n");
		return(NOTOK);
	}

	return(OK);
}

/* == makeL2 ============================================================== */
/* this function builds a L matrix describing a rotation or a traslation
   about axis a with pitch pitch. It is jtype==Rev (rotation) or jtype==Pri
   (traslation). If jtype==Pri, pitch is ignored */

int makeL2(int jtype, int a, real pitch, POINT P, MAT4 L)
{
	int x,y,z;

	if (a<X || a>Z)
	{
		fprintf(stderr,"Error (makeL2) : illegal axis\n");
		return(NOTOK);
	}

	clear4(L);

	if (jtype==Pri)      /* traslation only */
	{
		L[a][U]=1;
	}
	else if (jtype==Rev) /* revolute or screw */
	{
		x=(a+1) % 3; y=(a+2) % 3; z=a;

		L[y][x]=1; L[x][y]=-1;
		L[x][U] =  P[y];
		L[y][U] = -P[x];
		L[z][U] = pitch;
	}
	else
	{
		fprintf(stderr,"Error (makeL2) : illegal joint type\n");
		return(NOTOK);
	}

	return(OK);
}

/* == Wtovel ============================================================== */
/* extracts from a speed matrix W the screw parameters:
      u     : axis of rotation
      omega : angular velocity
      vel   : linear velocity along u
      P     : a point of the axis
*/

void Wtovel(MAT4 W, AXIS u, real *omega, real *vel, POINT P)
{
	AXIS v;

	P[U]=1;

	u[X]=W[Z][Y];
	u[Y]=W[X][Z];
	u[Z]=W[Y][X];
	*omega=unitv(u,u);
	if(*omega!=0)
	{
		v[X]=W[X][U]; v[Y]=W[Y][U]; v[Z]=W[Z][U];
		cross(u,v,P);
		P[X] /= *omega;
		P[Y] /= *omega;
		P[Z] /= *omega;
		*vel=dot(u,v);
	}
	else
	{
		u[X]=W[X][U]; u[Y]=W[Y][U]; u[Z]=W[Z][U];
		P[X]=P[Y]=P[Z]=0;
		*vel=unitv(u,u);
	}
}


/* --- Conversion from Cardan Angles to Matrices, Position --- */

/* == cardantor =========================================================== */
/* Cardan or Euler angles to rotation matrix. Builds a rotation matrix
   starting from the Cardan or Euler angles and stores it in the 3*3 upper-
   left sub-matrix of a dim*dim matrix A. It is
      i,j,k : rotation axes (their value must be the constant X=0, Y=1 or Z=2)
      q : vector containing 1st, 2nd and 3rd angle
*/

int cardantor(real *q, int i, int j, int k, MAT A, int dim)

#define a(ii,jj) *(A+(ii*dim)+jj)
#define alpha  *(q)
#define beta  *(q+1)
#define gamma *(q+2)
{
	int sig,l;
	real sa,sb,sc, ca,cb,cc;

	sa=(real) sin(alpha); sb=(real) sin(beta); sc=(real) sin(gamma);
	ca=(real) cos(alpha); cb=(real) cos(beta); cc=(real) cos(gamma);

	if(i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (cardantor) : illegal axis\n");
		return(NOTOK);
	}

	if((j-i+3)%3==1)  /* cyclic */
		sig=1;
	else              /* anti-cyclic */
		sig=-1;

	if(i!=k)          /* Cardan convention */
	{
		a(i,i)=cb*cc;
		a(i,j)=-sig*cb*sc;
		a(i,k)=sig*sb;

		a(j,i)=sa*sb*cc+sig*ca*sc;
		a(j,j)=-sig*sa*sb*sc+cc*ca;
		a(j,k)=-sig*sa*cb;

		a(k,i)=-sig*ca*sb*cc+sa*sc;
		a(k,j)=ca*sb*sc+sig*sa*cc;
		a(k,k)=ca*cb;
	}
	else              /* Euler convention */
	{
		l=3-i-j;
		a(i,i)=cb ;
		a(i,j)=sb*sc;
		a(i,l)=sig*sb*cc;

		a(j,i)=sa*sb;
		a(j,j)=-sa*cb*sc+ca*cc;
		a(j,l)=-sig*( ca*sc+sa*cb*cc );

		a(l,i)=-sig*ca*sb;
		a(l,j)=sig*( ca*cb*sc+sa*cc );
		a(l,l)=-sa*sc+ca*cb*cc;
	}

#undef a
#undef alpha
#undef beta
#undef gamma

	return(OK);
}

/* == rtocardan =========================================================== */
/* rotation matrix to Cardan or Euler angles. Extracts the Cardan or Euler
   angles from the 3*3 upper-left sub-matrix of a dim*dim matrix. It is
      i,j,k : rotation axes (their value must be the constant X=0, Y=1 or Z=2)
      q1[3] : first solution;
      q2[3] : second solution;
*/

int rtocardan(MAT R, int dim, int i, int j, int k, real q1[3], real q2[3])

#define r(ii,jj) (*((R)+(ii)*dim+(jj)))
{
	int sig,l;

	if(i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (rtocardan) : illegal axis\n");
		return(NOTOK);
	}
	if((j-i+3)%3==1)  /* cyclic */
		sig=1;
	else              /* anti-cyclic */
		sig=-1;
	if(i!=k)          /* i!=k Cardan convention */
	{
		q1[0]=(real) atan2(-sig*r(j,k),r(k,k));
		q1[1]=(real)  asin( sig*r(i,k));
		q1[2]=(real) atan2(-sig*r(i,j),r(i,i));

		q2[0]=(real) atan2(sig*r(j,k),-r(k,k));
		q2[1]=(real) fmod(PIG-asin(sig*r(i,k))+PIG,PIG2)-PIG;
		q2[2]=(real) atan2(sig*r(i,j),-r(i,i));

	}
	else              /* i==k Euler convention */
	{
		l=3-i-j;
		q1[0]=(real) atan2(r(j,i),-sig*r(l,i));
		q1[1]=(real) acos(r(i,i));
		q1[2]=(real) atan2(r(i,j), sig*r(i,l));

		q2[0]=(real) atan2(-r(j,i), sig*r(l,i));
		q2[1]=(real) -acos(r(i,i));
		q2[2]=(real) atan2(-r(i,j),-sig*r(i,l));
	}

#undef r

	return(OK);
}

/* == cardantoM =========================================================== */
/* build the position matrix m of a frame whose origin is O and whose
   orientation is specified by a Euler/Cardanic convention:
      * 3 rotations around axis i, j, k
      * the angles value are contained in array q
*/

void cardantoM(real *q, int i, int j, int k, POINT O, MAT4 m)
{
	int l;

	cardantor4(q,i,j,k,m);
	for(l=0;l<4;l++)
		m[l][U]=O[l];
	m[U][X]=m[U][Y]=m[U][Z]=0;
}


/* --- Conversion from Cardan Angles to Matrices, Vel. and Acc. --- */

/* == cardantoW =========================================================== */
/* builds the velocity matrix W of a frame whose origin is O and whose
   orientation is specified by a Euler/Cardanic convention:
      * 3 rotations around axis i, j, k
      * the angle values are contained in array q
      * the angular speeds are contained in array qp
*/

void cardantoW(real *q, real *qp, int i, int j, int k, POINT O, MAT4 W)
{
	int l,m;

	cardanto_omega4(q,qp,i,j,k,W);
	for(l=X;l<=Z;l++)
	{
		W[l][U]=0;
		/*W[i][U]=0; bug fixed 13/05/97 AmGhiPu */
		for(m=X;m<=Z;m++)
			W[l][U] -= W[l][m]*O[m];
	}
	W[U][X]=W[U][Y]=W[U][Z]=W[U][U]=0;
}

/* == cardanto_OMEGA ====================================================== */

void cardanto_OMEGA(real *q, real *qp, int i, int j, int k, real *omega)
{
	real mat[9];

	cardantol(q,i,j,k,mat,3);
	molt(mat,qp,omega,3,3,1);
}

/* == cardanto_omega ====================================================== */

void cardanto_omega(real *q, real *qp, int i, int j, int k, MAT A, int dim)
{
	real mat[9], omega[3];

	cardantol(q,i,j,k,mat,3);
	molt(mat,qp,omega,3,3,1);    /* evaluates angular velocity */
	vtom(omega,A,dim);           /* store a.v. in 3*3 skew simm. matrix */
}

/* == cardantoH =========================================================== */
/* builds the acceleration matrix H of a frame whose origin is O and whose
   orientation is specified by a Euler/Cardanic convention:
      * 3 rotations, each around an axis i, j, k
      * the angles value are contained in array q
      * the angular speeds are contained in array qp
      * the angular accelerations are contained in array qpp
*/

void cardantoH(real *q,real *qp,real *qpp, int i, int j, int k, POINT O,
	       MAT4 H)
{
	int l,m;

	cardanto_G4(q,qp,qpp,i,j,k,H);
	for(l=X;l<=Z;l++)
	{
		H[l][U]=0;
		/* H[i][U]=0; bug fixed 13/05/97 AmGhiPu */
		for(m=X;m<=Z;m++)
			H[l][U] -= H[l][m]*O[m];
	}
	H[U][X]=H[U][Y]=H[U][Z]=H[U][U]=0;
}

/* == cardanto_G ========================================================== */

void cardanto_G(real *q, real *qp, real *qpp, int i, int j, int k, MAT A,
		int dim)

#define a(ii,jj) (*((A)+(ii)*dim+(jj)))
{
	real mat[9],mat1[9];
	real omega[3],omegapto[3],wprod[3];
	real buffer[3],buffer1[3];

	cardantol(q,i,j,k,mat,3);
	molt(mat,qp,omega,3,3,1);    /* evaluates angular velocity */
	cardantoWPROD(q,i,j,k,mat1,3);
	wprod[X] = qp[Y]*qp[Z];
	wprod[Y] = qp[X]*qp[Z];
	wprod[Z] = qp[X]*qp[Y];
	molt(mat,qpp,buffer,3,3,1);
	molt(mat1,wprod,buffer1,3,3,1);
	sum(buffer,buffer1,omegapto,3,1);
	vtom(omegapto,A,dim);        /* ang. acc. into 3*3 skew-simm. mat. */
	#ifdef _BORLAND_
		*A += -omega[Y]*omega[Y]-omega[Z]*omega[Z];
		a(Y,Y) += -omega[X]*omega[X]-omega[Z]*omega[Z];
		a(Z,Z) += -omega[X]*omega[X]-omega[Y]*omega[Y];
	#else
		a(X,X) += -omega[Y]*omega[Y]-omega[Z]*omega[Z];
		a(Y,Y) += -omega[X]*omega[X]-omega[Z]*omega[Z];
		a(Z,Z) += -omega[X]*omega[X]-omega[Y]*omega[Y];
	#endif

	for(i=X;i<=Z;i++)
		for(j=X;j<=Z;j++)
			if(i!=j)
				a(i,j) += omega[i]*omega[j];
#undef a
}

/* == cardanto_OMEGAPTO =================================================== */

void cardanto_OMEGAPTO(real *q, real *qp, real *qpp, int i, int j, int k,
		       real *omegapto)
{
	real mat[9],mat1[9];
	real wprod[3];
	real buffer[3],buffer1[3];

	cardantol(q,i,j,k,mat,3);
	cardantoWPROD(q,i,j,k,mat1,3);
	wprod[X] = qp[Y]*qp[Z];
	wprod[Y] = qp[X]*qp[Z];
	wprod[Z] = qp[X]*qp[Y];
	molt(mat,qpp,buffer,3,3,1);
	molt(mat1,wprod,buffer1,3,3,1);
	sum(buffer,buffer1,omegapto,3,1);
}

/* == cardantol =========================================================== */

int cardantol(real *q, int i, int j, int k, MAT R, int dim)

#define r(ii,jj) *(R+(ii*dim)+jj)
#define alpha  *(q)
#define beta  *(q+1)
#define gamma *(q+2)
{
	int sig,l;
	real sa,sb, ca,cb;

	sa=(real) sin(alpha); sb=(real) sin(beta);
	ca=(real) cos(alpha); cb=(real) cos(beta);

	if(i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (cardantol) : illegal axis\n");
		return(NOTOK);
	}
	if((j-i+3)%3==1)  /* cyclic */
		sig=1;
	else              /* anti-cyclic */
		sig=-1;

	if(i!=k) 	  /* Cardan convention */
	{
		r(i,X)=1; r(i,Y)=0;      r(i,Z)=sig*sb;
		r(j,X)=0; r(j,Y)=ca;     r(j,Z)=-sig*sa*cb;
		r(k,X)=0; r(k,Y)=sig*sa; r(k,Z)=ca*cb;
	}
	else     	  /* Euler convention */
	{
		l=3-i-j;
		r(i,X)=1; r(i,Y)=0; 	 r(i,Z)=cb;
		r(j,X)=0; r(j,Y)=ca; 	 r(j,Z)=sa*sb;
		r(l,X)=0; r(l,Y)=sig*sa; r(l,Z)=-sig*ca*sb;
	}

#undef r
#undef alpha
#undef beta
#undef gamma

	return(OK);
}

/* == cardantoWPROD ======================================================= */

int cardantoWPROD(real *q, int i, int j, int k, MAT R, int dim)

#define r(ii,jj) *(R+(ii*dim)+jj)
#define alpha  *(q)
#define beta  *(q+1)
#define gamma *(q+2)
{
	int sig,l;
	real sa,sb, ca,cb;

	sa=(real) sin(alpha); sb=(real) sin(beta);
	ca=(real) cos(alpha); cb=(real) cos(beta);

	if(i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (cardantoWPROD): illegal axis\n");
		return(NOTOK);
	}

	if((j-i+3)%3==1)  /* cyclic */
		sig=1;
	else              /* anti-cyclic */
		sig=-1;

	if(i!=k) 	  /* Cardan convention */
	{
		r(i,X)=sig*cb;    r(i,Y)=0;          r(i,Z)=0;
		r(j,X)=sig*sa*sb; r(j,Y)=-sig*ca*cb; r(j,Z)=-sa;
		r(k,X)=-ca*sb;    r(k,Y)=-sa*cb;     r(k,Z)=sig*ca;
	}
	else              /* Euler convention */
	{
		l=3-i-j;
		r(i,X)=-sb;        r(i,Y)=0;         r(i,Z)=0;
		r(j,X)=sa*cb;      r(j,Y)=ca*sb;     r(j,Z)=-sa;
		r(l,X)=-sig*ca*cb; r(l,Y)=sig*sa*sb; r(l,Z)=sig*ca;
	}

#undef r
#undef alpha
#undef beta
#undef gamma

	return(OK);
}


/* --- Construction of Frames Attached to Points or Vectors --- */

/* == frameP ============================================================== */
/* builds a rotation matrix from three points. Axis a1 from P1 toward point
   P2, axis a2 from P1 toward point P3. Axis a1 has priority. a1,a2 must be
   either the constant X=0, Y=1 or Z=2 (a1!=a2). The rotation matrix is stored
   in the 3*3 upper-left part of the dim*dim matrix A */

int frameP(POINT P1, POINT P2, POINT P3, int a1, int a2, MAT A, int dim)
{
	AXIS a,b;

	if((a1==a2) || (min(a1,a2)<X) || (max(a1,a2)>Z))
	{
		fprintf(stderr,"Error (frameP) : illegal axis\n");
		return(NOTOK);
	}

	a[X]=P2[X]-P1[X];  a[Y]=P2[Y]-P1[Y]; a[Z]=P2[Z]-P1[Z];
	b[X]=P3[X]-P1[X];  b[Y]=P3[Y]-P1[Y]; b[Z]=P3[Z]-P1[Z];
	frameV(a, b, a1, a2, M A, dim);

/* #undef RR 	Bug fixed 03.04.97 AmGhiPu */

	return(OK);
}

/* == frame4P ============================================================= */
/* builds a frame from three points. Origin is in point P1. Axis a1 is toward
   point P2. Axis a2 is toward point P3. Axis a1 has priority. a1,a2 must be
   either the constant X=0, Y=1 or Z=2 (a1!=a2) */

void frame4P(POINT P1, POINT P2, POINT P3, int a1, int a2, MAT4 m)
{
	int i;

	frameP(P1,P2,P3, a1,a2, M m,4);
	for(i=X;i<=U;i++)
		m[i][U]=P1[i];
	m[U][X]=m[U][Y]=m[U][Z]=0;
}

/* == frameV ============================================================== */
/* builds a rotation matrix from two vectors. Axis a1 directed as v1, axis a2
   directed as v2. Axis a1 has priority. Third axis as v1xv2. a1,a2 must be
   either the constant X=0, Y=1 or Z=2 (a1!=a2). The rotation matrix is
   stored in the 3*3 upper-left part of the dim*dim matrix A */

int frameV(VECTOR v1, VECTOR v2, int a1, int a2, MAT A, int dim)

#define r(i,j) (*((A)+(i)*dim+(j)))
{
	int a3;
	AXIS a,b,c;
	int i;

	a3=3-a1-a2;

	if((a1==a2) || (min(a1,a2)<X) || (max(a1,a2)>Z))
	{
		fprintf(stderr,"Error (frameV) : illegal axes\n");
		return(NOTOK);
	}

	unitv(v1,a);
	if((a1==X && a2==Y) || (a1==Y && a2==Z) || (a1==Z && a2==X))
	{
		cross(a,v2,c); unitv(c,c); cross(c,a,b);
	}
	else
	{
		cross(v2,a,c); unitv(c,c); cross(a,c,b);
	}
	for(i=0;i<3;i++)
	{
		r(i,a1)=a[i]; r(i,a2)=b[i]; r(i,a3)=c[i];
	}

#undef r

	return(OK);
}

/* == frame4V ============================================================= */
/* builds a frame from two vectors and one point. Origin is in point P1.
   Axis a1 is parallel to v1. Axis a2 is directed as v2. Axis a1 has priority.
   a1,a2 must be either the constant X=0, Y=1 or Z=2 (a1!=a2) */

void frame4V(POINT P1, VECTOR v1, VECTOR v2, int a1, int a2, MAT4 m)
{
	int i;

	frameV(v1,v2, a1, a2, M m,4);
	for(i=X;i<=U;i++)
		m[i][U]=P1[i];
	m[U][X]=m[U][Y]=m[U][Z]=0;
}


/* --- Working with Points Lines and Planes, Operations on Points --- */

/* == angle =============================================================== */
/* angle between three points */

real angle(POINT P1, POINT P2, POINT P3)
{
	VECTOR a,b;

	vect(P1,P2,a);  unitv(a,a);
	vect(P3,P2,b);  unitv(b,b);
	return acos(dot(a,b));
}

/* == dist ================================================================ */
/* distance between two points */

real dist(POINT P1, POINT P2)
{
	real d;

	d= (real) sqrt((P1[X]-P2[X])*(P1[X]-P2[X])
	  +(P1[Y]-P2[Y])*(P1[Y]-P2[Y])
	  +(P1[Z]-P2[Z])*(P1[Z]-P2[Z]));

	return d;
}

/* == middle ============================================================== */
/* middle point between two points */

void middle(POINT P1, POINT P2, POINT P)
{
	P[X]=(real) 0.5 * (P1[X]+P2[X]);
	P[Y]=(real) 0.5 * (P1[Y]+P2[Y]);
	P[Z]=(real) 0.5 * (P1[Z]+P2[Z]);
	P[U]=1;
}

/* == vect ================================================================ */
/* vector between two points */

void vect(POINT P1, POINT P2, VECTOR v)
{
	v[X]=P1[X]-P2[X];
	v[Y]=P1[Y]-P2[Y];
	v[Z]=P1[Z]-P2[Z];
}


/* --- Operations on Matrices and Vectors, Gen. Operations on Matrices --- */

/* == clear =============================================================== */
/* clears a id*jd matrix */

void clear(MAT A, int id, int jd)

#define a(i,j) (*((A)+(i)*jd+(j)))
{
	int i,j;

	for(i=0;i<id;i++)
		for(j=0;j<jd;j++)
			a(i,j)=0;

#undef a
}

/* == idmat =============================================================== */
/* makes unitary a ijd*ijd matrix A */

void idmat(MAT A, int ijd)

#define a(i,j) (*((A)+(i)*ijd+(j)))
{
	int i,j;

	for(i=0;i<ijd;i++)
		for(j=0;j<ijd;j++)
			if(i!=j)
				a(i,j)=0;
			else
				a(i,j)=1;

#undef a
}


/* --- Copy Functions --- */

/* == mmcopy ============================================================== */

void mmcopy(MAT A, MAT B, int d1, int d2, int im, int jm)

#define a(i,j)  (*((A)+(d1)*(i)+(j)))      /* attention to the define */
#define b(i,j)  (*((B)+(d2)*(i)+(j)))      /* attention to the define */
{
	int i,j;

	for (i=0;i<im;i++)
		for (j=0;j<jm;j++)
			b(i,j)=a(i,j);

#undef a
#undef b
}


/* --- Print Functions --- */

/* == prmat =============================================================== */
/* writes a position matrix with the convention of GRP_MAN graphics post
   processor */

void prmat(FILE *grpout, char str[], MAT4 m)
{
	int i,j;
	real rr;

	fprintf(grpout,"\n%s\n",str);

	for(j=0;j<4;j++)
	{
		if (j!=2)
		{
			for(i=0;i<3;i++)
			{
				rr = m[i][j];
				/*if(j==3) rr = rr*100/1.8;*/ /* debug only */
				fprintf(grpout,"%8.3f ",rr);
			}
		}
		fprintf(grpout,"\n");
	}
	fprintf(grpout,"\n\n");
}