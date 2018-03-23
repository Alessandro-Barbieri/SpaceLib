/*
------------------------------------------------------------------------------
		SPACELIB.C ©
		copyright by

		G. Legnani and R. Faglia

		preliminary version v0.1 - Oct. 1990
		revision (bug fixed) 25.4.1993
		revision (bug fixed) 03.4.1997
                last revision July 1998

		Università degli Studi di Brescia
		Dipartimento di Ingegneria Meccanica ed Industriale
		Via Branze 38
		25123 BRESCIA - ITALY

		giovanni.legnani @ unibs.it

------------------------------------------------------------------------------
		Do not remove this copyright notice.
------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "SPACELIB.H"


/* --- Rotation and Position Matrices --- */

/* == dhtom =============================================================== */
/* builds D-H transformation matrix from the extended D-H parameters & the
   joint type */

int dhtom(int jtype, real theta, real d, real b, real a, real alpha, real q,
	  MAT4 m)
{
	real ca,sa,ct,st;

	if (jtype==Pri)                        /* Prismatic joint */
		d+=q;
	else                                   /* Revolute joint */
	{
		if (jtype!=Rev)
		{
			fprintf(stderr,"Error (dhtom) : illegal joint type\n");
			return (NOTOK);
		}
		theta += q;
	}

	ca=(real) cos(alpha);  sa=(real)sin(alpha);
	ct=(real) cos(theta); st=(real)sin(theta);

	m[X][X] = ct;      m[X][Y] = -st*ca;  m[X][Z] = st*sa;
	m[Y][X] = st;      m[Y][Y] = ct*ca;   m[Y][Z] = -ct*sa;
	m[Z][X] =(real)0.; m[Z][Y] = sa;      m[Z][Z] = ca;
	m[U][X] =(real)0.; m[U][Y] =(real)0.; m[U][Z] = (real)0.;

	m[X][U] = a*ct-b*st;
	m[Y][U] = a*st+b*ct;
	m[Z][U] = d;
	m[U][U] =(real)1.;

	return (OK);
}

/* == extract ============================================================= */
/* extracts the unit vector u and the rotational angle fi from the upper-left
   3*3 submatrix of a dim*dim matrix A. fi is assumed >= 0 and <= PIG */

void extract(MAT A, AXIS u, real *fi, int dim)

#define r(i,j) (*((A)+(j)+(i)*(dim)))      /* attention to the define */
{
	real a,b,c,s,co,t,v;
	int x,y,z;

	a = (real) 0.5*(r(Z,Y)-r(Y,Z));
	b = (real) 0.5*(r(X,Z)-r(Z,X));
	c = (real) 0.5*(r(Y,X)-r(X,Y));
	s = (real) sqrt(a*a+b*b+c*c);

	#ifdef _BORLAND_
		co = (real) (*A+r(Y,Y)+r(Z,Z)-1.)/2;
		/* r(X,X) substituted with *R ric 4.11.93*/
	#else
		co = (real) (r(X,X)+r(Y,Y)+r(Z,Z)-1.)/2;
	#endif

	*fi = (real) atan2(s,co);

	co= (real) min(1.,max(-1.,co));
	v= (real) 1.-co;

	if (aabs(s)>0.1)                            /* normal case */
	{
		u[X] = a/s; u[Y] = b/s; u[Z] = c/s;
	}
	else if ((*fi!=0) && (co>0))                /* fi nearly zero */
	{
		t=(1./v);

		#ifdef _BORLAND_
			u[X] = (real) sign(r(Z,Y)-r(Y,Z))*
			       (real) sqrt(aabs((*A-co)*t));
			/* r(X,X) substituted ric 4.11.93 */
			u[Y] = (real) sign( *(A+Z) - *(A+dim*Z) )*
			       (real) sqrt(aabs((r(Y,Y)-co)*t));
			u[Z] = (real) sign( *(A+Y) - *(A+dim*Y))*
			       (real) sqrt(aabs((r(Z,Z)-co)*t));
		#else
			u[X] = (real) sign(r(Z,Y)-r(Y,Z))*
			       (real) sqrt(aabs((r(X,X)-co)*t));
			u[Y] = (real) sign(r(X,Z)-r(Z,X))*
			       (real) sqrt(aabs((r(Y,Y)-co)*t));
			u[Z] = (real) sign(r(Y,X)-r(X,Y))*
			       (real) sqrt(aabs((r(Z,Z)-co)*t));
		#endif

	}
	else if ((*fi!=0) && (co<0))                /* fi nearly PIG (3.14) */
	{
		t= ((real)1./v);

		#ifdef _BORLAND_
			u[X]=(real) sqrt(aabs((*A-co)*t));
			u[Y]=(real) sqrt(aabs((r(Y,Y)-co)*t));
			u[Z]=(real) sqrt(aabs((r(Z,Z)-co)*t));
		#else
			u[X]=(real) sqrt(aabs((r(X,X)-co)*t));
			u[Y]=(real) sqrt(aabs((r(Y,Y)-co)*t));
			u[Z]=(real) sqrt(aabs((r(Z,Z)-co)*t));
		#endif

		if ((u[X]>=u[Y])&&(u[X]>=u[Z]))
			x=X;
		else if ((u[Y]>=u[X])&&(u[Y]>=u[Z]))
			x=Y;
		else
			x=Z;
		y=(x+1)%3; z=(x+2)%3;
		s=((r(z,y)-r(y,z))>=0 ? 1 : -1);
		u[x] *= s;
		u[y] *= (real) sign(r(y,x)+r(x,y))*s;
		u[z] *= (real) sign(r(z,x)+r(x,z))*s;
	}
	else                                        /* fi == 0 */
	{
		u[X] = u[Y] = u[Z] = (real) 0.;
	}

#undef r
}

/* == mtoscrew ============================================================ */

int mtoscrew (MAT4 Q, AXIS u, real *fi, real *h, POINT P)
{
	real a,b,c,s,co,d,v,t,t1;
	int i,j;
	int x,y,z;

	a = (real) 0.5*(Q[Z][Y]-Q[Y][Z]);
	b = (real) 0.5*(Q[X][Z]-Q[Z][X]);
	c = (real) 0.5*(Q[Y][X]-Q[X][Y]);
	s = (real) sqrt((double) (a*a+b*b+c*c));
	co = (Q[X][X]+Q[Y][Y]+Q[Z][Z]-1.)/2.;
	*fi = (real) atan2(s,co);

	co = (real) min(1.,max(-1.,co));
	v = (real) 1.-co;

	if (aabs(s)>0.1)                            /* normal case */
	{
		u[X] = a/s; u[Y] = b/s; u[Z] = c/s;
	}
	else if ((*fi!=0) && (co>0))                /* fi small */
	{
		t= ((real)1./v);
		u[X] = (real)sign(Q[Z][Y]-Q[Y][Z])* (real) sqrt(aabs((Q[X][X]-co)*t));
		u[Y] = (real)sign(Q[X][Z]-Q[Z][X])* (real) sqrt(aabs((Q[Y][Y]-co)*t));
		u[Z] = (real)sign(Q[Y][X]-Q[X][Y])* (real) sqrt(aabs((Q[Z][Z]-co)*t));
	}
	else if ((*fi!=0) && (co<0))                /* fi nearly PIG (3.14) */
	{
		t= ((real)1./v);
		u[X]=(real) sqrt(aabs((Q[X][X]-co)*t));
		u[Y]=(real) sqrt(aabs((Q[Y][Y]-co)*t));
		u[Z]=(real) sqrt(aabs((Q[Z][Z]-co)*t));
		if ((u[X]>=u[Y])&&(u[X]>=u[Z]))
			x=X;
		else if ((u[Y]>=u[X])&&(u[Y]>=u[Z]))
			x=Y;
		else
			x=Z;
		y=(x+1)%3; z=(x+2)%3;
		s=((Q[z][y]-Q[y][z])>=0 ? 1 : -1);
		u[x] *= s;
		u[y] *= (real) sign(Q[y][x]+Q[x][y])*s;
		u[z] *= (real) sign(Q[z][x]+Q[x][z])*s;
	}
	else                                        /* fi==0 */
	{
		d=(real)sqrt(Q[X][U]*Q[X][U]+Q[Y][U]*Q[Y][U]+Q[Z][U]*Q[Z][U]);
		if (d==0)
		{
			u[X] = u[Y] = u[Z]= (real) 0.;
		}
		else
		{
			u[X]=Q[X][U]/d; u[Y]=Q[Y][U]/d; u[Z]=Q[Z][U]/d;
		}
	}

	*h = u[X]*Q[X][U] + u[Y]*Q[Y][U] + u[Z]*Q[Z][U];

	P[U]= (real) 1.;
	if (v== (real) 0.)
	{
		P[X]=P[Y]=P[Z]=(real) 0.;
		if (d!=(real)0.)
			return(OK);
		else
			return(NOTOK);
	};
	t1 = (real).5/v;
	for(i=X;i<U;i++)
	{
		for(j=X,t=0;j<U;j++)
			t += -Q[j][i] * Q[j][U];
		P[i] = (t + Q[i][U])*t1;
	}
	return(OK);
}

/* == screwtom ============================================================ */
/* builds rototraslation matrix Q from: unit vector u, rotation angle fi,
   traslation h and a point P of the axis. fi is assumed >= 0 and <= PIG */

void screwtom (AXIS u, real fi, real h, POINT P, MAT4 Q)
{
	int i,j;
	double t;

	rotat(u,fi,M Q,4);  /* builds rotational sub-matrix */

	for(i=X;i<U;i++)
	{
		for(j=X,t=0;j<U;j++)
			t += -Q[i][j] * P[j];
		Q[i][U] = (real) t + P[i] + u[i] * h;
	}
	Q[U][X]=Q[U][Y]=Q[U][Z]= 0; Q[U][U]= 1;
}

/* == rotat =============================================================== */
/* builds rotational matrix from unit vector u and rotation angle fi.
   Stores the matrix in the 3*3 upper-left part of a dim*dim matrix A */

void rotat (AXIS u, real fi, MAT A, int dim)

#define a(i,j) (*((A)+(i)*(dim)+(j)))
{
	real s,v;

	s = (real) sin(fi);
	v =  1 - (real) cos(fi);


	#ifdef _BORLAND_
		*A=1+(u[X]*u[X]-1)*v;       *(A+1)=-u[Z]*s+u[X]*u[Y]*v;
		a(Y,X)=u[Z]*s+u[X]*u[Y]*v;  a(Y,Y)=1 +(u[Y]*u[Y]-1)*v;
		a(Z,X)=-u[Y]*s+u[X]*u[Z]*v; a(Z,Y)=u[X]*s+u[Y]*u[Z]*v;

		*(A+2)=u[Y]*s+u[X]*u[Z]*v;
		a(Y,Z)=-u[X]*s+u[Y]*u[Z]*v;
		a(Z,Z)=1+(u[Z]*u[Z]-1)*v;
	#else
		a(X,X)=1+(u[X]*u[X]-1)*v;   a(X,Y)=-u[Z]*s+u[X]*u[Y]*v;
		a(Y,X)=u[Z]*s+u[X]*u[Y]*v;  a(Y,Y)=1 +(u[Y]*u[Y]-1)*v;
		a(Z,X)=-u[Y]*s+u[X]*u[Z]*v; a(Z,Y)=u[X]*s+u[Y]*u[Z]*v;

		a(X,Z)=u[Y]*s+u[X]*u[Z]*v;
		a(Y,Z)=-u[X]*s+u[Y]*u[Z]*v;
		a(Z,Z)=1+(u[Z]*u[Z]-1)*v;
	#endif

#undef a
}


/* --- Speed and Acceleration matrix --- */

/* == gtom ================================================================ */
/* builds the gravity acceleration matrix starting from the three components
   of gravity acceleration with respect to the absolute frame */

void gtom(real gx, real gy, real gz, MAT4 Hg)
{
	clear4(Hg);
	Hg[X][U]=gx;
	Hg[Y][U]=gy;
	Hg[Z][U]=gz;
}

/* == WtoL ================================================================ */
/* extracts L matrix from the correspondent W matrix */

void WtoL(MAT4 W, MAT4 L)
{
	VECTOR u;
	real m,tol;

	tol=zerom();
	mtov(M W, 4, u);
	m=mod(u);
	if (m<tol)
		m = (real) sqrt(W[X][U]*W[X][U] +W[Y][U]*W[Y][U] +
				W[Z][U]*W[Z][U]);
	if (m<tol)
		clear4(L);
	else
		rmolt4(W,(real)1./m,L);
}

/* == velacctoWH ========================================================== */
/* builds both speed and acceleration matrices in local frame starting from
   the speed and the acceleration of the joint */

int velacctoWH(int jtype, real qp, real qpp, MAT4 W, MAT4 H)
{
	clear4(W); clear4(H);
	if (jtype==Pri)                        /* Prismatic joint */
	{
		W[Z][U]=qp;
		H[Z][U]=qpp;
	}
	else                                   /* Revolute joint */
	{
		if (jtype!=Rev)
		{
			fprintf(stderr,"Error (velacctoWH) : illegal joint type\n");
			return (NOTOK);
		}
		W[X][Y]= -(W[Y][X]=qp);
		H[X][X]=H[Y][Y]= -qp*qp;
		H[X][Y]= -(H[Y][X]=qpp);
	}
	return(OK);
}

/* == velacctoWH2 ========================================================= */
/* builds both speed and acceleration matrices in local frame starting from
   the speed and the acceleration of the joint and the axis around which the
   rotation is performed */

int velacctoWH2(int jtype, int a, real qp, real qpp, MAT4 W, MAT4 H)
{
	int x, y, z;

	if (a<X || a>Z)
	{
		fprintf(stderr,"Error (velacctoWH2) : illegal axis");
		return (NOTOK);
	}
	if (jtype != Rev && jtype != Pri)
	{
		fprintf(stderr,"Error (velacctoWH2) : illegal joint type");
		return (NOTOK);
	}

	x=(a+1) % 3; y=(a+2) % 3; z=a;

	clear4(W);  clear4(H);

	if (jtype == Pri)                      /* Prismatic joint */
	{
		W[z][U] = qp;                  /* W[a][U] = qp;  ric 4.11.93*/
		H[z][U] = qpp;                 /* H[a][U] = qpp; ric 4.11.93*/
	}
	else                                   /* Revolute joint */
	{
		W[x][y] = -(W[y][x] = qp);
		H[x][y] = -(H[y][x] = qpp); H[x][x] = H[y][y] = -qp*qp;
	}
	return(OK);
}


/* --- Inertial and Actions Matrices --- */

/* == actom =============================================================== */
/* builds action matrix FI starting from forces fx, fy, fz and couples cx, cy
   and cz */

void actom(real fx, real fy, real fz, real cx, real cy, real cz, MAT4 FI)
{
	FI[X][X]= (real) 0.; FI[X][Y]= -cz;      FI[X][Z]=  cy;
	FI[Y][X]=  cz;       FI[Y][Y]=(real) 0.; FI[Y][Z]= -cx;
	FI[Z][X]= -cy;       FI[Z][Y]=  cx;      FI[Z][Z]=(real) 0.;
	FI[U][X]= -fx;       FI[U][Y]= -fy;      FI[U][Z]= -fz;

	FI[X][U]=fx;
	FI[Y][U]=fy;
	FI[Z][U]=fz;
	FI[U][U]=(real)0.;
}

/* == jtoJ ================================================================ */
/* builds inertia matrix J of a body starting from the mass, the inertia
   moments and the position of the centre of mass */

void jtoJ(real mass, real jxx, real jyy, real jzz, real jxy, real jyz,
	  real jxz, real xg, real yg, real zg, MAT4 J)
{
	J[X][X]=(-jxx+jyy+jzz)/2+mass*xg*xg;
	J[Y][Y]=( jxx-jyy+jzz)/2+mass*yg*yg;
	J[Z][Z]=( jxx+jyy-jzz)/2+mass*zg*zg;

	J[X][Y]=J[Y][X]= -jxy+mass*xg*yg;
	J[X][Z]=J[Z][X]= -jxz+mass*xg*zg;
	J[Y][Z]=J[Z][Y]= -jyz+mass*yg*zg;

	J[U][U]=mass;
	J[U][X]=J[X][U]=mass*xg;
	J[U][Y]=J[Y][U]=mass*yg;
	J[U][Z]=J[Z][U]=mass*zg;
}


/* --- Matrix Transformations, Normalization --- */

/* == normal ============================================================== */
/* makes ortogonal a rotation matrix (if n==3) or the 3*3 upper-left
   sub-matrix of a n*n square matrix

   Note : a square matrix is ortogonal if its transposte and its invers are
	 equal
*/

int normal(MAT R, int n)

#define r(i,j)  (*((R)+(n)*(i)+(j)))      /* attention to the define */
{
	MAT3 RR;
	real maxx,det,t,t1,z;
	real kz=1.; /* update by joe July 1998 */
	real a,b;
	int i,j,nit=0;

	z=zerom(); /* tolerance */

	for(;;)
	{
		nit++;

		det =   r(0,0)*(r(1,1)*r(2,2)-r(2,1)*r(1,2))
		      - r(0,1)*(r(1,0)*r(2,2)-r(2,0)*r(1,2))
		      + r(0,2)*(r(1,0)*r(2,1)-r(2,0)*r(1,1));

		a = (real) 1./ pow(aabs(det),1./3.);
		b = (a*a);

		/* printf("det=%f\n",det);
		   printf("a=%f\n",a); printf("b=%f\n",b); */ /* debug only */

		RR[X][X] =  r(1,1)*r(2,2)-r(1,2)*r(2,1);
		RR[X][Y] =  r(1,2)*r(2,0)-r(1,0)*r(2,2);
		RR[X][Z] =  r(1,0)*r(2,1)-r(1,1)*r(2,0);

		RR[Y][X] = -r(0,1)*r(2,2)+r(0,2)*r(2,1);
		RR[Y][Y] = -r(0,2)*r(2,0)+r(0,0)*r(2,2);
		RR[Y][Z] = -r(0,0)*r(2,1)+r(0,1)*r(2,0);

		RR[Z][X] =  r(0,1)*r(1,2)-r(0,2)*r(1,1);
		RR[Z][Y] =  r(0,2)*r(1,0)-r(0,0)*r(1,2);
		RR[Z][Z] =  r(0,0)*r(1,1)-r(0,1)*r(1,0);

		maxx=(real)0.;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				t=(real)0.5*(r(i,j)*a+RR[i][j]*b);
				t1=aabs(t-r(i,j));
				if(t1>maxx)
					maxx=t1;
				r(i,j)=t;
			}
		}

		/* printf(" max= %f\n",maxx); */              /* debug only */
		if(nit>=8 && nit <=12) kz *= 3.;    /* if fails to converge */
											/* increase tolerence - july 1998 */
		if(maxx<=z*kz)
			return(nit);
	}

#undef r
}

/* == norm_simm_skew ====================================================== */
/* makes symmetric/antisymmetric (normalize) the n*n upper-left submatrix of
   a dim*dim square matrix. It is
		sign=1,-1  1=symmetric -1=anti-sym. */

int norm_simm_skew(MAT A, int n, int dim, int sign)

#define a(i,j)  *((A)+(dim)*(i)+(j))       /* attention to the define */
{
	real t;
	int i,j;

	if(sign!=1 && sign!=-1)   /* error in input data */
	{
		fprintf(stderr,"Error (norm_simm_skew) : illegal matrix type");
		return (NOTOK);
	}
	for(i=0;i<n;i++)          /* bug fixed on 25.4.1993 */
	{
		for(j=i+1;j<n;j++)
		{
			t=(real)0.5*(a(i,j)+sign*a(j,i));
			a(i,j)=t; a(j,i)=sign*t;
		}
	}
	if (sign==-1)
	{
		for(i=0;i<dim;i++)
			a(i,i)=(real)0.;
	}

#undef a

	return (OK);
}


/* --- Matrix transformations, Change of Reference --- */

/* == trasf_mami ========================================================== */
/* change of reference for 4*4 matrices which are contra-variant with respect
   to the row index and co-variant with respect to the column index and
   therefore transform as:
		A2 = m * A1 * m(i)       (i)=invers
   It is assumed m[3][0]==m[3][1]==m[3][2]==0 and m[3][3]==1 */

void trasf_mami(MAT4 A1, MAT4 m, MAT4 A2)
{
	int i,j,h,k;
	double t;

	for (i=0;i<3;i++)             /* 3*3 sub matrix */
		for (j=0;j<3;j++)
		{
			t=0.;
			for (h=0;h<3;h++)
				for (k=0;k<3;k++)
					t += m[i][h]*A1[h][k]*m[j][k];
			A2[i][j]=(real)t;
		}
	for (i=0;i<3;i++)             /* elements 1,2,3 of last column */
	{
		for (t=0.,j=0;j<3;j++)
			t += -A2[i][j]*m[j][U] +m[i][j]*A1[j][U];
		A2[i][U] =(real) t+ m[i][U]*A1[U][U];
	}

	A2[U][X]=A2[U][Y]=A2[U][Z]=(real)0.; A2[U][U]=A1[U][U]; /* last row */
}

/* == trasf_miam ========================================================== */
/* inverse change of reference for 4*4 transformation matrices which are
   contra-variant with respect to the row index and co-variant with respect to
   the column index. These matrices transform as:
		A2 = m(i) * A1 * m        (i)=invers
   It is assumed that  m[3][0]==m[3][1]==m[3][2]==0 and m[3][3]==1 */

void trasf_miam(MAT4 A1, MAT4 m, MAT4 A2)
{
	int i,j,h,k;
	double t;
	real tmp[3];

	for (i=0;i<3;i++)             /* 3*3 rotational sub-matrix*/
		for (j=0;j<3;j++)
		{
			t=0.;
			for (h=0;h<3;h++)
				for (k=0;k<3;k++)
					t += m[h][i]*A1[h][k]*m[k][j];
			A2[i][j]=(real)t;
		}
	for (i=0;i<3;i++)             /* temporary vector */
	{
		for (t=0.,j=0;j<3;j++)
			t += A1[i][j]*m[j][U];
		tmp[i] =(real) t;
	}
	for (i=0;i<3;i++)             /* 1,2,3 last column elements */
	{
		for (t=0.,j=0;j<3;j++)
			t += m[j][i]*(tmp[j]-A1[U][U]*m[j][U]+A1[j][U]);
		A2[i][U] = (real) t;
	}

	A2[U][X]=A2[U][Y]=A2[U][Z]=(real)0.; A2[U][U]=A1[U][U]; /* last row */
}

/* == trasf_mamt ========================================================== */
/* change of reference for (dim*dim) matrices which are contra-variant and
   therefore which transform as:
		A2 = m * A1 * m(t)  (t)=transposte
*/

void trasf_mamt(MAT A1, MAT m, MAT A2, int dim)

#define a1(i,j)  (*((A1)+(dim)*(i)+(j)))   /* attention to the define */
#define m(i,j)   (*((m )+(dim)*(i)+(j)))   /* attention to the define */
#define a2(i,j)  (*((A2)+(dim)*(i)+(j)))   /* attention to the define */
{
	int i,j,h,k;
	double t;

	for (i=0;i<dim;i++)
		for (j=0;j<dim;j++)
		{
			t=0.;
			for (h=0;h<dim;h++)
				for (k=0;k<dim;k++)
					t += m(i,h)*a1(h,k)*m(j,k);
				a2(i,j)=(real)t;
		}

#undef a1
#undef m
#undef a2
}

/* == trasf_miamit ======================================================== */
/* inverse change of reference for 4*4 contravariant matrices which transform
   as :
		A2 = m(i) * A1 * m(i)(t)   (t)=transposte, (i)=invers
   It is assumed that m[3][0]==m[3][1]==m[3][2]==0 and m[3][3]==1 */

void trasf_miamit(MAT4 A1, MAT4 m, MAT4 A2)
{
	int i,j,h,k;
	double t,t1,t2;
	double tmp[3];

	for (i=0;i<3;i++)             /* build temporary vector */
	{
		for (t=0.,j=0;j<3;j++)
			t += -m[j][i]*m[j][U];
		tmp[i] = t;
	}
	for (i=0;i<3;i++)             /* 3*3 rotational sub-matrix */
		for (j=0;j<3;j++)
		{
			t=0.; t1=0.; t2=0.;
			for (h=0;h<3;h++)
			{
				t1 += A1[U][h] * m[h][j];
				t2 += A1[h][U] * m[h][i];
				for (k=0;k<3;k++)
					t += m[h][i]*A1[h][k]*m[k][j];
			}
			A2[i][j] =(real) t +t1*tmp[i] +t2*tmp[j] +tmp[i]*tmp[j]*A1[U][U];
		}

	for (i=0;i<3;i++)             /* 3*1 & 1*3 parts */
	{
		t1=0.; t2=0;
		for (j=0;j<3;j++)
		{
			t1 += m[j][i]*A1[j][U];
			t2 += m[j][i]*A1[U][j];
		}
		A2[i][U]=(real)t1+tmp[i]*A1[U][U];
		A2[U][i]=(real)t2+tmp[i]*A1[U][U];
	}
	A2[U][U] = A1[U][U];
}


/* --- Matrix Transformations, General Operations --- */

/* == coriolis ============================================================ */
/* Coriolis' theorem:  H = H0 + H1 + 2 W0 * W1. Last line of W & H matrices is
   always null */

void coriolis(MAT4 H0, MAT4 H1, MAT4 W0, MAT4 W1, MAT4 H)
{

	double t;
	int i,j,k;

	for (i=0;i<3;i++)
	{
		for (j=0;j<4;j++)
		{
			t=0;
			for (k=0;k<3;k++)
			{
				t += W0[i][k]*W1[k][j];
			}
			H[i][j] =(real) 2.*t+H0[i][j]+H1[i][j];
		}
	}
	H[U][X]=H[U][Y]=H[U][Z]=H[U][U]=0;
}

/* == invers ============================================================== */
/* inverts transformation matrix m : mi = inv(m) */

void invers(MAT4 m, MAT4 mi)
{
	int i,j;
	real t;

	if (m==mi){
		runtime_error("invers(m,mi)","the input and the output operands 'm' and 'mi' must be different", 1);
	}

	for (i=0;i<3;i++)
	{
		for (t=(real)0.,j=0;j<3;j++)
		{
			mi[i][j]=m[j][i];
			t -= m[j][i]*m[j][U];
		}
		mi[i][U]=t;
	}
	mi[U][X] = mi[U][Y] = mi[U][Z] = (real) 0.; mi[U][U]=(real)1.;
}

/* == mtov ================================================================ */
/* 3*3 submatrix to 3 element vector. Extracts a 3 element vector from a 3*3
   skew-symmetric sub-matrix in the upper-left part of a dim*dim matrix A */

void mtov(MAT A, int dim, VECTOR v)

#define a(i,j) (*((A)+(i)*dim+(j)))
{
	v[X]=(real)0.5*(a(2,1)-a(1,2));
	v[Y]=(real)0.5*(a(0,2)-a(2,0));
	v[Z]=(real)0.5*(a(1,0)-a(0,1));

#undef a
}

/* == vtom ================================================================ */
/* 3 element vector to 3*3 submatrix. Creates a 3*3 skew-symmetric sub-matrix
   in the upper-left part of a dim*dim matrix A */

void vtom(VECTOR v, MAT A, int dim)

#define a(i,j) (*((A)+(i)*dim+(j)))
{

	#ifdef _BORLAND_
		*A=	   0; *(A+1)=-v[Z]; *(A+2)= v[Y];
		a(1,0)= v[Z]; a(1,1)=    0; a(1,2)=-v[X];
		a(2,0)=-v[Y]; a(2,1)= v[X]; a(2,2)=    0;
	#else
		a(0,0)=	   0; a(0,1)=-v[Z]; a(0,2)= v[Y];
		a(1,0)= v[Z]; a(1,1)=    0; a(1,2)=-v[X];
		a(2,0)=-v[Y]; a(2,1)= v[X]; a(2,2)=    0;
	#endif

#undef a
}

/* == skew ================================================================ */
/* evaluates C=skew(A*B) for dim*dim matrices */

void skew(MAT A, MAT B, MAT C, int dim)

#define a(i,j)   (*((A)+(dim)*(i)+(j)))    /* attention to the define */
#define b(i,j)   (*((B)+(dim)*(i)+(j)))    /* attention to the define */
#define c(i,j)   (*((C)+(dim)*(i)+(j)))    /* attention to the define */
{
	real t;
	int i,j,k;

	for(i=0;i<dim;i++)
	{
		for(j=i+1;j<dim;j++)
		{
			t=(real)0.;
			for(k=0;k<dim;k++)
				t += a(i,k)*b(k,j) -a(j,k)*b(k,i);
			c(i,j) = t; c(j,i)= -t;
		}
		c(i,i)=0;
	}

#undef a
#undef b
#undef c
}

/* == trac_ljlt4 ========================================================== */
/* evaluates the trace of L1*J*L2(t) (t) = transposte. L1,L2 are 4*4 matrices
   whose 4th row is null. J is a 4*4 square matrix */

real trac_ljlt4(MAT4 L1, MAT4 J, MAT4 L2)
{
	int i,h,k;
	double t;

	t=0.;
	for (i=0;i<3;i++) /* 4th row of L1 and 4th column of L2 assumed null */
		for (h=0;h<4;h++)
			for (k=0;k<4;k++)
				t += L1[i][h]*J[h][k]*L2[i][k];
	return ((real) t);
}


/* --- Operations on Matrices and Vectors, Matrices and Vectors Algebra --- */

/* == molt ================================================================ */
/*  matrix product : C[d1][d3]=A[d1][d2]*B[d2][d3] */

void molt(MAT A, MAT B, MAT C, int d1, int d2, int d3)

#define a(i,j)  (*((A)+(d2)*(i)+(j)))      /* attention to the define */
#define b(i,j)  (*((B)+(d3)*(i)+(j)))      /* attention to the define */
#define c(i,j)  (*((C)+(d3)*(i)+(j)))      /* attention to the define */
{
	double t;
	int i,j,k;

	if (A==C || B==C){
		runtime_error("molt(A,B,C)","the third operand C must be different from A and B", 1);
	}

	for (i=0;i<d1;i++)
	{
		for (j=0;j<d3;j++)
		{
			t=0;
			for (k=0;k<d2;k++)
			{
				t += a(i,k)*b(k,j);
			}
			c(i,j) = (real) t;
		}
	}

#undef a
#undef b
#undef c
}

/* == rmolt =============================================================== */
/* multiplies a d1*d2 matrix A by scalar r. B = r * A. It can be A == B */

void rmolt(MAT A, real r, MAT B, int d1, int d2)

#define a(i,j)  (*((A)+(d2)*(i)+(j)))      /* attention to the define */
#define b(i,j)  (*((B)+(d2)*(i)+(j)))      /* attention to the define */
{
	int i,j;

	for (i=0;i<d1;i++)
		for (j=0;j<d2;j++)
			b(i,j)=a(i,j)*r;

#undef a
#undef b
}

/* == sum ================================================================= */
/* matrix sum C[d1][d2]=A[d1][d2]+B[d1][d2] */

void sum(MAT A, MAT B, MAT C, int d1, int d2)

#define a(i,j)  (*((A)+(d2)*(i)+(j)))      /* attention to the define */
#define b(i,j)  (*((B)+(d2)*(i)+(j)))      /* attention to the define */
#define c(i,j)  (*((C)+(d2)*(i)+(j)))      /* attention to the define */
{
	int i,j;

	for (i=0;i<d1;i++)
		for (j=0;j<d2;j++)
			c(i,j)=a(i,j)+b(i,j);

#undef a
#undef b
#undef c
}

/* == sub ================================================================= */
/* matrix substraction C[d1][d2]=A[d1][d2]-B[d1][d2] */

void sub(MAT A, MAT B, MAT C, int d1, int d2)

#define a(i,j)  (*((A)+(d2)*(i)+(j)))      /* attention to the define */
#define b(i,j)  (*((B)+(d2)*(i)+(j)))      /* attention to the define */
#define c(i,j)  (*((C)+(d2)*(i)+(j)))      /* attention to the define */
{
	int i,j;

	for (i=0;i<d1;i++)
		for (j=0;j<d2;j++)
			c(i,j)=a(i,j)-b(i,j);

#undef a
#undef b
#undef c
}


/* --- Operations on Matrices and Vectors, Gen. Operations on Matrices --- */

/* == transp ============================================================== */
/* transposes a d1*d2 matrix A */

void transp(MAT A, MAT At, int d1, int d2)

#define a(i,j)  (*((A)+(d2)*(i)+(j)))      /* attention to the define */
#define at(i,j)  (*((At)+(d1)*(i)+(j)))    /* attention to the define */
{
	int i,j;

	for(i=0;i<d1;i++)
		for(j=0;j<d2;j++)
			at(j,i)=a(i,j);

#undef a
#undef at
}


/* --- Operations on Matrices and Vectors, Gen. Operations on Vectors --- */

/* == cross =============================================================== */
/* cross product c = a x b */

void cross(VECTOR a, VECTOR b, VECTOR  c)
{
	c[X] = a[Y] * b[Z] - a[Z] * b[Y];
	c[Y] = a[Z] * b[X] - a[X] * b[Z];
	c[Z] = a[X] * b[Y] - a[Y] * b[X];
}

/* == dot ================================================================= */
/* dot product */

real dot(VECTOR a, VECTOR b)
{
	return(a[X]*b[X] +a[Y]*b[Y] +a[Z]*b[Z]);
}

/* == mod ================================================================= */
/* module of a vector */

real mod(VECTOR a)
{
	return((real)sqrt(a[X]*a[X] +a[Y]*a[Y] +a[Z]*a[Z]));
}

/* == norm =============================================================== */
/* norm of a matrix (the element with max abs value) */
real norm(MAT A, int d1, int d2)

#define a(i,j)  (*((A)+(d2)*(i)+(j)))    /* attention to the define */
{
	int i,j;
	real rmax=0.;

	for (i=0;i<d1;i++)
		for (j=0;j<d2;j++)
			rmax=max(rmax,aabs(a(i,j)));
return rmax;
#undef a
}


/* == unitv =============================================================== */
/* evaluates the unit vector u of a vector v and it also returns its module */

real unitv(VECTOR v, AXIS u)
{
	double t;

	t = sqrt ((double) v[X]*v[X] +v[Y]*v[Y] +v[Z]*v[Z]);
	if(t!=0)
	{
		u[X] = v[X]/t; u[Y] = v[Y]/t; u[Z] = v[Z]/t;
	}
	else
	{
		u[X] = u[Y] = u[Z] =  0;
		return (0.); /* bug fixed 03/04/1997 AmGhiPu */
	}
	if (t>10. || t<.1)
		unitv(u,u);  /* to avoid numerical problems */
	return (real) t;
}

/* == vector ============================================================== */
/* builds a vector v whose module is m and whose unit vector is u */

void vector(AXIS u, real mdl, VECTOR v)
{
	v[X]=u[X]*mdl; v[Y]=u[Y]*mdl; v[Z]=u[Z]*mdl;
}


/* --- Copy Functions --- */

/* == mcopy =============================================================== */

void mcopy(MAT A1, MAT A2, int d1, int d2)

#define a1(i,j)  (*((A1)+(d2)*(i)+(j)))    /* attention to the define */
#define a2(i,j)  (*((A2)+(d2)*(i)+(j)))    /* attention to the define */
{
	int i,j;

	for (i=0;i<d1;i++)
		for (j=0;j<d2;j++)
			a2(i,j)=a1(i,j);

#undef a1
#undef a2
}


/* --- Print Functions --- */

/* == fprintm3 ============================================================ */

void fprintm3(FILE *out, char *s, MAT3 A)
{
	fprintf(out,"\n\n %s \n",s);
	fprintf(out,"             %8.3f  %8.3f  %8.3f\n",A[X][X],A[X][Y],A[X][Z]);
	fprintf(out,"             %8.3f  %8.3f  %8.3f\n",A[Y][X],A[Y][Y],A[Y][Z]);
	fprintf(out,"             %8.3f  %8.3f  %8.3f\n",A[Z][X],A[Z][Y],A[Z][Z]);
	fprintf(out,"\n");
}

/* == fprintm4 ============================================================ */

void fprintm4(FILE *out, char *s, MAT4 A)
{
	fprintf(out,"\n\n %s \n",s);
	fprintf(out,"            %8.3f  %8.3f  %8.3f  %8.3f\n"
					    ,A[X][X],A[X][Y],A[X][Z],A[X][U]);
	fprintf(out,"            %8.3f  %8.3f  %8.3f  %8.3f\n"
					    ,A[Y][X],A[Y][Y],A[Y][Z],A[Y][U]);
	fprintf(out,"            %8.3f  %8.3f  %8.3f  %8.3f\n"
					    ,A[Z][X],A[Z][Y],A[Z][Z],A[Z][U]);
	fprintf(out,"            %8.3f  %8.3f  %8.3f  %8.3f\n"
					    ,A[U][X],A[U][Y],A[U][Z],A[U][U]);
	fprintf(out,"\n");
}

/* == fprintv ============================================================= */

void fprintv(FILE *out, char *s, real *v, int n)
{
	int i;

	fprintf(out,"\n %s \n   ",s);
	for (i=0;i<n;i++)
		fprintf(out,"%8.3f",v[i]);
	fprintf(out,"\n");
}


/* --- Machine Precision Functions --- */

/* == zerom =============================================================== */

real zerom(void)
{
	real  a=(real)1., b=(real)-1;
	static real dzero= (real) -1.;

	if(dzero>(real)0.)
		return((real)dzero);
	loop:
		if (1.==b)
		{
			dzero=2*a;
			return(dzero);
		}

	a /=2.; b=a+1.;

	goto loop;
}

/* == fzerom ============================================================== */

float fzerom(void)
{
	float a=(float)1.,b=(float)-1.;
	static float dzero= (float)-1.;

	if(dzero>(float)0)
		return(dzero);

	loop:
		if (1.==b)
		{
			dzero=2*a;
			return(dzero);
		}
		a /=2; b=a+1;

	goto loop;
}

/* == dzerom ============================================================== */

double dzerom(void)
{
	double a=1.,b=-1.;
	static double dzero= -1.;

	if(dzero>0.)
		return(dzero);

	loop:
		if (1.==b)
		{
			dzero=2*a;
			return(dzero);
		}
		a /=2.; b=a+1.;

	goto loop;
}

/* == runtime_error ======================================================== october 2005.*/
/*
     prints an error messages and exit if requested (ex>0)

     name: name of the routine in which the error is produced
	 msg:  description of the error
	 ex:   exit? (exit if ex>0)
*/
void runtime_error(char *name, char *msg, int ex)
{
	fprintf(stderr,"\n\n*** runtime error spacelib/%s\n",name);
	fprintf(stderr,"%s\n\n",msg);
	if (ex){
   		fprintf(stderr,"hit any key to exit\a\a"); getch(); 
		exit(ex);
	}
	printf("hit any key to continue"); getch(); 	
}