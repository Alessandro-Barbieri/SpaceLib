/*
------------------------------------------------------------------------------
			SPACELI4.C
			(C) copyright by

			D.Amadori, P.Ghislotti, G.Pugliese

			version v1.0 - April 1997
                        revision July 1998
			patch (function pseudo_inv) by g. legnani - November 2001
			patch (function interslpl) by g. legnani - January 2003
			patch (function intermediate) by d. manara - January 2005

			Developed for :

			Prof. G. LEGNANI
			University of Brescia
			Mechanical Eng. Department
			Via Branze 38
			25123 BRESCIA - ITALY

			giovanni.legnani @ unibs.it

------------------------------------------------------------------------------
			Don't remove this copyright notice
------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPACELIB.H"


/* --- Rotation and Position Matrices --- */

/* == rotat34 ============================================================= */
/*  this function builds a position matrix of a frame whose origin is
    initially stored in point P and rotated of angle q about axis a. It is
      P : starting homogeneous coordinates of the rotating frame origin
      a : rotation axis (it must be either the constant X=0,Y=1, Z=2, U=3)
    If a=U the rotation is assumed null */

void rotat34(int a, real q, POINT P, MAT4 m)
{
	MAT3 R;
	POINT O;

	rotat24(a,q,P,m);
	mmcopy(M m,M R,4,3,3,3);
	moltmv3(M R,P,O);
	vmcopy(M O,3,4,Col,M m,4,4);
}


/* --- Speed and Acceleration Matrices --- */

/* == Gtomegapto ========================================================== */
/* function that extracts the angular velocity vector from the 3x3 upper-left
   submatrix G of the acceleration matrix H
    Input parameters
      G        : 3x3 upper-left submatrix of the acceleration matrix H
    Output parameters
      omegapto : angular acceleration vector
*/

void Gtomegapto(MAT3 G, VECTOR omegapto)
{
	MAT3 Gt;              /* transpose of matrix G */
	MAT3 OMEGAPTO;        /* angular velocity matrix */

	transp3(G,Gt);
	sub3(G,Gt,OMEGAPTO);
	rmolt3(OMEGAPTO,0.5,OMEGAPTO); /* OMEGAPTO=1/2(G-G(t)) (t)= transp. */
	mtov3(omegapto,OMEGAPTO);      /* a. v. matrix to vector */
}

/* == velacctoWH3 ========================================================= */
/* function that builds the velocity and position matricx of a body moving
   around an axis parallel to an axis of the reference frame
    Input parameters
      jtype = defines the joint type : prismatic or revolute (Rev or Pri)
      a     = defines the axis of reference frame parallel to the moving axis
	      (X,Y o Z)
      qp    = module of the angular velocity vector
      qpp   = module of the angular acceleration vector
      O     = point of the motion axis considered as origin of the frame put
	      on the moving body
    Output parameters
      W     = velocity matrix
      H     = acceleration matrix
*/

void velacctoWH3(int jtype, int a, real qp, real qpp, POINT O, MAT4 W,
		 MAT4 H)
{
	VECTOR vel,acc,omega,omegapto,a_norm;
	int x,y,z;

	clear4(W); clear4(H);
	velacctoWH2(jtype,a,qp,qpp,W,H);
	if (jtype == Rev)
	{
		x=(a+1)%3; y=(a+2)%3; z=a;
		clearv(omega); clearv(omegapto);
		omega[z]=qp; omegapto[z]=qpp;
		cross(omega,O,vel);
		cross(omega,vel,a_norm); /* norm. comp. of the acc. vector */
		cross(omegapto,O,acc); /* tang. comp. of the acc. vector */
		sumv(acc,a_norm,acc);
		W[x][U]=-vel[x]; W[y][U]=-vel[y]; W[z][U]=-vel[z];
		H[x][U]=-acc[x]; H[y][U]=-acc[y]; H[z][U]=-acc[z];
	}
}


/* --- Conversion form Cardan Angles to Matrices, Position --- */

/* == Mtocardan =========================================================== */
/* function that extracts the Euler/Cardan angles from a position matrix m
    Input parameters
      m     : position matrix
      i,j,k : constants that define the rotation sequence
	      (they must be X,Y or Z)
    Output parameters
      q1,q2 : arrays of 3 real that contain the two Euler/Cardan angles
	      solution
*/

int Mtocardan(MAT4 m, int i, int j, int k, real q1[3], real q2[3])
{
	MAT3 R;

	if (i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (Mtocardan) : illegal axis\n");
		return(NOTOK);
	}
	mcopy43(m,R);
	rtocardan3(R,i,j,k,q1,q2);
					/* sing. config. */
	if ( ((i!=k) && ((cos(q1[1])==0) || (cos(q2[1])==0))) ||
	     ((i==k) && ((sin(q1[1])==0) || (sin(q2[1])==0))) )
	{
		fprintf(stderr,"Error (Mtocardan) : there is a singularity\n");
		return(NOTOK);
	}
	else
		return(OK);
}


/* --- Conversion form Cardan Angles to Matrices, Vel. and Acc. --- */

/* == Wtocardan =========================================================== */
/* function that extracts the Euler/Cardan angles and their first time
   derivative from position matrix m and velocity matrix W respectively
    Input parameters
      m       : position matrix
      W       : velocity matrix
      i,j,k   : constants that define the rotation sequence
		(they must be X,Y or Z)
    Output parameters
      q1,q2   : arrays of 3 real that contain the Euler/Cardan angles
      qp1,qp2 : arrays of 3 real that contain the Euler/Cardan angles first
		time derivative

    Note : The first time derivative of Euler/Cardan angles is evaluated using
	   the relation
			  qpx=omega*A(i)   (i)=invers
	   where qpx can be either qp1 or qp2
*/

int Wtocardan(MAT4 m, MAT4 W, int i, int j, int k, real q1[3], real q2[3],
	      real qp1[3], real qp2[3])

#define alpha1 q1[0]
#define beta1  q1[1]
#define alpha2 q2[0]
#define beta2  q2[1]
{

	MAT3 OMEGA,Ai;
	VECTOR omega;
	int sig,test;

	if (i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (Wtocardan) : illegal axis\n");
		return(NOTOK);
	}
	if ((j-i+3)%3==1)               /* cyclic */
		sig=1;
	else                            /* anti-cyclic */
		sig=-1;

	test=Mtocardan(m,i,j,k,q1,q2);  /* Eul./Card. angles from m */
	if (test==OK)
	{
		mcopy43(W,OMEGA);	/* a. v. matrix */
		mtov3(omega,OMEGA);	/* a. v. vector */

		invA(alpha1,beta1,sig,i,j,k,Ai);
		moltmv3(Ai,omega,qp1);  /* angles first time derivative */

		invA(alpha2,beta2,sig,i,j,k,Ai);
		moltmv3(Ai,omega,qp2)   /* angles first time derivative */
	}
	else
	{
		fprintf(stderr,"Error (Wtocardan) : there is a singularity\n");
		return(NOTOK);
	}

#undef alpha1
#undef beta1
#undef alpha2
#undef beta2

	return(OK);
}

/* == Htocardan =========================================================== */
/* function that extracts the Euler/Cardan angles and their first and second
   time derivative from position matrix m velocity matrix W and acceleration
   matrix H respectively
    Input parameters
      m         : position matrix
      W         : velocity matrix
      H         : acceleration matrix
      i,j,k     : constants that define the rotation sequence
		  (they must be X,Y or Z)
    Output parameters
      q1,q2     : arrays of 3 real that contain the Euler/Cardan angles
      qp1,qp2   : arrays of 3 real that contain the Euler/Cardan angles first
		  time derivative
      qpp1,qpp2 : arrays of 3 real that contain the Euler/Cardan angles
		  second time derivative

    Note : The second time derivative of the Euler/Cardan angles is evaluated
	   using the relation
			  qppx=Ainverse*(omegapto-Atilde*qtilde)
	   where qppx is qpp1 or qpp2 and qtilde is [bp*gp,ap*gp,ap*bp]
*/

int Htocardan(MAT4 m, MAT4 W, MAT4 H, int i, int j, int k, real q1[3],
	      real q2[3], real qp1[3], real qp2[3], real qpp1[3],
	      real qpp2[3])

#define alpha1 q1[0]
#define beta1  q1[1]
#define alpha2 q2[0]
#define beta2  q2[1]
{

	MAT3 G,Ai,Atil;
	VECTOR omegapto,qtilde,prod,aux;
	int sig,test;

	if (i<X || i>Z || j<X || j>Z || k<X || k>Z || i==j || j==k)
	{
		fprintf(stderr,"Error (Htocardan) : illegal axis\n");
		return(NOTOK);
	}
	if ((j-i+3)%3==1)       	/* cyclic */
		sig=1;
	else
		sig=-1;			/* anti-cyclic */

	/* Euler/Cardan angles and first time derivative */
	test=Wtocardan(m,W,i,j,k,q1,q2,qp1,qp2);
	if (test==OK)
	{
		mcopy43(H,G);
		Gtomegapto(G,omegapto);

		/* evaluates qpp1 */
		WPRODtocardan(alpha1,beta1,sig,i,j,k,Atil);
		qtilde[i]=qp1[Y]*qp1[Z]; qtilde[j]=qp1[X]*qp1[Z];
		qtilde[3-i-j]=qp1[X]*qp1[Y];

		moltmv3(Atil,qtilde,prod);
		aux[i]=prod[X]; aux[j]=prod[Y];
		aux[3-i-j]=prod[Z];
		subv(omegapto,aux,aux);

		invA(alpha1,beta1,sig,i,j,k,Ai);
		moltmv3(Ai,aux,qpp1);

		/* evaluates qpp2 */
		WPRODtocardan(alpha2,beta2,sig,i,j,k,Atil);
		qtilde[i]=qp2[Y]*qp2[Z]; qtilde[j]=qp2[X]*qp2[Z];
		qtilde[3-i-j]=qp2[X]*qp2[Y];

		moltmv3(Atil,qtilde,prod);
		aux[i]=prod[X]; aux[j]=prod[Y];
		aux[3-i-j]=prod[Z];
		subv(omegapto,aux,aux);

		invA(alpha2,beta2,sig,i,j,k,Ai);
		moltmv3(Ai,aux,qpp2);
	}
	else
	{
		fprintf(stderr,"Error (Htocardan) : there is a singularity\n");
		return(NOTOK);
	}

#undef alpha1
#undef beta1
#undef alpha2
#undef beta2

	return(OK);
}

/* == invA ================================================================ */
/* function that builds the inverse of the matrix A. It is useful in order to
   evaluate the first time derivative of the Euler/Cardan angles
    Input parameters
      alpha, beta : the first two Euler/Cardan angles
      sig         : parameter that defines the sign of some elements in the
		    inverse of A
      i,j,k       : constants that define the rotation sequence
		    (they must be X,Y or Z)
    Output parameters
      Ai          : it's the matrix where the inverse of A is stored
*/

int invA(real alpha, real beta, int sig, int i, int j, int k, MAT3 Ai)
{
	int l;
	real sa,ca,sb,cb;     /* sine and cosine of two Euler/Cardan angles */

	sa = (real) sin(alpha); ca = (real) cos(alpha);
	sb = (real) sin(beta);  cb = (real) cos(beta);

	if ( ((i!=k) && (cb==0.)) || ((i==k) && (sb==0.)) ) /* singularity */
		return (NOTOK);
	if (i!=k)	      /* Cardan convention */
	{
		Ai[X][i]=1.; Ai[X][j]=sa*sb/cb;     Ai[X][k]=(-sig)*ca*sb/cb;
		Ai[Y][i]=0.; Ai[Y][j]=ca;           Ai[Y][k]=sig*sa;
		Ai[Z][i]=0.; Ai[Z][j]=(-sig)*sa/cb; Ai[Z][k]=ca/cb;
	}
	else	              /* Euler convention */
	{
		l=3-i-j;

		Ai[X][i]=1.; Ai[X][j]=-sa*cb/sb;    Ai[X][l]=sig*ca*cb/sb;
		Ai[Y][i]=0.; Ai[Y][j]=ca;           Ai[Y][l]=sig*sa;
		Ai[Z][i]=0.; Ai[Z][j]=sa/sb;	    Ai[Z][l]=(-sig)*ca/sb;
	}
	return(OK);
}

/* == WPRODtocardan ======================================================= */
/* function that builds the matrix A tilde useful in order to evaluate the
   second time derivative of the Euler/Cardan angles
    Input parameters
      alpha, beta : the first two Euler/Cardan angles
      sig         : parameter that defines the sign of some elements in
		    matrix A tilde
      i,j,k       : constants that define the rotation sequence
		    (they must be X,Y or Z)
    Output parameters
      Atil        : it's the matrix where A tilde is stored
*/

void WPRODtocardan(real alpha, real beta, int sig, int i, int j, int k,
		   MAT3 Atil)
{
	int l;
	real sa,ca,sb,cb;     /* sine and cosine of two Euler/Cardan angles */

	sa = (real) sin(alpha); ca = (real) cos(alpha);
	sb = (real) sin(beta);  cb = (real) cos(beta);

	if (i!=k)	      /* Cardan convention */
	{
		Atil[X][i]=sig*cb;    Atil[X][j]=0.;
		Atil[Y][i]=sig*sa*sb; Atil[Y][j]=(-sig)*ca*cb;
		Atil[Z][i]=-ca*sb;    Atil[Z][j]=-sa*cb;

		Atil[X][k]=0.;
		Atil[Y][k]=-sa;
		Atil[Z][k]=sig*ca;
	}
	else		      /* Euler convention */
	{
		l=3-i-j;
		Atil[X][i]=-sb;          Atil[X][j]=0.;
		Atil[Y][i]=sa*cb;        Atil[Y][j]=ca*sb;
		Atil[Z][i]=(-sig)*ca*cb; Atil[Z][j]=sig*sa*sb;

		Atil[X][l]=0.;
		Atil[Y][l]=-sa;
		Atil[Z][l]=sig*ca;
	}
}


/* --- Working with Points Lines and Planes, Operations on Points --- */

/* == intermediate ======================================================== */
/* function that evalutes the middle point between two given points
    Input parameters
      P1 : first point
      a  : weight of P1
      P2 : second point
      b  : weight of P2
    Output parameters
      P3 : middle point

    Note : We have P3=([P1]*a+[P2]*b)/sqrt(a*a+b*b). If a==b==1 P3 is
	   exactly the middle point between P1 and P2
*/

void intermediate(POINT P1, real a, POINT P2, real b, POINT P3)
{
	int i;
	real c;

	c=a+b;
	if (c==0)
	{
		fprintf(stderr,"Error (intermediate) : c==0\n");
		/*c=sqrt(a*a+b*b);
		if (c==0)*/                                         /* by dino 2005 */
			c=1.;             /* if c==0 it is assumed c = 1 */
	}
	for (i=X;i<=Z;i++)
		P3[i]=(P1[i]*a+P2[i]*b)/c;
	P3[U]=1.;
}


/* --- Working with Points Lines and Planes, Op. on Lines and Planes --- */

/* == line2p ============================================================== */
/* function that builds a line passing through point P1 and P2. The unit
   vector of l has the same direction of vector P1->P2. The origin of line l
   is P1
    Input parameters
      P1,P2 : points that define the direction of line
    Output parameters
      l     : pointer to LINE structure which is filled with the
	      resulting values
*/

void line2p(POINT P1, POINT P2, LINEP l)
{
	pcopy(P1,l->P);
	vect(P2,P1,l->dir);
	unitv(l->dir,l->dir); /* evaluates the unit vector of line */
}

/* == linepvect =========================================================== */
/* function that builds a line passing through point P1. The direction of its
   unit vector is given by vector v
    Input parameters
      P1 : point that defines the origin of line
      v  : vector that defines the unit vector associated with the
	   direction of l
    Output parameters
      l  : pointer to LINE structure which is filled with the resulting values
*/

void linepvect(POINT P1, VECTOR v, LINEP l)
{
	pcopy(P1,l->P);
	unitv(v,l->dir);     /* unit vector of line */
}

/* == projponl ============================================================ */
/* function that finds the projection I of the point P1 on the line l
    Input parameters
      P1  : point that has to be projected on line
      l   : line that contains the projection point I
     Output parameters
      I   : projection point of the given point P1

     Note : Returns the distance between point P1 and line l
*/

real projponl(POINT P1, LINE l, POINT I)
{
	real r;
	PLANE pl;            /* plane trough P1 and orthogonal to line l */

	plane2(P1,l.dir,pl); /* builds plane pl */
	r=project(l.P,pl,I); /* finds I */
	r=dist( P1, I);      /* this line was added by Br1 */
	return r;
}

/* == distpp ============================================================== */
/* function that evaluates the distance (with a sign) of point P from plane pl
    Input parameters
      pl : plane whose distance from P must be evaluated
      P  : point whose distance from pl must be evaluated

    Note : Returns the distance between point P and plane pl
*/

real distpp(PLANE pl, POINT P)
{
	return dot2(pl,P,4);
}

/* == project ============================================================= */
/* function that finds the projection I of the point P on the plane pl
    Input parameters
      P  : point that has to be projected on plane
      pl : plane that contains the projection point I
    Output parameters
      I  : projection point of the given point P

    Note : Returns the distance between point P and plane pl
*/

real project(POINT P, PLANE pl, POINT I)
{
	real r;

	I[U]=1;
	r=distpp(pl,P);
	rmoltv(pl,r,I);
	subv(P,I,I);
	return r;
}

/* == plane =============================================================== */
/* function that builds the plane pl containing three points P1,P2,P3
    Input parameters
      P1,P2,P3 : points that lie on the plane
    Output parameters
      pl       : plane that contains the three points
*/

void plane(POINT P1, POINT P2, POINT P3, PLANE pl)
{
	VECTOR v12,v32;

	vect(P1,P2,v12); vect(P3,P2,v32);
	cross(v12,v32,pl); unitv(pl,pl);
	pl[U]=-dot(P2,pl);
}

/* == plane2 ============================================================== */
/* function that builds a plane pl which contains a point P1 and whose unit
   vector is v */

void plane2(POINT P1, VECTOR v, PLANE pl)
{
	unitv(v,pl);
	pl[U] = -dot(P1,pl); /* distance from the origin */
}

/* == inters2pl =========================================================== */
/* function that finds the intersection of two planes
    Input parameters
     pl1 : first plane
     pl2 : second plane
    Output parameters
     l   : line structure pointer where the intersection is stored

    Note : Returns NOTOK when pl1 is parallel to pl2 (no intersection)
*/

int inters2pl(PLANE pl1, PLANE pl2, LINEP l)
{
	real zero,d;
	real alpha,beta;
	VECTOR par;
	int i;

	zero=zerom();        /* machine precision */
	cross(pl1,pl2,par);
			     /* planes are not parallel */
	if (aabs(par[X]) > zero || aabs(par[Y]) > zero || aabs(par[Z]) > zero)
	{
		l->dir[X]=par[X]; l->dir[Y]=par[Y]; l->dir[Z]=par[Z];
		d=dot(pl1,pl2);
		alpha=(-pl1[U]+d*pl2[U])/(1-pow(d,2));
		beta=(d*pl1[U]-pl2[U])/(1-pow(d,2));
		for (i=X;i<U;i++)
			l->P[i]=alpha*pl1[i]+beta*pl2[i];
		l->P[U]=1;
		return (OK);
	}
	else                 /* planes are parallel */
		return(NOTOK);
}

/* == interslpl =========================================================== */
/* function that finds the intersection of a line l with a plane pl
    Input parameters
      l       : line whose intersection must be found
      pl      : plane whose intersection must be found
    Output parameters
      I       : point where the intersection is stored
      inttype : defines which kind of intersection was found
		 1 = the line is coincident with the plane (the intersection
		     is the line itself)
		-1 = the line is parallel to the plane (no intersection)
		 0 = there is only one intersection
*/

void interslpl(LINE l, PLANE pl, POINT I, int *inttype)
{
	real par,lie,zero,alpha;
	VECTOR v;

	zero=zerom();
	par=dot(pl,l.dir);
	lie=dot2(pl,l.P,4);
	if(fabs(par)<zero)	             /* l is parallel to pl */ /* fabs() bug fixed jan 2003 */
	{
		if(lie<zero)         /* pl passes through l.P */
			*inttype=1;  /* the intersection is the line itself */
		else		     /* pl does not pass through point l.P  */
			*inttype=-1; /* no intersection */
	}
	else			     /* line and plane are incident */
	{
		*inttype=0;	     /* there is only one intersection */
		alpha=(-lie) / par;
		vector(l.dir,alpha,v);
		sumv(v,l.P,I);       /* I=alpha*l.dir+l.P */
		I[U]=1.;
	}
}

/* == intersection ======================================================== */
/* function that finds the intersection of two lines. It builds also
   (if possible) the plane which contains the two lines. Builds the minimum
   distance line and evaluates the distance between the lines
    Input parameters
      l1,l2    : line structures whose intersection must be evaluated
    Output parameters
      lmindist : is a line orthogonal to l1 and l2 which has the minimum
		 distance from l1 and l2
      mindist  : the module of the minimum distance between l1 and l2
      pl       : plane that contains l1 and l2
      I        : intersection point of l1 and l2
      inttype  : defines which kind of intersection was found
		  1 = point I is not real intersection, only middle point
		      between two oblique lines. Plane pl does not contain
		      really l1 and l2
		  0 = point I is the intersection of l1 and l2. Plane pl
		      contains the two lines
		 -1 = there are infinite intersections (l1==l2). Plane pl is
		      undefined
		  2 = there aren't any intersections because l1 and l2 are
		      parallel . Plane pl can be found
*/

void intersection(LINE l1, LINE l2, LINEP lmindist, real * mindist, PLANE pl,
		  POINT I, int *inttype)
{
	VECTOR deltap,v1,v2,u3;
	POINT A,B;
	real a,c1,c2,alpha,beta,zero;

	zero = zerom();          /* machine zero */
	vect(l2.P,l1.P,deltap);  /* deltap=(P2-P1) */
	cross(l1.dir,l2.dir,u3); /* u3 is orthogonal to l1 and l2 */

				 /* l1 and l2 are not parallel */
	if (aabs(u3[X]) > zero || aabs(u3[Y]) > zero || aabs(u3[Z]) > zero)
	{
		c1 = dot(l1.dir,deltap); /* c1=U1t*(P2-P1) */
		c2 = dot(l2.dir,deltap); /* c2=U2t*(P2-P1) */
		a = dot(l1.dir,l2.dir);  /* a=U1t*U2 */
		alpha = (real) ( (a*c2-c1) / (pow(a,2)-1) );
		beta = (real) ( (c2-a*c1) / (pow(a,2)-1) );
		vector(l1.dir,alpha,v1); /* builds vector l1.P->A */
		vector(l2.dir,beta,v2);  /* builds vector l2.P->B */
		sumv(l1.P,v1,A);         /* point A in the reference frame */
		sumv(l2.P,v2,B);         /* point B in the reference frame */
		A[U]=B[U]=1;
		*mindist = dist(A,B);    /* minimum distance between l1, l2 */
		middle(A,B,I);           /* middle point between A and B */
		pcopy(I,lmindist->P);    /* origin of minimum distance line */
		vcopy(u3,lmindist->dir); /* u3 is the dir. of this line */
		unitv(lmindist->dir,lmindist->dir);
		plane2(I,u3,pl);         /* builds the minimum distance pl. */
		if ( (*mindist) > zero ) /* l1 and l2 are oblique lines */
			*inttype=1;      /* intersesction is middle point */
		else			 /* l1 and l2 are incident lines */
		{
			*mindist=0.;
			*inttype=0;      /* I is the real intersection */
		}
	}
	else                              /* l1 is parallel to l2 */
	{
		cross(deltap,l1.dir,v1);
		*mindist = mod(v1);	  /* min. dist. between l1, l2 */
		if ( (*mindist) < zero)   /* l1==l2 */
			*inttype = -1;
		else 			  /* l1!=l2 */
		{
			*inttype = 2; 	  /* no intersection between l1, l2 */
			cross(deltap,l1.dir,v1); /* dir. of min. dist. line */
			cross(v1,l1.dir,v2);
			unitv(v2,lmindist->dir);
			pcopy(l1.P,lmindist->P);
			plane2(l1.P,v1,pl); 	 /* pl. containing l1, l2 */
		}
	}
}


/* --- Operations on Matrices and Vectors, Gen. Operations on Matrices --- */

/* == pseudo_inv ========================================================= */
/* function that builds the pseudo-inverse matrix Api of a given matrix A. It
   performs Greville's algorithm by rows
    Input parameters
      real A[][]     : initial matrix whose pseudo-inverse must be found
      int  rows,cols : dimensions of matrix A
    Output parameters
      Api[][]        : matrix where the pseudo-inverse is stored

    Note : This algorithm can build also the pseudo-inverse of a rectangular
	   matrix (when rows!=cols). The algorithm's main variables are
	   dk,ck,bk. They are used to build the pseudo-inverse matrix Api(k)
	   at each step k of the algorithm (k=1..cols). aux is used for
	   general purposes
*/

int pseudo_inv(MAT A, MAT Api, int rows, int cols)

#define a(r,c)    (*((A)+(cols)*(r)+(c)))
{

	int i,k,test;
	real *dk,*ck,*bk,*aux;
	real temp;
	real toll=1.0e-5;
	int irank=0;	/* by joe July 1998 */

	toll=norm(M A,rows,cols)*zerom()*1000.;	/* by joe october 2000 */
	ck = (real *) calloc(rows,sizeof(real));
	bk = (real *) calloc(rows,sizeof(real));
	clear(M Api,cols,rows);
	for (i=0,test=0;i<rows && test==0;i++) /* algorithm initialization */
		if (aabs(a(i,0)) > toll)       /* tests if a1!=0 */
			test=1;
	if (test)                              /* a1!=0 : Api(1)=a1(t)/(a1(t)*a1) */
	{
		irank++;
		aux = (real *) calloc(rows,sizeof(real));
		mvcopy(M A,rows,cols,1,Col,M aux);
		temp = 1/(dot2(M aux,M aux,rows));
		rmolt(M aux,temp,M aux,rows,1);
		vmcopy(M aux,rows,1,Row,M Api,cols,rows);
		free(aux);
	}
	for (k=2;k<=cols;k++)                  /* algorithm body */
	{
					       /* dim. of dk depends on k */
		dk = (real *) calloc(k-1,sizeof(real));
		clear(M dk,k-1,1);
		clear(M ck,rows,1);
		clear(M bk,rows,1);
					       /* dk=Api(k-1)*ak */
		aux = (real *) calloc(rows,sizeof(real));
		mvcopy(M A,rows,cols,k,1,M aux);
		molt(M Api,M aux,M dk,k-1,rows,1);
		free(aux);
					       /* ck=ak-A(k-1)*dk */
		aux = (real *) calloc((rows*(k-1)),sizeof(real));
		clear(M aux,rows,k-1);
		mmcopy(M A,M aux,cols,k-1,rows,k-1);
		molt(M aux,M dk,M ck,rows,k-1,1);
		free (aux);
		aux = (real *) calloc(rows,sizeof(real));
		mvcopy(M A,rows,cols,k,Col,M aux);
		sub(M aux,M ck,M ck,1,rows);
		free(aux);
					       /* tests if ck!=0 */
		for (i=0,test=0;i<rows && test==0;i++)
			if (aabs(*(ck+i)) > toll)
				test=1;
		if (test) /* ck!=0 : bk(t)=ck(t)/(ck(t)*ck) */
		{
			irank ++;
			temp = dot2(M ck,M ck,rows);
			temp = 1/temp;
			mcopy(M ck,M bk,rows,1);
			rmolt(M bk,temp,M bk,rows,1);
		}
		else      /* ck==0 : bk(t)=(dk(t)*Api(k-1))/(1+dk(t)*dk) */
		{
			aux = (real *) calloc(k,sizeof(real));
			mcopy(M dk,M aux,k-1,1);
			molt(M aux,M Api,M bk,1,k-1,rows);
			temp = dot2(M dk,M dk,k-1);
			temp = 1/ (1+temp);
			rmolt(M bk,temp,M bk,rows,1);
			free(aux);
		}
					       /* builds Api(k) */
		aux = (real *) calloc(((k-1)*rows),sizeof(real));
		clear(M aux,k-1,rows);
		molt(M dk,M bk,M aux,k-1,1,rows);
		sub(M Api,M aux,M Api,k-1,rows);
		vmcopy(M bk,rows,k,Row,M Api,cols,rows);
		free(aux);
		free(dk);
	}
	free(bk); free(ck);

	return irank;

#undef a
}

/* == crossmtom =========================================================== */
/* function that evaluates the product C=bt*a-at*a where at and bt are the
   transpose of a and b respectively */

void crossmtom(real *a, real *b, MAT C, int dim)

#define c(l,m)    (*((C)+(dim)*(l)+(m)))
{
	int i,j;

	for (i=0;i<dim;i++)
		for (j=0;j<dim;j++)
			c(i,j)=b[i]*a[j]-a[i]*b[j];

#undef c
}


/* --- Operations on Matrices and Vectors, Gen. Operations on Vectors --- */

/* == dot2 ================================================================ */
/* function that evaluates the dot product of two dim elements vectors
    Input parameters
      v1,v2 : two dim elements vectors
      dim   : dimension of the two vectors

    Note : This function can be considered an extension of function dot
	   because it is not limited to the dot product of 3-elements vectors
	   but can evaluate dot product of two any-dimension vectors
*/

real dot2(real *v1, real *v2, int dim)
{
	int i;
	real aux;

	for (i=0,aux=0;i<dim;i++)
		aux += v1[i]*v2[i];
	return(aux);
}


/* --- Copy Functions --- */

/* == mvcopy ============================================================== */
/* function that extracts a row or a column from a matrix A. A is a rows*cols
   matrix while v is the vector where column or row is stored.
    Input parameters
      A         : matrix where the selected row or column is initially stored
      rows,cols : dimensions of matrix A
      val       : defines which row or column must be extracted
		  (range from 1 to rows/cols)
      type      : defines which kind of element must be extracted (Row or Col)
    Output parameters
      v         : the vector where the extracted row/column is stored
*/

int mvcopy(MAT A, int rows, int cols, int val, int type, real *v)

#define a(i,j)   (*((A)+(cols)*(i)+(j)))
#define v(i)	 (*(v+i))
{
	int i;

	if (type<Row || type>Col || val<1 || (type==Col && val>cols) || (type==Row && val>rows))
	{
		fprintf(stderr,"Error (mvcopy) : illegal value of parameters\n");
		return(NOTOK);
	}
	if (type==Col)        /* extracts a column */
		for (i=0;i<rows;i++)
			v(i) = a(i,val-1);
	else                  /* extracts a row */
		for (i=0;i<cols;i++)
			v(i) = a(val-1,i);

#undef a
#undef v

	return(OK);
}

/* == vmcopy ============================================================== */
/* function that copy a vector into a row or a column of a matrix A. A is a
   rows*cols matrix while v is the vector where column or row is initially
   stored.
    Input parameters
      v         : the vector where the row/column is initially stored
      rows,cols : dimensions of matrix A
      val       : defines into which row or column v must be copied
		  (range from 1 to rows/cols)
      type      : defines which kind of element must be extracted (Row or Col)
    Output parameters
      A         : matrix where the selected row or column is finally stored
*/

int vmcopy(real *v, int dim, int val, int type, MAT A, int rows, int cols)

#define a(i,j)   (*((A)+(cols)*(i)+(j)))
#define v(i)	 (*(v+i))
{
	int i;

	if (type<Row || type>Col || val<1 || (type==Col && val>cols) ||
	   (type==Row && val>rows))
	{
		fprintf(stderr,"Error (vmcopy) : illegal value of parameters\n");
		return(NOTOK);
	}
	if (type==Col)        /* copy a column */
		for (i=0;i<dim;i++)
			a(i,val-1)=v(i);
	else                  /* copy a row */
		for (i=0;i<dim;i++)
			a(val-1,i)=v(i);

#undef a
#undef v

	return(OK);
}


/* --- Print Functions --- */

/* == printmat ============================================================ */
/* function that prints on the standard output a imax*jmax submatrix of a
   idim*jdim matrix A preceded by the string s */

void printmat(char *s, MAT A, int idim, int jdim, int imax, int jmax)

#define a(ii,jj) (*(A+(jj)+jdim*(ii)))  /* attention to the define ! */
{
	int i,j,dummy;
	real ele;

	dummy=idim;           /* dummy statement */
	idim=dummy;           /* dummy statement */
	printf("%s\n",s);
	for (i=0;i<imax;i++)
	{
		for (j=0;j<jmax;j++)
		{
			ele=a(i,j);
			printf("%8.3f ",ele);
		}
		printf("\n");
	}

#undef a
}

/* == iprintmat =========================================================== */
/* function that prints on the standard output a imax*jmax submatrix of a
   idim*jdim matrix A preceded by the string s */

void iprintmat(char *s, int *A, int idim, int jdim, int imax, int jmax)

#define a(ii,jj) (*(A+(jj)+jdim*(ii)))  /* attention to the define ! */
{
	int i,j,dummy;
	int ele;

	dummy=idim;           /* dummy statement */
	idim=dummy;           /* dummy statement */
	printf("%s\n",s);
	for (i=0;i<imax;i++)
	{
		for (j=0;j<jmax;j++)
		{
			ele=a(i,j);
			printf("%8d ",ele);
		}
		printf("\n");
	}

#undef a
}
