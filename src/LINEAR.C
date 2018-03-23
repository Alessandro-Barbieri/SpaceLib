/*
------------------------------------------------------------------------------
		 LINEAR.C (c)
		 copyright by

		 G. Legnani and R. Faglia

		 preliminary version v0.1 - Oct. 1990

		 update and bugs fixed July 1998

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
#include "SPACELIB.H"
#include "LINEAR.H"


/* --- Speed and Acceleration Matrices --- */

/* == dyn_eq ============================================================== */

int   dyn_eq(MAT4 J, MAT4 Wp, MAT4 F, int var[2][6])
{
	real mat[6][7];
	real acc[6],t,arm;
	int i,j,k,it,nvar,irank,neq;
	int var_a[6],var_f[6],ivar[6],tappo;
	int var2[2][6];

	nvar=0;
	for(i=0;i<2;i++)
		for(j=0;j<6;j++)
			if (var[i][j])
				nvar++;
	if (nvar!=6)
	{
		fprintf(stderr,"Error (dyn_eq) : \a wrong number of unknown variables\n");
		return(NOTOK);
	}

	for(i=0;i<6;i++)
	{
		acc[i]=0.; mat[i][6]=0.; ivar[i]=i;
	}
	if (!var[0][X])		acc[X]=Wp[Z][Y];          /* angular acc. */
	if (!var[0][Y])		acc[Y]=Wp[X][Z];
	if (!var[0][Z])		acc[Z]=Wp[Y][X];
	if (!var[0][3+X])		acc[3+X]=Wp[X][U];        /* accelerations */
	if (!var[0][3+Y])		acc[3+Y]=Wp[Y][U];
	if (!var[0][3+Z])		acc[3+Z]=Wp[Z][U];
	if (!var[1][X])		mat[X][6]=F[Z][Y];        /* torques */
	if (!var[1][Y])		mat[Y][6]=F[X][Z];
	if (!var[1][Z])		mat[Z][6]=F[Y][X];
	if (!var[1][3+X])		mat[3+X][6]=F[X][U];      /* forces */
	if (!var[1][3+Y])		mat[3+Y][6]=F[Y][U];
	if (!var[1][3+Z])		mat[3+Z][6]=F[Z][U];

	mat[X][X]=J[Y][Y]+J[Z][Z]; mat[Y][Y]=J[X][X]+J[Z][Z];
	mat[Z][Z]=J[X][X]+J[Y][Y];

	mat[X][Y]=mat[Y][X]= -J[X][Y];
	mat[X][Z]=mat[Z][X]= -J[X][Z];
	mat[Y][Z]=mat[Z][Y]= -J[Y][Z];

	mat[X][3+X]=0.;mat[Y][3+Y]=0.;mat[Z][3+Z]=0.;
	mat[3+X][X]=0.;mat[3+Y][Y]=0.;mat[3+Z][Z]=0.;

	mat[Y][3+Z]=mat[3+Z][Y]= -(mat[Z][3+Y]=mat[3+Y][Z]=J[X][U]);
	mat[Z][3+X]=mat[3+X][Z]= -(mat[X][3+Z]=mat[3+Z][X]=J[Y][U]);
	mat[X][3+Y]=mat[3+Y][X]= -(mat[Y][3+X]=mat[3+X][Y]=J[Z][U]);
	
	mat[3+X][3+X]=J[U][U]; mat[3+X][3+Y]=0.;      mat[3+X][3+Z]=0.;
	mat[3+Y][3+X]=0.     ; mat[3+Y][3+Y]=J[U][U]; mat[3+Y][3+Z]=0.;
	mat[3+Z][3+X]=0.     ; mat[3+Z][3+Y]=0.;      mat[3+Z][3+Z]=J[U][U];


	//printm67("mat",mat);

	//getch();

	for(j=0;j<6;j++)
	{
		var_a[j]=j; var_f[j]=j;
	}
	
	for(i=0;i<2;i++)
		for(j=0;j<6;j++)
			var2[i][j]=var[i][j];

	neq=0;
	for (j=0;j<6;j++)
	{
		if(var[0][j]) {
			neq++;
			//printf("j: %d neq: %d\n",j,neq);
		} else {
			for(i=j+1;i<6;i++)
			{
				if(var[0][i])
				{
					neq++;
					//printf(" swabc i:%d j:%d neq:%d\n",i,j,neq); /* debug only */
					swabc(M mat,6,7,6,i,j,var_a);
					t=acc[j]; acc[j]=acc[i]; acc[i]=t;
					it=var[0][j]; var[0][j]=var[0][i]; var[0][i]=it;
					it=var[1][j]; var[1][j]=var[1][i]; var[1][i]=it;
					//printf(" swabr \n"); /* debug only */
					swabr(M mat,6,7,7,i,j);
					k=var_f[j]; var_f[j]=var_f[i]; var_f[i]=k;
					break;
				}
			}
		}
	}

	//printf("**neq: %d\n",neq);   getch();
	tappo= -1;
	for(i=0;i<6;i++)
	{
		for(j=neq,t=0;j<6;j++)
		{
			t += mat[i][j]*acc[j];
		}
		mat[i][neq]=mat[i][6]-t;
	}
	if(neq!=0)
	{
		linear(M mat,6,7,neq,neq,1,ivar,&irank,&arm,&tappo);
		//printm67("mat--",mat);
		//printf("irank: %d",irank);
		if (irank!=neq){
			fprintf(stderr,"\a DYN_EQ: error matrix is singular!\n");
			return(-irank);
			}
	}

	for (i=0;i<neq;i++)
		acc[var_a[ivar[i]]]=mat[i][neq];

	for (i=0;i<4;i++)
	{
		Wp[i][i]=Wp[U][i]=0.;
	}
	if(var2[0][X]) Wp[Y][Z]= -(Wp[Z][Y]=acc[X]);
	if(var2[0][Y]) Wp[Z][X]= -(Wp[X][Z]=acc[Y]);
	if(var2[0][Z]) Wp[X][Y]= -(Wp[Y][X]=acc[Z]);
	if(var2[0][3+X]) Wp[X][U]=acc[3+X]; 
	if(var2[0][3+Y]) Wp[Y][U]=acc[3+Y];
	if(var2[0][3+Z]) Wp[Z][U]=acc[3+Z];

	skew4(Wp,J,F);

	return(OK);
}


/* --- Functions solve minvers linear, linear --- */

/* == linear ============================================================== */
/* function that solves a linear system using the double pivoting algorithm
   (performed by rows and columns)
    Input parameters
      H          : square matrix where the coefficients are stored
      idim, jdim : matrix A physical dimensions (storage)
      imax, jmax : logical dimension of the matrix
      nsol       : number of right-hand side vectors (it is the number of
		   system that has to be handled at once)
    Output parameters
      ivet       : unknown terms vector (it gives useful informations in order
		   to reorder the system solution). If ivet(i)==k the value of
		   the k-th element of the unknown vector can be found at i-th
		   position
      irank      : it's the rank of matrix A
      arm        : greater element found during each step of the double
		   pivoting algorithm into the considered submatrix
    Additional input parameters
      vpr        : vector of the main variables (it is ended by -1)
		   es. vpr= 3 5 7 -1
*/


real   rmax(MAT, int, int, int, int, int, int, int *, int *);

void   linear(MAT H, int idim,int jdim, int imax, int jmax, int nsol,
	      int ivet[],int *irank, real *arm, int vpr[])

{
	static real eps=0.,toll=5;
	real rm,rmin;
	int i,j,k, im,jm, jmax1,j2max;

	if (eps == 0.)
		eps=toll*zerom();

	/* submatrix area where the pivoting step is performed */
	*irank=0;
	for(j=0;j<jmax;j++)
		ivet[j]=j;
	j2max=jmax+nsol;
	jmax1=jmax;
	for(i=0;vpr[i]!=-1;i++)
	{
		swabc(H,idim,jdim,imax,jmax1-1,vpr[i],ivet);
		jmax1=jmax1-1;
	}
	k=min(imax,jmax);
	for(i=0;i<k;i++)
	{
		rm=rmax(H,idim,jdim,i,imax,i,jmax1,&im,&jm);
		*arm=aabs(rm);
		if (i==0)
			rmin=eps*(*arm)*imax;
		if (rm == 0. || *arm < rmin)
			return;
		*irank=i+1;
		swabr(H,idim,jdim,j2max,i,im);
		swabc(H,idim,jdim,imax,i,jm,ivet);
		elimin(H,idim,jdim,imax,j2max,i);
		normalr(H,idim,jdim,j2max,i,rm);
	}
}

/* == swabr =============================================================== */

void   swabr(MAT H, int idim, int jdim, int j2max,int k, int im)

#define h(i,j) (*(H+(j)+jdim*(i)))  /* attention to the define ! */
{
	int j, dummy;
	real t;

	dummy=idim;                 /* dummy statement */
	idim=dummy;                 /* dummy statement */
	if(k==im)
		return;

	for(j=0;j<j2max;j++)
	{
		t=h(k,j);
		h(k,j)=h(im,j);
		h(im,j)=t;
	}

#undef h /* bug fixed 14/04/97 AmGhiPu */
}

/* == swabc =============================================================== */

void   swabc(MAT H, int idim, int jdim, int imax, int k, int jm, int vet[])

#define h(i,j) (*(H+(j)+jdim*(i)))  /* attention to the define ! */
{
	int i,it,dummy;
	real t;

	dummy=idim;                 /* dummy statement */
	idim=dummy;                 /* dummy statement */
	if (k==jm)
		return;

	for(i=0;i<imax;i++)
	{
		t=h(i,k);
		h(i,k)=h(i,jm);
		h(i,jm)=t;
	}
	it=vet[k];
	vet[k]=vet[jm];
	vet[jm]=it;

#undef h /* bug fixed 14/04/97 AmGhiPu */
}

/* == normalr ============================================================= */

void   normalr(MAT H, int idim, int jdim, int j2max, int k, real rm)

#define h(i,j) (*(H+(j)+jdim*(i)))  /* attention to the define ! */
{
	int j,dummy;

	dummy=idim;                 /* dummy statement */
	idim=dummy;                 /* dummy statement */
	for(j=0;j<j2max;j++)
		h(k,j)=h(k,j)/rm;

#undef h /* bug fixed 14/04/97 AmGhiPu */
}

/* == elimin ============================================================== */

void   elimin(MAT H, int idim, int jdim, int imax, int j2max, int k)

#define h(i,j) (*(H+(j)+jdim*(i)))  /* attention to the define ! */
{
	real r,t;
	int i,j, dummy ;

	dummy=idim;                 /* dummy statement */
	idim=dummy;                 /* dummy statement */
	r=h(k,k);
	for(i=0;i<imax;i++)
	{
		if(i!=k)
		{
			t=h(i,k)/r;
			for (j=k;j<j2max;j++)
				h(i,j)=h(i,j)-h(k,j)*t;
		}
	}

#undef h /* bug fixed 14/04/97 AmGhiPu */
}

/* == rmax ================================================================ */
/* this functions finds the greatest element (with sign) of a imax*jmax
   matrix stored in a idim*jdim array. The matrix is not scanned completely :
   only the elements whose coordinate are greater or equal to imin, jmin are
   considered. It returns the position of this element into the matrix by
   im, jm. If im==jm==0 an error has occurred */

real   rmax(MAT H, int idim, int jdim, int imin, int imax, int jmin, int jmax,
	    int *im, int *jm)

#define h(i,j) (*(H+(j)+jdim*(i)))  /* attention to the define ! */
{
	real rmaxx=0.;
	int i,j;

	*im=0; *jm=0;
	if (imin > imax || imax > idim)
		return(0.);
	if (jmin > jmax || jmax > jdim)
		return(0.);
	rmaxx=0.;
	for(i=imin;i<imax;i++)
		for(j=jmin;j<jmax;j++)
		{
			if (aabs(h(i,j)) > aabs(rmaxx))
			{
				rmaxx=h(i,j); *im=i; *jm=j;
			}
		}

#undef h /* bug fixed 14/04/97 AmGhiPu */

	return(rmaxx);
}
