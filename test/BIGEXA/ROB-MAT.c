/* ROB-MAT: program for direct kinematics and inverse dynamics of ANY serial robot v.2
   Developed on MS-DOS operative system with Microsoft C compiler V. 5.10 */

/* Note: To compile this program the type real must be set equivalent to the type float
	 (see also User's Manual). This is necessary because the formatting string of the
	 fscanf function have been written using %f as descriptor. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPACELIB.H"

int main(int argc,char *argv[])
{
	#define MAXLINK 10                                /* maximum number of links */

	int nlink,jtype[MAXLINK];                         /* number of links;  joint type */
	real theta[MAXLINK],d[MAXLINK];  /* Extended D.&H. parameters */
	real b[MAXLINK],a[MAXLINK],alpha[MAXLINK];
	real m,jxx,jxy,jxz,jyy,jyz,jzz,xg,yg,zg;/* dynamics parameters */
	real q,qp,qpp;                   /* joint variables */
	real gx,gy,gz;                   /* gravity acceleration */
	real fx,fy,fz,cx,cy,cz;    	 /* external forces and torques on end-effector */

	MAT4 *mrel,*T,*W,*W0,*WA,*H,*H0,*HA,*J,*J0,*FI,*ACT0;/* matrices */
	MAT4 EXT,Hg,Ht;
	MAT4 TMP;                        /* temporary matrix */
	FILE *data;                      /* file including robot parameters */
	FILE *motion;                    /* file including actuator motions */
	int i;                           /* counter */
	int ierr;                        /* error code */
	int size;                        /* for use of calloc function */
	float t,dt;

	if(argc!=3){
		printf("Usage: Rob-Mat  data_file  motion_file\n");
		exit(1);
	}
	data=fopen(argv[1],"r");         /* check of input data */
	if(data==NULL)
		exit(2);
	motion=fopen(argv[2],"r");
	if(motion==NULL)
		exit(3);
	size=sizeof(MAT4);               /* dynamic allocation of memory */
	mrel =(MAT4 *) calloc(MAXLINK,size);
	T =(MAT4 *) calloc(MAXLINK,size);
	W =(MAT4 *) calloc(MAXLINK,size);
	W0=(MAT4 *) calloc(MAXLINK,size);
	WA=(MAT4 *) calloc(MAXLINK,size);
	H =(MAT4 *) calloc(MAXLINK,size);
	H0=(MAT4 *) calloc(MAXLINK,size);
	HA=(MAT4 *) calloc(MAXLINK,size);
	J =(MAT4 *) calloc(MAXLINK,size);
	J0=(MAT4 *) calloc(MAXLINK,size);
	FI=(MAT4 *) calloc(MAXLINK,size);
	ACT0=(MAT4 *) calloc((MAXLINK+1),size);            /* step (1) */

	idmat4(T[0]);                    /* INITIALIZATION of matrices */
	clear4(WA[0]);
	clear4(HA[0]);
					 /* read robot description */
	fscanf(data,"%d",&nlink);        /* n. of links */
	for (i=1;i<=nlink;i++) 		 /* for each link */
	{                                /* D.&H. parameters */
					 /* dynamic data */
		fscanf(data,"%d %f %f %f %f %f",
		       &jtype[i],&theta[i],&d[i],&b[i],&a[i],&alpha[i]);
		fscanf(data,"%f %f %f %f %f %f %f",
		       &m,&jxx,&jxy,&jxz,&jyy,&jyz,&jzz);
		fscanf(data,"%f %f %f",&xg,&yg,&zg);
					   /* build inertia matrix */
		jtoJ(m,jxx,jyy,jzz,jxy,jyz,jxz,xg,yg,zg,J[i]);
	}
				       /* read gravity acceleration vector */
	fscanf(data,"%f %f %f",&gx,&gy,&gz);
	gtom(gx,gy,gz,Hg);             /* build gravity acceleration matrix */
	fscanf(motion,"%f",&dt);       /* read the range of time */
	for(t=0;;t+=dt)                /* for each instant of time */
	{
					   /* ***** KINEMATICS *****  */
		for (i=1;i<=nlink;i++) /* for each link */
		{
					   /* read motions (2) */
			ierr=fscanf(motion,"%f %f %f",&q,&qp,&qpp);
			if (ierr!=3)   /* end of data in file MOTION */
			{       fcloseall();
				exit(0);
			}
					  /* build relative position matrix (3) */
			dhtom(jtype[i],theta[i],d[i],b[i],a[i],alpha[i],q,
				  mrel[i]);/* build relative velocity and acceleration matrix in local frame (4) */
			velacctoWH(jtype[i],qp,qpp,W[i],H[i]);
			molt4(T[i-1],mrel[i],T[i]);   /* evaluate absolute position matrix (5) */
			trasf_mami(W[i],T[i-1],W0[i]);/* transform relative velocity matrix from local frame to base frame  (6) */
			trasf_mami(H[i],T[i-1],H0[i]);/* transform relative acceleration matrix from local frame to base frame (7) */
			sum4(WA[i-1],W0[i],WA[i]);    /* evaluate absolute velocity matrix  (8) */
				   /* evaluate absolute acceleration matrix (9) */
			coriolis(HA[i-1],H0[i],WA[i-1],W0[i],HA[i]);
		}
				   /* ***** DYNAMICS ***** */
				   /* initializations  (10) */
				   /* read external actions on end-effector */
		fscanf(data,"%f %f %f %f %f %f",&fx,&fy,&fz,&cx,&cy,&cz);
		actom(fx,fy,fz,cx,cy,cz,EXT);           /* build external action matrix */
		trasf_mamt4(EXT,T[nlink],ACT0[nlink+1]);/* transforms external actions from local to base frame */
		for(i=nlink;i>0;i--)                    /* for each link */
		{
			trasf_mamt4(J[i],T[i],J0[i]);   /* transform inertia matrix from local to base frame  (11) */
			rmolt4(HA[i],-1.,TMP);          /* change sign to find inertia action */
			sum4(TMP,Hg,Ht);                /* evaluate total acceleration matrix */
			skew4(Ht,J0[i],FI[i]);          /* evaluate the action matrix due to inertia and weight (12) */
			sum4(FI[i],ACT0[i+1],ACT0[i]);  /* evaluate total action matrix (13) */
		}
					      /* ***** OUTPUT RESULTS ***** */
		for(i=1;i<=nlink;i++)         /* for each link */
		{
			printf("\n\n Link %d \n\n",i);
			printm4("Rel. position matrix",mrel[i]);
			printm4("Absolute position matrix",T[i]);
			printm4("Rel. velocity matrix in frame (i)",W[i]);
			printm4("Rel. velocity matrix in frame (0)",W0[i]);
			printm4("Absolute velocity matrix in frame (0)",
				WA[i]);
			printm4("Rel. acceleration matrix in frame (i)",H[i]);
			printm4("Rel. acceleration matrix in frame (0)",H0[i]);
			printm4("Absolute acceleration matrix in frame (0)",
				HA[i]);
			printm4("Inertia matrix in frame (i)",J[i]);
			printm4("Inertia matrix in frame (0)",J0[i]);
			printm4("Total actions",FI[i]);
			printm4("Action on joint i",ACT0[i]);
		}
	}
}  			                      /*  end main  */
