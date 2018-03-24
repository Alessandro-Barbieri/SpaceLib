/* elb_I_DH.c program for the INVERSE kinematics of ELBOW robot.
   Frames assigned according to Denavit and Hartenberg conventions.
   the output of this program is compatible with the input  of elb_d_dh.c
   the input  of this program is compatible with tho output of elb_d_dh.c */

/* Note: To compile this program the type real must be set equivalent to the type float
	 (see also User's Manual). This is necessary because the formatting string of the
	 fscanf function have been written using %f as descriptor. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "SPACELIB.H"
#include "LINEAR.H"

#define MAXLINK 6

int main(int argc,char *argv[])
{
				       /* Denavit & Hartenberg's parameters */
	real theta[MAXLINK+1]={0.,0.,0.,0.,0.,0.,0.};
	real d[MAXLINK+1]={0.,0.,0.,0.,0.,0.,0.};
	real b[MAXLINK+1]={0.,0.,0.,0.,0.,0.,0.};
	real fi[MAXLINK+1]={0.,PIG_2,0.,0.,3*PIG_2,PIG_2,0.};
	real a[MAXLINK+1];

	VECTOR q1,                     /* Euler/Cardan angle configurations */
	       qp1,                    /* Euler/Cardan angle first time derivative */
	       qpp1,vel,acc;           /* Euler/Cardan angle second time derivative */

				       /* array of joint position variables */
	real q[MAXLINK];               /* joint angles */
	real qp[MAXLINK];              /* array of joint velocity variables */
	real qpp[MAXLINK];             /* array of joint acceleration variables */
	real ds[MAXLINK],              /* solution of the eq. J*dq=ds */
	     dq[MAXLINK],              /* solution of Newton/Raphson algorithm step */
	     buf[MAXLINK];
	real t,dt;
	real toll=0.0005;              /* precision of the solution */
	int maxiter=15;                /* maximum number of iterations in N-R algorithm */
	int ierr,i,k,p;
        int ii,jj,kk;                  /* Euler/Cardan convention (gripper orientation) */
	int rank;                      /* rank of linear system */
	real n;
	POINT O,orig=ORIGIN;
	real Jac[MAXLINK][MAXLINK];    /* Jacobian matrix */
	MAT4 mtar,                     /* target position matrix */
	     mrelp_1[MAXLINK+1],       /* array containing position matrix of frame
					  (p) seen in frame (p-1) */
	     mabs[MAXLINK+1],          /* array containing absolute position matrix of
                                          frame (p) in base frame */
	     mabsinv,                  /* inverse position matrix of the frame
					  positioned in the centre of the
					  gripper */
	     Lrelp,                    /* L relative matrix of p-th joint
					  seen in frame (p-1) */
	     Lrel0,                    /* L relative matrix of p-th joint
                                          seen in base frame */
	     dm,dS,
	     Wrelp_1[MAXLINK+1],       /* array containing relative velocity matrix of
					  frame (p) seen in frame (p-1) */
	     Wrel0[MAXLINK+1],         /* array containing relative velocity matrix of
                                          frame (p) seen in base frame */
	     Wabs[MAXLINK+1],          /* array containing absolute velocity matrix of
					  frame (i) in base frame */
	     Wtar,                     /* target velocity matrix */
	     Hrelp_1[MAXLINK+1],       /* array containing relative acceleration matrix of
					  frame (p) seen in frame (p-1) */
	     Hrel0[MAXLINK+1],         /* array containing relative acceleration matrix of
                                          frame (p) seen in base frame */
	     Habs[MAXLINK+1],          /* array containing absolute acceleration matrix of
					  frame (i) in base frame */
             Htar,                     /* target acceleration matrix */
             dH;                       /* Htar - H~    H~ is the acceleration
						       evaluated with qpp=0 */
        MAT4 Last ={ {0.,1.,0.,0.},    /* transformation matrix from */
                     {0.,0.,1.,0.},    /* frame (6) to gripper       */
                     {1.,0.,0.,0.},    /* element Z-U is in a[6]     */
                     {0.,0.,0.,1.} };  /*                            */
	MAT4 gripper;       /* absolute frame of gripper */
	POINT first=ORIGIN; /* origin of frame 0 with respect to base,
			       Z value is in a[1] */
	MAT4 Aux,                  /* Auxiliary frame, origin in gripper, parallel to base */
	     Waux, Haux;           /* absolute gripper velocity and acceleration in Auxiliary frame */

        FILE *data;                /* file containing robot description */
        FILE *motion;              /* file containing motion description */
        FILE *out;                 /* output file */
        FILE *guess;               /* 1st guess for q */

        if (argc!=5)               /* testing parameters */
	{
		fprintf(stderr,"Usage: Elb_i_dh data_file motion_file q_guess.1st output_file\n");
		exit(1);
	}
	data=fopen(argv[1],"r");
	if(data==NULL)
	{
		fprintf(stderr,"Error: Unable to open data_file\n");
		exit(2);
	}
	motion=fopen(argv[2],"r");
	if(motion==NULL)
	{
		fprintf(stderr,"Error: Unable to open motion_file\n");
		exit(3);
	}
	guess=fopen(argv[3],"r");
	if(guess==NULL)
	{
		fprintf(stderr,"Error: Unable to open 1st_guess_file\n");
		exit(4);
	}
	out=fopen(argv[4],"w");
	if(out==NULL)
	{
		fprintf(stderr,"Error: Unable to open output_file\n");
		exit(5);
	}
	for(p=1;p<=MAXLINK;p++)
		fscanf(data,"%f",&a[p]); /* read robot description */
	for(p=0;p<MAXLINK;p++)
		fscanf(guess,"%f",&q[p]); /* 1st guess for q */

				       /* matrices initialization */
	clear(M Jac,MAXLINK,MAXLINK);

	first[Z]=a[1];
	rotat24(Z,PIG_2,first,mabs[0]); /* position matrix of frame 0 from base frame */
	Last[Z][U]=a[6];                /* gripper position in frame 6 */
        idmat4(Aux);
	clear4(Waux); clear4(Haux);
        clear4(Wabs[0]); /* abs velocity of base and of frame 0 */
        clear4(Habs[0]); /*     acceleration                    */

	a[1]=a[6]=0;    /* D&H parameter 'a' of link 1 and link 6 are zero */

        fscanf(motion,"%f",&dt);                /* read time step */
	fscanf(motion,"%d %d %d",&ii,&jj,&kk);  /* read Cardan convention */
	fprintf(out,"%f\n\n",dt);

	for(t=0;;t+=dt)                /* main loop */
	{
		ierr=fscanf(motion,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
			    &q1[0],&q1[1],&q1[2],&qp1[0],&qp1[1],&qp1[2],
			    &qpp1[0],&qpp1[1],&qpp1[2],&O[X],&O[Y],&O[Z],
			    &vel[0],&vel[1],&vel[2],&acc[0],&acc[1],&acc[2]);
		if (ierr!=18)
			exit(0);
		O[U]=1.;
				       /* builds target position matrix */
		cardantoM(q1,ii,jj,kk,O,mtar);
		vmcopy(M O,3,4,Col,M mtar,4,4);

		for (k=0;k<maxiter;k++)
		{
			for (p=1;p<=MAXLINK;p++)
			{
				       /* builds relative position matrix */
				dhtom(Rev,theta[p],d[p],b[p],a[p],fi[p],
				      q[p-1],mrelp_1[p]);
				       /* builds absolute position matrix */
				molt4(mabs[p-1],mrelp_1[p],mabs[p]);
                                       /* builds relative L matrix in base frame */
				makeL2(Rev,Z,0.,orig,Lrelp);
				       /* builds relative L matrix in frame (p) */
				trasf_mami(Lrelp,mabs[p-1],Lrel0);
				buf[0]=Lrel0[X][U];
				buf[1]=Lrel0[Y][U];
				buf[2]=Lrel0[Z][U];
				buf[3]=Lrel0[Z][Y];
				buf[4]=Lrel0[X][Z];
				buf[5]=Lrel0[Y][X];
				vmcopy(M buf,6,p,Col,M Jac,MAXLINK,
				       MAXLINK);
			}
			molt4(mabs[MAXLINK],Last,gripper);
			sub4(mtar,gripper,dm);
			n=norm4(dm);
				       /* tests if solution has been reached */
			if (n>toll)
			{
				inverse(gripper,mabsinv);
				molt4(dm,mabsinv,dS);
				ds[0]=dS[X][U];
				ds[1]=dS[Y][U];
				ds[2]=dS[Z][U];
				ds[3]=dS[Z][Y];
				ds[4]=dS[X][Z];
				ds[5]=dS[Y][X];
				rank=solve(M Jac,M ds,M dq,MAXLINK);
				       /* builds the joint var. at next step */
				if(rank!=MAXLINK) printf("*** rank is %d: singular position!\a",rank);
				sum(M q,M dq,M q,MAXLINK,1);
			}
			else
				break;
		}
		if (k<maxiter)
		{
		 Aux[X][U]=gripper[X][U];
		 Aux[Y][U]=gripper[Y][U];
		 Aux[Z][U]=gripper[Z][U];

				       /* builds target velocity matrix */
		 cardantoW(q1,qp1,ii,jj,kk,O,Waux);
		 vmcopy(M vel,3,4,Col,M Waux,4,4);
		 trasf_mami(Waux,Aux,Wtar); /* transform velocity from auxiliary frame to base frame*/
				       /* builds target acceleration matrix */
		 cardantoH(q1,qp1,qpp1,ii,jj,kk,O,Haux);
		 vmcopy(M acc,3,4,Col,M Haux,4,4);
		 trasf_mami(Haux,Aux,Htar); /* transform acceleration from auxiliary frame to base frame */

				       /* builds joint velocity array */
		 buf[0]=Wtar[X][U];
		 buf[1]=Wtar[Y][U];
		 buf[2]=Wtar[Z][U];
		 buf[3]=Wtar[Z][Y];
		 buf[4]=Wtar[X][Z];
		 buf[5]=Wtar[Y][X];
		 rank=solve(M Jac,M buf,M qp,MAXLINK);
		 if(rank!=MAXLINK) printf("*** rank is %d: singular position!\a",rank);

		 for(i=1;i<=MAXLINK;i++) /* acceleration */
		 { velacctoWH(Rev,qp[i-1],0.,Wrelp_1[i],Hrelp_1[i]);
				   /* W and H matrices in frame 0 */
		   trasf_mami(Wrelp_1[i],mabs[i-1],Wrel0[i]);
		   trasf_mami(Hrelp_1[i],mabs[i-1],Hrel0[i]);
				   /* absolute velocity and acceleration matrices */
		   sum4(Wabs[i-1],Wrel0[i],Wabs[i]);
		   coriolis(Habs[i-1],Hrel0[i],Wabs[i-1],Wrel0[i],Habs[i]);
		 }
		 sub4(Htar,Habs[MAXLINK],dH);
				       /* builds joint acceleration array */
		 buf[0]=dH[X][U];
		 buf[1]=dH[Y][U];
		 buf[2]=dH[Z][U];
		 buf[3]=dH[Z][Y];
		 buf[4]=dH[X][Z];
		 buf[5]=dH[Y][X];
		 rank=solve(M Jac,M buf,M qpp,MAXLINK);
		 if(rank!=MAXLINK) printf("*** rank is %d: singular position!\a",rank);
		 }
		 else
		 {
			fprintf(stderr,"Newton-Raphson method doesn't converge\n");
			exit(1);
		}

				       /* output results */
		printf("\nTime=%f",t);
                printv("The joint angles q are",q,6);
                printv("The joint velocity qp are",qp,6);
                printv("The joint acceleration qpp are",qpp,6);
		char a; scanf(" %c",&a);

		for (p=0;p<MAXLINK;p++)
		{
			fprintf(out,"%15.5f",q[p]);
			fprintf(out,"%15.5f",qp[p]);
			fprintf(out,"%15.5f\n",qpp[p]);
		}
		fprintf(out,"\n");
	}
	fcloseall();
}