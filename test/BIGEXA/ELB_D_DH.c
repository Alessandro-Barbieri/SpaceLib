/* elb_d_DH.c program for the DIRECT kinematics of ELBOW robot.
   Frames assigned according to Denavit and Hartenberg conventions.
   the output of this program is compatible with the input  of elb_i_dh.c
   the input  of this program is compatible with the output of elb_i_dh.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPACELIB.H"

#define MAXLINK 6

int main(int argc,char *argv[])
{
	int i,j,k,ierr;
	int ii=X;                  /* Euler/Cardan convention */
	int jj=Y;                  /* for gripper             */
	int kk=Z;                  /* angular position        */

				   /* Denavit & Hartenberg's parameters (D&H) */
	real theta[MAXLINK+1]={0.,0.,0.,0.,0.,0.,0.};
	real d[MAXLINK+1]={0.,0.,0.,0.,0.,0.,0.};
	real b[MAXLINK+1]={0.,0.,0.,0.,0.,0.,0.};
	real fi[MAXLINK+1]={0.,PIG_2,0.,0.,3*PIG_2,PIG_2,0.};
	real a[MAXLINK+1];         /* to be read from data_file */

	real q,qp,qpp;             /* joint variables : position, velocity, acceleration */
	real t,dt;

	VECTOR q1,q2;              /* Euler/Cardan angle configurations */
	VECTOR qp1,qp2;            /* Euler/Cardan angle first time derivative */
	VECTOR qpp1,qpp2;          /* Euler/Cardan angle second time derivative */

	MAT4 mreli_1[MAXLINK+1],   /* array containing position matrix of frame (i) seen in frame (i-1) */
	     mabs[MAXLINK+1],      /* array containing absolute position matrix of frame (i) in base frame */
	     Wreli_1[MAXLINK+1],   /* array containing relative velocity matrix of frame (i) seen in frame (i-1) */
	     Wrel0[MAXLINK+1],     /* array containing relative velocity matrix of frame (i) seen in base frame */
	     Wabs[MAXLINK+1],      /* array containing absolute velocity matrix of frame (i) in base frame */
	     Hreli_1[MAXLINK+1],   /* array containing relative acceleration matrix of frame (i) seen in frame (i-1) */
	     Hrel0[MAXLINK+1],     /* array containing relative acceleration matrix of frame (i) seen in base frame */
	     Habs[MAXLINK+1];      /* array containing absolute velocity matrix of frame (i) in base frame */
	MAT4 Last ={ {0.,1.,0.,0.},      /* transformation matrix from */
		     {0.,0.,1.,0.},      /* frame (6) to gripper       */
		     {1.,0.,0.,0.},      /* element Z-U is in a[6]     */
		     {0.,0.,0.,1.} };    /*                            */
	MAT4 gripper;              /* absolute frame of gripper */
	POINT first=ORIGIN;        /* origin of frame 0 with respect to base, Z value is in a[1] */
	MAT4 Aux,                  /* Auxiliary frame, origin in gripper, parallel to base */
	     Waux, Haux;           /* absolute gripper velocity and acceleration in Auxiliary frame */
	FILE *data;                /* file containing robot description */
	FILE *motion;              /* file containing motion description */
	FILE *out;                 /* output file */

	if (argc!=4)               /* testing parameters */
	{
		fprintf(stderr,"Usage: elb_d_dh data_file joint_motion_file gripper_output_file\n");
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
	out=fopen(argv[3],"w");
	if(out==NULL)
	{
		fprintf(stderr,"Error: Unable to open output_file\n");
		exit(4);
	}

	a[0]=0;
        for(i=1;i<=MAXLINK;i++)          /* read link lengths */
	{
		fscanf(data,"%f",&a[i]);
	}
				   /* matrices initialization */
	first[Z]=a[1];
	rotat24(Z,PIG_2,first,mabs[0]); /* position matrix of frame 0 from base frame */
	Last[Z][U]=a[6];                /* gripper position in frame 6 */
	idmat4(Aux);
	clear4(Wabs[0]);               /* absolute velocity of base and of frame 0 */
	clear4(Habs[0]);               /*      acceleration */

	a[1]=a[6]=0;    /* D&H parameter 'a' of link 1 and link 6 are zero */

	fscanf(motion,"%f",&dt);       /* read time step */
	fprintf(out,"%f\n",dt);        /* write dt to out file */
	fprintf(out,"%d %d %d\n\n",ii,jj,kk);  /* write Cardan convention to out file */

	for(t=0;!feof(motion);t+=dt) /* main loop */
	{
		for(i=1;i<=MAXLINK;i++)
		{                  /* read joint motion */
			ierr=fscanf(motion,"%lf %lf %lf\n",&q,&qp,&qpp);
			if(ierr!=3)
				exit(0);

				   /* builds relative position matrix */
			dhtom(Rev,theta[i],d[i],b[i],a[i],fi[i],q,mreli_1[i]);

				   /* builds relative velocity and acceleration matrices */
			velacctoWH(Rev,qp,qpp,Wreli_1[i],Hreli_1[i]);

				   /* absolute position matrix of frame (i) */
			molt4(mabs[i-1],mreli_1[i],mabs[i]);

				   /* W and H matrices in base frame */
			trasf_mami(Wreli_1[i],mabs[i-1],Wrel0[i]);
			trasf_mami(Hreli_1[i],mabs[i-1],Hrel0[i]);

				   /* absolute velocity and acceleration matrices */
			sum4(Wabs[i-1],Wrel0[i],Wabs[i]);
			coriolis(Habs[i-1],Hrel0[i],Wabs[i-1],Wrel0[i],Habs[i]);
		}
		molt4(mabs[MAXLINK],Last,gripper); /* gripper position */
		Aux[X][U]=gripper[X][U];
		Aux[Y][U]=gripper[Y][U];
		Aux[Z][U]=gripper[Z][U];
		trasf_miam(Wabs[MAXLINK],Aux,Waux); /* transform velocity */
		trasf_miam(Habs[MAXLINK],Aux,Haux); /* and acceleration in auxiliary frame */
					  /* extracts Cardan angles (and their time derivatives) of gripper  */
		Htocardan(gripper,Waux,Haux,
			  ii,jj,kk,q1,q2,qp1,qp2,qpp1,qpp2);

				   /* output results */
		printf("Time=%f\n",t);
		printm4("The position matrix of the gripper is:",gripper);
		printm4("The velocity matrix of the gripper is:",Waux);
		printm4("The acceleration matrix of the gripper is:",Haux);

		fprintf(out,"%f %f %f\n",q1[0],q1[1],q1[2]);
		fprintf(out,"%f %f %f\n",qp1[0],qp1[1],qp1[2]);
		fprintf(out,"%f %f %f\n",qpp1[0],qpp1[1],qpp1[2]);

		fprintf(out,"%f %f %f\n",gripper[X][U],gripper[Y][U],
						       gripper[Z][U]);
		fprintf(out,"%f %f %f\n",Waux[X][U],Waux[Y][U],Waux[Z][U]);
		fprintf(out,"%f %f %f\n",Haux[X][U],Haux[Y][U],Haux[Z][U]);
	}
	exit(0);
}
