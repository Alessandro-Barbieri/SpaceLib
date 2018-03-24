/* elb_d_PA.c program for the DIRECT kinematics of ELBOW robot.
   the output of this program is compatible with the input  of elb_i_dh.c
   the input  of this program is compatible with tho output of elb_i_dh.c */

/* Note: To compile this program the type real must be set equivalent to the type float
	 (see also User's Manual). This is necessary because the formatting string of the
	 fscanf function have been written using %f as descriptor. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPACELIB.H"

#define MAXLINK 6

int main(int argc,char *argv[])
{
	int axis[MAXLINK+2]={U,Z,X,X,X,Z,X,U};
	int i,j,k,ierr;
        int ii=X;                  /* Euler/Cardan convention */
        int jj=Y;                  /* for gripper             */
        int kk=Z;                  /* angular position        */

	real a[MAXLINK+1];         /* array of link lengths */
	POINT O[MAXLINK+2];        /* array containing orig. of ref. frames */
	real q,qp,qpp;             /* joint variables : position, velocity, acceleration */
	real t,dt;

	VECTOR q1,q2,              /* Euler/Cardan angle configurations */
	       qp1,qp2,            /* angles first time derivative */
	       qpp1,qpp2;          /* angles second time derivative*/

	MAT4 mreli_1[MAXLINK+2],   /* array containing position matrix of frame
				      (i) seen in frame (i-1) */
	     mabs[MAXLINK+2],      /* array containing absolute position matrix of
				      frame (i) in frame (0) */
	     Wreli_1[MAXLINK+2],   /* array containing relative velocity matrix of
				      frame (i) seen in frame (i-1) */
	     Wrel0[MAXLINK+2],     /* array containing relative velocity matrix of
				      frame (i) seen in frame (0) */
	     Wabs[MAXLINK+2],      /* array containing absolute velocity matrix of
				      frame (i) in frame (0) */
	     Hreli_1[MAXLINK+2],   /* array containing relative acceleration matrix of
				      frame (i) seen in frame (i-1) */
	     Hrel0[MAXLINK+2],     /* array containing relative acceleration matrix of
				      frame (i) seen in frame (0) */
	     Habs[MAXLINK+2];      /* array containing absolute velocity matrix of
				      frame (i) in frame (0) */
        MAT4 Aux,                  /* Auxiliary frame, origin in gripper, parallel to base */
	     Waux, Haux;           /* absolute gripper velocity and acceleration in Aux frame */

	FILE *data;                /* file containing robot description */
	FILE *motion;              /* file containing motion description */
	FILE *out;                 /* output file */

	if (argc!=4)               /* testing parameters */
	{
		fprintf(stderr,"Usage: Elb_d_pa data_file joint_motion_file griper_output_file\n");
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

				   /* matrices initialization */
	idmat4(mabs[0]); idmat4(Aux);
	clear4(Wabs[0]);
	clear4(Habs[0]);
	clearv(O[0]);
	O[0][U]=1.;

        for(i=1;i<=MAXLINK;i++)         /* read link lengths */
	{
		fscanf(data,"%f",&a[i]);
	}

				   /* relative origin of frame (i) in (i-1) */
	O[1][X]=0.;   O[1][Y]=0.;   O[1][Z]=0.;   O[1][U]=1.;
	O[2][X]=0.;   O[2][Y]=0.;   O[2][Z]=a[1]; O[2][U]=1.;
	O[3][X]=0.;   O[3][Y]=a[2]; O[3][Z]=0.;   O[3][U]=1.;
	O[4][X]=0.;   O[4][Y]=a[3]; O[4][Z]=0.;   O[4][U]=1.;
	O[5][X]=0.;   O[5][Y]=a[4]; O[5][Z]=0.;   O[5][U]=1.;
	O[6][X]=0.;   O[6][Y]=0.;   O[6][Z]=0.;   O[6][U]=1.;
	O[7][X]=a[6]; O[7][Y]=0.;   O[7][Z]=0.;   O[7][U]=1.;

        fscanf(motion,"%f",&dt);       /* read time step */
        fprintf(out,"%f\n",dt);        /* write dt to out file */
        fprintf(out,"%d %d %d\n\n",ii,jj,kk);  /* write gripper Cardan convention to out file */
	for(t=0;eof(motion);t+=dt) /* main loop */
	{
		for(i=1;i<=MAXLINK;i++)
		{
			ierr=fscanf(motion,"%f %f %f",&q,&qp,&qpp);
			if(ierr!=3)
				exit(0);

				   /* builds relative position matrix */
			rotat24(axis[i],q,O[i],mreli_1[i]);

				   /* builds relative velocity and acceleration matrices */
			velacctoWH3(Rev,axis[i],qp,qpp,O[i],Wreli_1[i],
				    Hreli_1[i]);

				   /* absolute position matrix of frame (i) */
			molt4(mabs[i-1],mreli_1[i],mabs[i]);
				   /* W and H matrices in frame 0 */
			trasf_mami(Wreli_1[i],mabs[i-1],Wrel0[i]);
			n_skew34(Wrel0[i]);
			trasf_mami(Hreli_1[i],mabs[i-1],Hrel0[i]);

				   /* absolute velocity and acceleration matrices */
			sum4(Wabs[i-1],Wrel0[i],Wabs[i]);
			coriolis(Habs[i-1],Hrel0[i],Wabs[i-1],Wrel0[i],
				 Habs[i]);
		}
				   /* relative position matrix of frame (7) */
		idmat4(mreli_1[MAXLINK+1]);
		mreli_1[MAXLINK+1][X][U]=O[7][X];
				   /* relative W and H matrices of frame (7) */
		clear4(Wreli_1[MAXLINK+1]);
		clear4(Hreli_1[MAXLINK+1]);
		molt4(mabs[MAXLINK],mreli_1[MAXLINK+1],mabs[MAXLINK+1]);

				   /* W and H matrices in frame (0) */
		trasf_mami(Wreli_1[MAXLINK+1],mabs[MAXLINK],Wrel0[MAXLINK+1]);
		n_skew34(Wrel0[MAXLINK]);
		trasf_mami(Hreli_1[MAXLINK+1],mabs[MAXLINK],Hrel0[MAXLINK+1]);

				   /* absolute velocity and acceleration matrices */
		sum4(Wabs[MAXLINK],Wrel0[MAXLINK+1],Wabs[MAXLINK+1]);
		coriolis(Habs[MAXLINK],Hrel0[MAXLINK+1],Wabs[MAXLINK],
			 Wrel0[MAXLINK+1],Habs[MAXLINK+1]);

				   /* extracts Cardan angles */
		Htocardan(mabs[MAXLINK+1],Wabs[MAXLINK+1],Habs[MAXLINK+1],
			  ii,jj,kk,q1,q2,qp1,qp2,qpp1,qpp2);

		Aux[X][U]=mabs[MAXLINK+1][X][U];
		Aux[Y][U]=mabs[MAXLINK+1][Y][U];
		Aux[Z][U]=mabs[MAXLINK+1][Z][U];
		trasf_miam(Wabs[MAXLINK],Aux,Waux); /* transform velocity */
		trasf_miam(Habs[MAXLINK],Aux,Haux); /* and acceleration in auxiliary frame */
					  /* extracts Cardan angles (and their
					     time derivatives) of gripper  */
		Htocardan(mabs[MAXLINK+1],Waux,Haux,
			  ii,jj,kk,q1,q2,qp1,qp2,qpp1,qpp2);

				   /* output results */
		printf("Time=%f\n",t);
		printm4("The position matrix of the gripper is:",mabs[MAXLINK+1]);
		printm4("The velocity matrix of the gripper is:",Waux);
		printm4("The acceleration matrix of the gripper is:",Haux);
		printf("\nPress any key to continue\n");
		char a; scanf(" %c",&a);

		fprintf(out,"%f %f %f\n",q1[0],q1[1],q1[2]);
		fprintf(out,"%f %f %f\n",qp1[0],qp1[1],qp1[2]);
		fprintf(out,"%f %f %f\n",qpp1[0],qpp1[1],qpp1[2]);

		fprintf(out,"%f %f %f\n",mabs[MAXLINK+1][X][U],mabs[MAXLINK+1][Y][U],
						       mabs[MAXLINK+1][Z][U]);
		fprintf(out,"%f %f %f\n",Waux[X][U],Waux[Y][U],Waux[Z][U]);
		fprintf(out,"%f %f %f\n",Haux[X][U],Haux[Y][U],Haux[Z][U]);
	}
	fcloseall();
}
