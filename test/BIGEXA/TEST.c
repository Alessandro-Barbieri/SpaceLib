/*                                 TEST
		     Program for the trajectory prediction
                  of a two-link system floating in the space

                                October 1990
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "SPACELIB.H"
#include "LINEAR.H"

void delta_m(MAT4 W, MAT4 H, real dt, MAT4 dm);

	POINT O=ORIGIN;
	AXIS  Zax=Zaxis;
	MAT4 m1, m12, m2, W1, W12, W120, W2, H1, H12, H120, H2;
	MAT4 Wp, dm, TMP, UNIT=UNIT4;
	MAT4 J1, J2, J10, J20, Jtot;
	MAT4 F1,F2;

	int var[2][6] = {{1, 1, 1, 1, 1, 1}, {0, 0, 0, 0, 0, 0}};

int main(void)
{
	real q, qp, qpp;
	real t,dt;
	float mass,jx,jy,jz,jxy,jyz,jxz,xg,yg,zg;
	int i,j,n,excode;
	FILE *fil;

	fil=fopen("TEST.DAT","r"); /* open description file */
	if (fil==NULL)
	{
		printf("Error on input file TEST.DAT");
		exit(1);
	}
			     /* read description of the links --- step (1) */
	fscanf(fil,"%f %f %f %f %f %f %f %f %f %f",
	       &mass, &jx,&jy,&jz, &jxy,&jyz,&jxz, &xg,&yg,&zg); /* link 1 */
	jtoJ(mass,jx,jy,jz,jxy,jyz,jxz,xg,yg,zg,J1); /*build inertia matrix*/

	fscanf(fil,"%f %f %f %f %f %f %f %f %f %f",
	       &mass, &jx,&jy,&jz, &jxy,&jyz,&jxz, &xg,&yg,&zg); /* link 2 */
	jtoJ(mass,jx,jy,jz,jxy,jyz,jxz,xg,yg,zg,J2); /*build inertia matrix*/

			     /* read initial condition of the system */
	for(i=0;i<4;i++)     /* read velocity matrix of link 1 */
		for(j=0;j<4;j++)
			fscanf(fil,"%lf",&W1[i][j]);
	for(i=0;i<4;i++)     /* read position matrix of link 1 */
		for(j=0;j<4;j++)
			fscanf(fil,"%lf",&m1[i][j]);
	fclose(fil);

	fil=fopen("TEST.MOT","r"); /* open motion file */
	if (fil==NULL)
	{
		printf("error on input file TEST.MOT");
		exit(1);
	}
	fscanf(fil,"%lf",&dt);      /* read integration step "dt" */
	for(t=0;;t+=dt)            /* loop for each instant t --- step (2) */
	{
		n=fscanf(fil,"%lf %lf %lf ",&q,&qp,&qpp); /* read motion of motor (a) */
		printf("--- t: %f ---- q: %f   qp: %f   qpp: %f",t,q,qp,qpp);

		screwtom(Zax, q, 0., O, m12);     /* relative position links 1&2 (b) */
		molt4(m1,m12,m2);                      /* absolute position	of link 2 (c) */
		molt4(W1,W1,H1);                       /* partial acceleration of link 1 (d) */
		velacctoWH2(Rev,Z,qp,qpp,W12,H12);     /* relative velocity & acceleration of link 1&2 (e) */
		trasf_mami(W12,m1,W120);
		trasf_mami(H12,m1,H120);               /* (f) */
		norm_simm_skew(M W120,3,4,SKEW);       /* normalization reducing numerical error */
			     /* absolute velocity and partial acceleration of link 2 (g) */
		sum4(W1,W120,W2);
		coriolis(H1,H120,W1,W120,H2);

			     /* refer inertia moment to absolute frame (h) */
		trasf_mamt4(J1,m1,J10);
		trasf_mamt4(J2,m2,J20);
		norm_simm_skew(M J10,4,4,SYMM);	       /* normalization reducing numerical errors */
		norm_simm_skew(M J20,4,4,SYMM);	       /* normalization reducing numerical errors */

		mcopy4(J10,Jtot);                      /* total inertia (i)*/
		sum4(Jtot,J20,Jtot);

		skew4(H1,J10,F1);                      /* evaluate inertia actions (j) */
		skew4(H2,J20,F2);
		sum4(F1,F2,F1);

		excode=dyn_eq(Jtot,Wp,F1,var);
		rmolt4(Wp,-1,Wp); 		       /* evaluate Wp (k) */

		if(excode==NOTOK)
		{
			printf(" excode=%d",excode);
			exit(1);
		}

		sum4(Wp,H1,H1);                        /* absolute acceleration of links 1 & 2 (l) */
		sum4(Wp,H2,H2);

		printm4("Position matrix of link 1",m1);
		printm4("Absolute position matrix of link 2",m2);
		printm4("Velocity matrix of link 1",W1);
		printm4("Absolute velocity matrix of link 2",W2);
		printm4("Acceleration matrix of link 1",H1);
		printm4("Absolute acceleration matrix of link 2",H2);

		if(n!=3) break;               	       /*if motion file empty -> end of loop */
		delta_m(W1,H1,dt,dm);                  /* builds matrix dm (m) */
		molt4(dm,m1,TMP);                      /* new position of link 1 (n) */
		mcopy4(TMP,m1);

		rmolt4(Wp,dt,Wp);                      /* new velocity of link 1 (o) */
		sum4(Wp,W1,W1);
	}
	exit(0);
}


/* --- Function delta_m: builds matrix dm where dm=[1]+Wdt+0.5Hdt^2 --- */

void delta_m(MAT4 W, MAT4 H, real dt, MAT4 dm)
{
	int i,j;
	real dt2;

	dt2=dt*dt;
	for (i=0; i<3; i++)
		for (j=0; j<4; j++)
			dm[i][j]=UNIT[i][j]+W[i][j]*dt+.5*H[i][j]*dt2;
	dm[U][X]=dm[U][Y]=dm[U][Z]=0.;
	dm[U][U]=1;
	normal4(dm);
}
