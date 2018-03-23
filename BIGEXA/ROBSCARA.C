/*
----------------------------------------------------------------------------
ROBSCAR.C:  Sample program for direct kinematics of Scara robot
	    (See User's Manual)

		University of Brescia
		Mechanical Eng. Department
		Via Branze 38
		25123 BRESCIA - ITALY

		giovanni.legnani @ bsing.unibs.it
-----------------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#include "spacelib.h"

void main(void)
{
	real q[3] ={PIG/4,PIG/6,-0.5};      /* joint variables array */
	real qp[3]={PIG*5/4,PIG*5/4,-0.5};  /* joint var. first time der. */
	POINT O1={0.,0.,1.5,1.};          /* origin of frame (1) in frame (0) */
	POINT O2={0.33,0.,0.,1.};         /* origin of frame (2) in frame (1) */
	POINT O3={0.33,0.,0.,1.};         /* origin of frame (3) in frame (2) */
	POINT O4={0.,0., 0., 1.};         /* origin of frame (4) in frame (3) */
	POINT Oa;                         /* origin of frame (a) in frame (0) */
	POINT O=ORIGIN;
	MAT4 m01,m12,m23,m34;           /* mi-1,i  : pos. of frame (i) in frame (i-1) */
	MAT4 m02,m03,m04,m0a;           /* m0,i    : pos. of frame (i) in frame (0) */
	MAT4 L12r,L23r,L34r;            /* Li,i+1r : L mat. of joint i in frame (i-1) */
	MAT4 L12f,L23f,L34f;            /* Li,i+1f : L mat. of joint i in frame (0) */
	MAT4 W01,W12,W23,W34;           /* Wi,i+1  : W mat. of frame (i+1) in frame (i) */
	MAT4 W04;                       /* velocity matrix of the gripper in frame (0) */
	MAT4 W04a;                      /* velocity matrix of the gripper in frame (a) */


					/* builds relative position matrices*/
	idmat4(m01);
	idmat4(m34);
	vmcopy(M O1,4,4,Col,M m01,4,4);
	rotat34(Z,q[0],O2,m12);
	rotat34(Z,q[1],O3,m23);
	vmcopy(M O4,4,4,Col,M m34,4,4); m34[Z][U]+=q[2];

	molt4(m01,m12,m02);             /* builds absolute position matrices */
	molt4(m02,m23,m03);
	molt4(m03,m34,m04);

	idmat4(m0a);                    /* builds pos. matrix of frame (a) in frame (0) */
	mvcopy(M m04,4,4,4,Col,M Oa);
	vmcopy(M Oa,4,4,Col,M m0a,4,4);

	makeL2(Rev,Z,0.,O,L12r);        /* builds relative L matrices */
	makeL2(Rev,Z,0.,O,L23r);
	makeL2(Pri,Z,0.,O,L34r);

	trasf_mami(L12r,m01,L12f);      /* evaluates L matrices in frame (0) */
	trasf_mami(L23r,m02,L23f);
	trasf_mami(L34r,m03,L34f);

	clear4(W01);                    /* builds relative W matrices */
	rmolt4(L12f,qp[0],W12);
	rmolt4(L23f,qp[1],W23);
	rmolt4(L34f,qp[2],W34);

	sum4(W01,W12,W04);              /* builds abs. W matrix of frame (4) in frame (0) */
	sum4(W04,W23,W04);
	sum4(W04,W34,W04);

	trasf_miam(W04,m0a,W04a);       /* evaluates W matrix of frame (4) in frame (a) */

					/* output results */
	printm4("The velocity matrix of the gripper in frame (0) is:",W04);
	printm4("The velocity matrix of the gripper in frame (a) is:",W04a);

}
