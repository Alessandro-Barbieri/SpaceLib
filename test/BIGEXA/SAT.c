#include <stdio.h>
#include <math.h>
#include "SPACELIB.H"

int main (void)

#define sb      (real) sin(beta)
#define cb      (real) cos(beta)
#define sb_a    (real) sin(beta-alpha)
#define cb_a    (real) cos(beta-alpha)
{
	POINT P1={0.875, 2.100, 1.500, 1. }; /* values of initial configuration */
	POINT P2={1.750, 2.100, 3.200, 1. };
	MAT4 mi={ {1.,  0., 0., 0.875},
			  {0.,  0., 1., 2.100},
			  {0., -1., 0., 1.500},
			  {0.,  0., 0., 1.   } };
	real alpha=atan2( P2[Z]-P1[Z], P2[X]-P1[X] );
	real beta=atan2( P2[Y], P2[X] );
	real d=dist(P1,P2);

	MAT4 m4, Q5, mf, miinv;              /* elements used in computation */
	AXIS u5;

	MAT4 Qtot;                           /* global elements which are evaluated */
	AXIS utot;
	real phi,h;
	POINT P;

   /* from initial configuration to step 4 */
	m4[X][X] = -cb_a; m4[X][Y] = -sb_a; m4[X][Z] =  0.; m4[X][U] = P2[X]+d*cb;
	m4[Y][X] = -sb_a; m4[Y][Y] =  cb_a; m4[Y][Z] =  0.; m4[Y][U] = P2[Y]+d*sb;
	m4[Z][X] =  0.;   m4[Z][Y] =  0.;   m4[Z][Z] = -1.; m4[Z][U] = P2[Z];
	m4[U][X] =  0.;   m4[U][Y] =  0.;   m4[U][Z] =  0.; m4[U][U] = 1.;

                                             /* step 5 */
    u5[X]=sb; u5[Y]=-cb; u5[Z]=0.;
	screwtom(u5,rad(26.),0.,P2,Q5);

                                             /* final configuration */
	molt4(Q5,m4,mf);

                                             /* evaluate global elements */
	invers(mi,miinv);
	molt4(mf,miinv,Qtot);
	mtoscrew(Qtot,utot,&phi,&h,P);

                                             /* output results */
	printf("\n*** Results ***\n");
	printv("The rototranslation axis is :",utot,3);
	printf("\nThe rotation angle about this axis is :"PRIr" [deg]\n",deg(phi));
	printf("\nThe translation about this axis is :"PRIr"\n",h);
    printv("The point P of the axis is :",P,4);

	#undef sb
	#undef cb
	#undef sb_a
	#undef cb_a
	return(0);
}
