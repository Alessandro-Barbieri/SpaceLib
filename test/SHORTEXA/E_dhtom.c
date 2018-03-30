/* E_dhtom.c
  sample program to test the functions 'dhtom' and 'DHtoMstd'
  they must perform the same result

  direct kinematic of the Stanford Arm
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "SPACELIB.H"

#define MAXLINK 6

int jtype[]={Rev, Rev, Pri, Rev, Rev, Rev};   // joint type
							    /* Denavit & Hartenberg's parameters */
real alpha[]={-PIG_2, PIG_2, 0., -PIG_2, PIG_2, 0.};
real a[]={0., 0., 0., 0., 0., 0};
real d[]={0., 0.2, 0., 0., 0., 0.};
real theta[]={0., 0., 0., 0., 0., 0.};

real q[]={0.11,0.22, 0.33, 0.44, 0.55, 0.66}; /* array of joint position variables */

#define RandData 0     /* 1 = random data, 0 = fixed data */

int main(int argc,char *argv[])
{
	int i;

	MAT4 Ma, Mb, MM, dM;
	MAT4 tmp;

	if (RandData==1) {  // choose fixed or random data
	   srand( (unsigned)time( NULL ) );
	   for (i=0;i<MAXLINK;i++) {
	   	    q[i]=rand()/(1.*RAND_MAX);
	   }
	}

	idmat4(Ma);					/* direct kinematics with 'dhtom' */
	for (i=0;i<MAXLINK;i++) {
		dhtom(jtype[i],theta[i],d[i],0.,a[i],alpha[i],q[i],MM);
		molt4(Ma,MM,tmp);
		mcopy4(tmp,Ma);
	}
	printm4("Ma: result with 'dhtom'",Ma);

	idmat4(Mb);					/* direct kinematics with 'DHtoMstd' */
	for (i=0;i<MAXLINK;i++) {
		if (jtype[i]==Rev){
			DHtoMstd(q[i],d[i],a[i],alpha[i],MM);
		} else {
			DHtoMstd(theta[i],q[i],a[i],alpha[i],MM);
		}
		molt4(Mb,MM,tmp);
		mcopy4(tmp,Mb);
	}
	printm4("Mb: result with 'DHtoMstd",Mb);

	sub4(Ma,Mb,dM);
	printm4("the results must be identical and so dM=Ma-Mb=[0]",dM);
}
