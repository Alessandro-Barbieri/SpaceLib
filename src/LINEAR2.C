/*
------------------------------------------------------------------------------
			LINEAR2.C
			© copyright by

			D.Amadori, P.Ghislotti, G.Pugliese

			version v1.0 - April 1997
			Developed for :

			Prof. G.LEGNANI
			Università degli Studi di Brescia
			Dipartimento di Ingegneria Meccanica ed Industriale
			Via Branze 38
			25123 BRESCIA - ITALY

			giovanni.legnani @ unibs.it

------------------------------------------------------------------------------
			Do not remove this copyright notice
------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include "SPACELIB.H"
#include "LINEAR.H"


/* --- Functions solve minvers linear, solve --- */

/* == solve =============================================================== */
/* function for the solution of a "standard" linear system in the form
			  A*x=b
    Input paremeters
      real A[][] : square matrix containing the system coefficients
      real b[]   : array of the right hand coefficients
      int  dim   : system dimension
    Output parameters
      real x[]   : array where the solution is stored

    Note : this function is used in order to solve a "standard" linear system
	   using the more complex function linear (function linear can deal
	   with more general cases). The function returns the rank of matrix A
*/

int solve(MAT A, real *b, real *x, int dim)
{
	int solution=1;
	real *H;
	int *ivet;
	int vpr=-1;
	int irank;
	real arm;
	int i,j;

	H = (real *) calloc(dim*(dim+solution),sizeof(real));
	ivet = (int *)  calloc(dim,sizeof(int));

	for(i=0;i<dim;i++)          /* copies matrix A into matrix H */
		for(j=0;j<dim;j++)
			*(H+i*(dim+solution)+j) = *(A+i*dim+j);

	for(i=0;i<dim;i++)	    /* copies right hand vector b into H */
		H[i*(dim+solution)+dim] = b[i];

	linear(H,dim,dim+solution,dim,dim,solution,ivet,&irank,&arm,&vpr);

	for(i=0;i<dim;i++)          /* finds solution x reordering matrix H */
		x[ivet[i]] = H[i*(dim+solution)+dim];

	free(H); free(ivet);

	return irank;
}


/* --- Functions solve minvers linear, minvers --- */

/* == minvers ========================================================= */
/* function that builds the inverse of a matrix A
    Input parameters
      real A[][]  : initial square matrix whose inverse must be found
      int dim     : dimension of matrix A
    Output parameters
      real Ai[][] : matrix where the inverse is stored

    Note : This function calls the subroutine linear contained in file
	   linear.c. It uses a general property of the system A*x=I where A is
	   the starting matrix and I is the identity matrix. The only solution
	   of this system, whenever the rank of A is dim, is the inverse of
	   the matrix A. This function also returns the rank of matrix A
*/

int minvers(MAT A, MAT Ai, int dim)
{

	real *H;                /* matrix where A and I must be copied */
	int *ivet;              /* vector filled with reordering info. */
	int vpr = -1;
	int irank;
	real arm;
	int i,j;

	H = (real *) calloc (2*dim*dim,sizeof(real));
	ivet = (int *) calloc (dim,sizeof(int));

	for (i=0;i<dim;i++)     /* copies matrix A into matrix H */
		for (j=0;j<dim;j++)
			H[i*dim*2+j] = *(A+i*dim+j);
	for (i=0;i<dim;i++)     /* copies identity matrix I into matrix H */
		for (j=dim;j<2*dim;j++)
			if (i!=(j-dim))
				H[i*dim*2+j] = 0;
			else
				H[i*dim*2+j] = 1;
	linear(H,dim,2*dim,dim,dim,dim,ivet,&irank,&arm,&vpr);
	for (j=dim;j<2*dim;j++) /* extracts solution x reordering matrix H */
		for (i=0;i<dim;i++)
			*(Ai+ivet[i]*dim+j-dim) = *(H+i*2*dim+j);

	free (H); free (ivet);

	return irank;
}
