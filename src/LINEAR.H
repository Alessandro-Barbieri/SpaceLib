#ifndef _SPACE_LIN_
#define _SPACE_LIN_
                                          /* --- Speed and Acceleration Matrices --- */
int    dyn_eq(MAT4 J, MAT4 Wp, MAT4 F, int var[2][6]);
                                          /* --- Functions solve minvers linear, solve --- */
int    solve(MAT A, real *b, real *x, int dim); /* solves linear syst. A*x=b */
                                          /* --- Functions solve minvers linear, minvers --- */
int    minvers(MAT A, MAT Ai, int dim);   /* inverse of A using linear */
                                          /* --- Functions solve minvers linear, linear --- */
void   linear(MAT H, int idim,int jdim, int imax, int jmax, int nsol,
              int ivet[],int *irank, real *arm, int vpr[]);
void   swabr(MAT H, int idim, int jdim, int j2max,int k, int im);
void   swabc(MAT H, int idim, int jdim, int imax, int k, int jm, int vet[]);
void   normalr(MAT H, int idim, int jdim, int j2max, int k, real rm);
void   elimin(MAT H, int idim, int jdim, int imax, int j2max, int k);
real   rmax(MAT H, int idim, int jdim, int imin, int imax, int jmin, int jmax,
            int *im, int *jm);
#endif
