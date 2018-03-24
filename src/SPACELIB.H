        /* ============= MATH MACRO, CONSTANTS AND TYPE DEFINITION ================ */

#ifndef _SPACELIB_                                             /* revision february 2006*/
#define _SPACELIB_
#ifdef __TURBOC__
#define _BORLAND_
#endif
#ifdef FLOAT
typedef float real;
#else 
typedef double real;
#endif
#define aabs(a)   ((double)(a)>0.  ? (a) : -(a))
                   /* conditional compilation to avoid conflict with standard definitions of min
                      and maximum SPACELIB © definitions coincide with ANSI definitions */
#ifndef max
#define max(a,b) ((a)>(b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a)<(b) ? (a) : (b))
#endif
#define sign(x) ((x)<0 ? -1 : ((x)>0 ? 1 :0))
#define  PIG_2 1.5707963267948966192313216916397514420
#define  PIG   3.1415926535897932384626433832795028841
#define  PIG2  6.2831853071795864769252867665590057683
#define rad(x)  ((x)/180.*PIG)           /* degrees to radians */
#define deg(x)  ((x)*180./PIG)           /* radians to degrees */
#define M (real *)                       /* to recast pointers when calling functions
                                           (used for adjustable dimension matrices) */
typedef real *MAT;                       /* generic matrix                 */
typedef real MAT3[3][3];                 /* rotation submatrices           */
typedef real MAT4[4][4];                 /* transform matrices (& others)  */
typedef real POINT[4];                   /* homogeneous points coordinates */
typedef real AXIS[3];                    /* unit vectors (axes)            */
typedef real VECTOR[3];                  /* vectors                        */
typedef real (*MAT4P)[4];                /* used to point to 4×4 matrices  */
typedef real (*MAT3P)[3];                /* used to point to 3×3 matrices  */
typedef struct {            /* Definition of the LINE type following the parametric definition*/
                  POINT P;               /*  X = Xp+t*alpha ** Y = Yp+t*beta ** Z = Zp+t*gamma*/
                  VECTOR dir;            /* where P=(Xp,Yp,Zp) and dir=[alpha,beta,gamma]*/
                } LINE, * LINEP;         /* LINEP is a pointer to LINE structure type */
typedef real PLANE[4];   /* PLANE type is defined as a four real element array. The first three
                            define a unit vector orthogonal to the plane itself. The fourth real
                            is the distance of plane  from the origin of reference frame.
                            The distance is positive when the vector directed from plane to
                            origin has the same direction of the plane unit vector */
#define OK      1                        /*    successful completation of a function */
#define NOTOK   0                        /* un-successful completation of a function */
#define SYMM    1                        /* denote symmetrical or skew-symmetrical matrices */
#define SKEW   -1
#define Row     0                        /* denote a matrix row */
#define Col     1                        /* denote a matrix column */
#define Rev     0                        /* Revolute (and screw) pairs */
#define Pri     1                        /* Prismatic (Sliding) pairs */
#define Tor     0                        /* Torque */
#define For     1                        /* Force  */
#define X  0                             /* homogeneous coordinates */
#define Y  1
#define Z  2
#define U  3
#define Xaxis  {1.,0.,0.}                /* positive axes unit vectors */
#define Yaxis  {0.,1.,0.}
#define Zaxis  {0.,0.,1.}
#define Xaxis_n  {-1.,0.,0.}             /* negative axes unit vectors */
#define Yaxis_n  {0.,-1.,0.}
#define Zaxis_n  {0.,0.,-1.}
#define ORIGIN {0,0,0,1}
                                         /* identity and null matrices 3×3 and 4×4 */
#define NULL3 {{0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.}}
#define UNIT3 {{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}
#define NULL4 {{0.,0.,0.,0.}, {0.,0.,0.,0.}, {0.,0.,0.,0.}, {0.,0.,0.,0.}}
#define UNIT4 {{1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}}

/* ============== MACRO FOR ADJUSTABLE-DIMENSIONS MATRICES ================ */
            /* Note : 3/4 suffixes refer to 3×3 or 4×4 matrices */
                                                    /* --- norm of a matrix --- */
#define norm3(a)           norm(M a,3,3)
#define norm4(a)           norm(M a,4,4)
                                                    /* --- Rotation and Position Matrices --- */
#define rotat23(a,q,R)     rotat2(a,q,M R,3);       /* rotation matrix */
                                             /* --- Matrix Transformations, Normalization --- */
#define normal3(a)         normal(M a,3)            /* normalizes rotation matrices */
#define normal4(a)         normal(M a,4)  /* normalizes rotation part. of transformation matrices */
                                                  /* normalizes symmetrical and skew symmetrical matrices */
#define n_simm3(m)         norm_simm_skew(M m,3,3, SYMM);
#define n_simm34(m)        norm_simm_skew(M m,3,4, SYMM);
#define n_simm4(m)         norm_simm_skew(M m,4,4, SYMM);
#define n_skew3(m)         norm_simm_skew(M m,3,3, SKEW);
#define n_skew34(m)        norm_simm_skew(M m,3,4, SKEW);
#define n_skew4(m)         norm_simm_skew(M m,4,4, SKEW);
                                                                    /* standard D&H matrix */
#define DHtoMstd(theta,d,a,alpha,MM)      dhtom(0,theta,d,0.,a,alpha,0.,MM)   

                                      /* --- Matrix Transformations, Change of Reference --- */
#define trasf_mamt4(a,b,c) trasf_mamt(M a,M b,M c,4) /* change of reference for 4×4
                                                        contra-variant matrices */
                                        /* --- Matrix transformations, General Operations --- */
#define mtov3(v,m)         mtov(M m,3,v)             /* skew matrix 3×3 to vector */
#define mtov4(v,m)         mtov(M m,4,v)             /* skew matrix 4×4 to vector */
#define vtom3(v,m)         vtom(v,M m,3)             /* vector[3] to 3×3 skew matrix */
#define vtom4(v,m)         vtom(v,M m,4)             /* vector[3] to 4×4 skew matrix */
#define skew4(a,b,c)       skew(M a,M b,M c,4)       /* performs c=skew(a,b) for 4×4 matrices */
                               /* --- Conversion from Cardan Angles to Matrices, Position --- */
#define cardantor3(a,i,j,k,m)   cardantor(a,i,j,k,M m,3)
#define cardantor4(a,i,j,k,m)   cardantor(a,i,j,k,M m,4)
#define rtocardan3(m,i,j,k,a,b) rtocardan(M m,3,i,j,k,a,b)
#define rtocardan4(m,i,j,k,a,b) rtocardan(M m,4,i,j,k,a,b)
#define cardtoM(q,O,m)          cardantoM(q,X,Y,Z,O,m);    /* position matrix from T-B angles */
#define eultoM(q,O,m)           cardantoM(q,Z,X,Z,O,m);    /* position matrix from Euler angles */
#define nauttoM(q,O,m)          cardantoM(q,Z,Y,X,O,m);    /* position matrix from nautical angles */
#define Mtocard(m,q1,q2)        Mtocardan(m,X,Y,Z,q1,q2);  /* T-B angles from position matrix */
#define Mtoeul(m,q1,q2)         Mtocardan(m,Z,X,Z,q1,q2);  /* Euler angles from position matrix */
#define Mtonaut(m,q1,q2)        Mtocardan(m,Z,Y,X,q1,q2);  /* nautical angles from position matrix */
                          /* --- Conversion from Cardan Angles to Matrices, Velocity and Acceleration --- */
                                                          /* W from T-B angles and derivative */
#define cardtoW(q,qp,O,W)          cardantoW(q,qp,X,Y,Z,O,W);
                                                        /* W from Euler angles and derivative */
#define eultoW(q,qp,O,W)           cardantoW(q,qp,Z,X,Z,O,W);
                                                     /* W from nautical angles and derivative */
#define nauttoW(q,qp,O,W)          cardantoW(q,qp,Z,Y,X,O,W);
                                                       /* T-B angles and derivative from m, W */
#define Wtocard(m,W,q1,q2,qp1,qp2)    Wtocardan(m,W,X,Y,Z,q1,q2,qp1,qp2);
                                                     /* Euler angles and derivative from m, W */
#define Wtoeul(m,W,q1,q2,qp1,qp2)     Wtocardan(m,W,Z,X,Z,q1,q2,qp1,qp2);
                                                  /* nautical angles and derivative from m, W */
#define Wtonaut(m,W,q1,q2,qp1,qp2)    Wtocardan(m,W,Z,Y,X,q1,q2,qp1,qp2);
#define cardanto_omega3(a,ap,i,j,k,m) cardanto_omega(a,ap,i,j,k,M m,3)
#define cardanto_omega4(a,ap,i,j,k,m) cardanto_omega(a,ap,i,j,k,M m,4)
                                            /* H from T-B angles, first and second derivative */
#define cardtoH(q,qp,qpp,O,H)         cardantoH(q,qp,qpp,X,Y,Z,O,H);
                                          /* H from Euler angles, first and second derivative */
#define eultoH(q,qp,qpp,O,H)          cardantoH(q,qp,qpp,Z,X,Z,O,H);
                                       /* H from nautical angles, first and second derivative */
#define nauttoH(q,qp,qpp,O,H)         cardantoH(q,qp,qpp,Z,Y,X,O,H);
                                      /* T-B angles, first and second derivative from m, W, H */
#define Htocard(m,W,H,q1,q2,qp1,qp2,qpp1,qpp2) Htocardan(m,W,H,X,Y,Z,q1,q2,qp1,qp2,qpp1,qpp2);
                                    /* Euler angles, first and second derivative from m, W, H */
#define Htoeul(m,W,H,q1,q2,qp1,qp2,qpp1,qpp2)  Htocardan(m,W,H,Z,X,Z,q1,q2,qp1,qp2,qpp1,qpp2);
                                 /* nautical angles, first and second derivative from m, W, H */
#define Htonaut(m,W,H,q1,q2,qp1,qp2,qpp1,qpp2) Htocardan(m,W,H,Z,Y,X,q1,q2,qp1,qp2,qpp1,qpp2);
#define cardanto_G3(a,ap,app,i,j,k,m) cardanto_G(a,ap,app,i,j,k,M m,3)
#define cardanto_G4(a,ap,app,i,j,k,m) cardanto_G(a,ap,app,i,j,k,M m,4)
                              /* --- Construction of Frames Attached to Points or Vectors --- */
#define frameP3(p1,p2,p3,a1,a2,r)     frameP(p1,p2,p3,a1,a2, M r,3)
#define frameP4(p1,p2,p3,a1,a2,r)     frameP(p1,p2,p3,a1,a2, M r,4)
#define frameV3(v1,v2,a1,a2,r)        frameV(v1,v2,a1,a2, M r,3)
#define frameV4(v1,v2,a1,a2,r)        frameV(v1,v2,a1,a2, M r,4)
                  /* --- Operations on Matrices and Vectors, Matrices and Vectors Algebra --- */
#define molt3(a,b,c)       molt(M a,M b,M c,3,3,3)  /* multiplication square matrices */
#define molt4(a,b,c)       molt(M a,M b,M c,4,4,4)
#define moltp(a,b,c)       molt(M a,M b,M c,4,4,1)  /* multiplication position matrix and point */
#define moltmv3(A,v1,v2)   molt(M A,M v1,M v2,3,3,1); /* multiplication 3x3 matrix and vector[3] */
#define rmolt3(a,r,b)      rmolt(M a,r,M b,3,3)     /* multiplication scalar and square matrix */
#define rmolt4(a,r,b)      rmolt(M a,r,M b,4,4)
#define rmoltv(v1,r,v2)    rmolt(M v1,r,M v2,3,1);  /* multiplication vector[3] by a scalar r */
#define sum3(a,b,c)        sum(M a,M b,M c,3,3)     /* sum of square matrices */
#define sum4(a,b,c)        sum(M a,M b,M c,4,4)
#define sumv(a,b,c)        sum(M a,M b,M c,3,1);    /* sum of vectors */
#define sub3(a,b,c)        sub(M a,M b,M c,3,3)     /* difference of square matrices */
#define sub4(a,b,c)        sub(M a,M b,M c,4,4)
#define subv(a,b,c)        sub(M a,M b,M c,3,1);    /* difference of vectors */
                   /* --- Operations on Matrices and Vectors, General Operations on Matrices --- */
#define clear3(m)          clear(M m,3,3)           /* clears a 3×3 matrix */
#define clear4(m)          clear(M m,4,4)           /* clears a 4×4 matrix */
#define idmat3(m)          idmat(M m,3)             /* makes unitary 3×3 matrices */
#define idmat4(m)          idmat(M m,4)             /* makes unitary 4×4 matrices */
#define transp3(a,b)       transp(M a,M b,3,3)      /* transpose of a matrix */
#define transp4(a,b)       transp(M a,M b,4,4)
                    /* --- Operations on Matrices and Vectors, General Operations on Vectors --- */
#define clearv(v)          clear(M v,3,1);          /* clears vector */
                                                                    /* --- Copy Functions --- */
#define mcopy3(a,b)        mcopy(M a,M b,3,3)       /* copy of square matrix */
#define mcopy4(a,b)        mcopy(M a,M b,4,4)
#define mcopy34(a,b)       mmcopy(M a,M b,3,4,3,3)  /* copy 3×3 submatrix */
#define mcopy43(a,b)       mmcopy(M a,M b,4,3,3,3)
#define vcopy(v1,v2)       mcopy(M v1,M v2,3,1);    /* copy vectors */
#define pcopy(P1,P2)       mcopy(M P1,M P2,4,1);    /* copy points */
                                                                   /* --- Print Functions --- */
#define printm(str,m)   fprintm3(stdout,str,m)      /* for compatibility */
#define printm3(str,m)  fprintm3(stdout,str,m)      /* print 3×3 matrices to stdout */
#define printm4(str,m)  fprintm4(stdout,str,m)      /* print 4×4 matrices to stdout */
#define printv(str,m,n) fprintv(stdout,str,m,n)     /* print vector to stdout */

        /* ============================ FUNCTIONS PROTOTYPES ====================== */

                                                    /* --- Rotation and Position Matrices --- */
int    dhtom(int jtype, real theta, real d, real b, real a, real alpha,real q, MAT4 m);
void   extract(MAT A, AXIS u, real *fi, int dim);
int    mtoscrew(MAT4 Q, AXIS u, real *fi, real *h, POINT P);
void   screwtom(AXIS u, real fi, real h, POINT P, MAT4 Q);
void   rotat(AXIS u, real fi, MAT A, int dim);                       /* build rotation matrix */
int    rotat2(int a, real q, MAT R, int dim);
int    rotat24(int a, real q, POINT O, MAT4 m);
void   rotat34(int a, real q, POINT P, MAT4 m);
                                                   /* --- Speed and Acceleration Matrices --- */
void   gtom(real gx, real gy, real gz, MAT4 Hg);
void   Gtomegapto (MAT3 G, VECTOR omegapto);/* extracts the angular velocity vector from the 3x3
                                          upper left submatrix G of the acceleration matrix H */
int    makeL(int jtype, AXIS u, real pitch, POINT P, MAT4 L);
int    makeL2(int jtype, int a, real pitch, POINT P, MAT4 L);
void   WtoL(MAT4 W, MAT4 L);                         /* W matrix to L matrix */
void   Wtovel(MAT4 W, AXIS u, real *omega, real *vel, POINT P);
int    velacctoWH(int jtype, real qp, real qpp, MAT4 W, MAT4 H);
int    velacctoWH2(int jtype, int a, real qp, real qpp, MAT4 W, MAT4 H);
                                      /* W and H matrices from axis, point and joint variable */
void   velacctoWH3(int jtype, int a, real qp, real qpp, POINT O, MAT4 W, MAT4 H);
                                                     /* --- Inertial and Actions Matrices --- */
void   actom(real fx, real fy, real fz, real cx, real cy, real cz, MAT4 FI);
void   jtoJ(real mass, real jxx, real jyy, real jzz, real jxy, real jyz,
            real jxz, real xg, real yg, real zg, MAT4 J);
                                             /* --- Matrix Transformations, Normalization --- */
int    normal(MAT R, int n);
int    norm_simm_skew(MAT A, int n, int dim, int sign);
                                       /* --- Matrix Transformations, Change of Reference --- */
void   trasf_mami(MAT4 A1, MAT4 m, MAT4 A2);
void   trasf_miam(MAT4 A1, MAT4 m, MAT4 A2);
void   trasf_mamt(MAT A1, MAT m, MAT A2, int dim);
void   trasf_miamit(MAT4 A1, MAT4 m, MAT4 A2);
                                        /* --- Matrix transformations, General Operations --- */
void   coriolis(MAT4 H0, MAT4 H1, MAT4 W0, MAT4 W1, MAT4 H);
void   inverse(MAT4 m, MAT4 mi);
void   mtov(MAT A, int dim, VECTOR v);      /* 3×3 skew matrix to vector[3] */
void   vtom(VECTOR v, MAT A, int dim);      /* vector[3] to 3×3 skew matrix */
void   skew(MAT A, MAT B, MAT C, int dim);
real   trac_ljlt4(MAT4 L1, MAT4 J, MAT4 L2);
                               /* --- Conversion from Cardan Angles to Matrices, Position --- */
int    cardantor(real *q, int i, int j, int k, MAT A, int dim);
int    rtocardan(MAT R, int dim, int i, int j, int k, real q1[3], real q2[3]);
void   cardantoM(real *q, int i, int j, int k, POINT O, MAT4 m);
                                              /* Euler/Cardan angles from a position matrix m */
int    Mtocardan(MAT4 m, int i, int j, int k, real q1[3], real q2[3]);
                          /* --- Conversion from Cardan Angles to Matrices, Velocity and Acceleration --- */
void   cardantoW(real *q, real *qp, int i, int j, int k, POINT O, MAT4 W);
    /* extracts the Euler/Cardan angles and their first time derivative from matrices m and W */
int    Wtocardan(MAT4 m, MAT4 W, int i, int j, int k, real q1[3], real q2[3], real qp1[3],
                 real qp2[3]);
void   cardanto_OMEGA(real *q, real *qp, int i, int j, int k, real *omega);
void   cardanto_omega(real *q, real *qp, int i, int j, int k, MAT A, int dim);
void   cardantoH(real *q, real *qp, real *qpp, int i, int j, int k, POINT O, MAT4 H);
                  /* extracts the Euler/Cardan angles and their first and second time derivative
                     from matrices m, W and H */
int    Htocardan(MAT4 m, MAT4 W, MAT4 H, int i, int j, int k, real q1[3], real q2[3],
                 real qp1[3], real qp2[3], real qpp1[3], real qpp2[3]);
void   cardanto_G(real *q, real *qp, real *qpp, int i, int j, int k, MAT A, int dim);
void   cardanto_OMEGAPTO(real *q, real *qp, real *qpp, int i, int j, int k, real *omegapto);
int    cardantol(real *q, int i, int j, int k, MAT R, int dim);
int    cardantoWPROD(real *q, int i, int j, int k, MAT R, int dim);
int    invA(real alpha, real beta, int sig, int i, int j, int k, MAT3 Ai);
void   WPRODtocardan(real alpha, real beta, int sig, int i, int j, int k, MAT3 Atil);
                              /* --- Construction of Frames Attached to Points or Vectors --- */
int    frameP(POINT P1, POINT P2, POINT P3, int a1, int a2, MAT A, int dim);
void   frame4P(POINT P1, POINT P2, POINT P3, int a1, int a2, MAT4 m);
int    frameV(VECTOR v1, VECTOR v2, int a1, int a2, MAT A, int dim);
void   frame4V(POINT P1, VECTOR v1, VECTOR v2, int a1, int a2, MAT4 m);
                        /* --- Working with Points Lines and Planes, Operations on Points --- */
real   angle(POINT P1, POINT P2, POINT P3);
real   dist(POINT P1, POINT P2);
void   intermediate(POINT P1, real a, POINT P2, real b, POINT P3);
void   middle(POINT P1, POINT P2, POINT P);                /* middle point with weights */
void   vect(POINT P1, POINT P2, VECTOR v);
                     /* --- Working with Points Lines and Planes, Op. on Lines and Planes --- */
void   line2p(POINT P1, POINT P2, LINEP l);                /* Builds line through P1,P2 */
void   linepvect(POINT P1, VECTOR v, LINEP l);    /* Builds line through P1, with direction v */
real   projponl(POINT P1, LINE l, POINT I);                /* projection of P on line l */
real   distpp(PLANE pl, POINT P);                          /* distance of P from plane */
real   project(POINT P, PLANE pl, POINT I);                /* projection of P on plane */
void   plane(POINT P1, POINT P2, POINT P3, PLANE pl);      /* plane through P1,P2,P3 */
void   plane2(POINT P1, VECTOR v, PLANE pl);         /* plane pl through P1 and unit vector v */
int    inters2pl(PLANE pl1, PLANE pl2, LINEP l);           /* intersection of two planes */
void   interslpl(LINE l, PLANE pl, POINT I, int *inttype); /*intersection of a line and plane */
void   intersection(LINE l1, LINE l2, LINEP lmindist, real * mindist,
                    PLANE pl, POINT I,int *inttype);       /* intersection of two lines */
                  /* --- Operations on Matrices and Vectors, Matrices and Vectors Algebra --- */
void   molt(MAT A, MAT B, MAT C, int d1, int d2, int d3);
void   rmolt(MAT A, real r, MAT B, int d1, int d2);
void   sum(MAT A, MAT B, MAT C, int d1, int d2);
void   sub(MAT A, MAT B, MAT C, int d1, int d2);
real   norm(MAT A, int d1, int d2);
                   /* --- Operations on Matrices and Vectors, General Operations on Matrices --- */
void   clear(MAT A, int id, int jd);
void   idmat(MAT A, int ijd);
void   transp(MAT A, MAT At, int d1, int d2);
int    pseudo_inv(MAT A, MAT Api, int rows, int cols);
                    /* --- Operations on Matrices and Vectors, General Operations on Vectors --- */
void   cross(VECTOR a, VECTOR b, VECTOR  c);   /* cross product c = a × b */
real   dot(VECTOR a, VECTOR b);                /* dot product of 3 element vector */
real   dot2(real *v1, real *v2, int dim);      /* dot product of dim element vector */
real   mod(VECTOR a);                          /* module of a vector */
real   unitv(VECTOR v, AXIS u);                /* evaluate the unit vector u of a vector v
                                                 and it also returns its module */
void   vector(AXIS u, real m, VECTOR v);       /* evaluate a vector v whose module is
                                                  m and whose unit vector is u */
                                          /* --- Copy Functions --- */
void   mcopy(MAT A1, MAT A2, int d1, int d2);
void   mmcopy(MAT A, MAT B, int d1, int d2, int im, int jm);
int    mvcopy(MAT A, int rows, int cols, int val, int type, real *v); /* from row/column of A
                                                                         to vector */
int    vmcopy(real *v, int dim, int val, int type, MAT A, int rows,
              int cols);                       /* from vector to row/column of A */
                                               /* --- Print Functions --- */
void   fprintm3(FILE *out, char *s, MAT3 A);
void   fprintm4(FILE *out, char *s, MAT4 A);
void   fprintv(FILE *out, char *s, real *v, int n);
void   prmat(FILE *grpout,char str[], MAT4 m);
void   printmat(char *s, MAT A, int idim, int jdim, int imax, int jmax);
void   iprintmat(char *s, int *A, int idim, int jdim, int imax, int jmax);
                                          /* --- Machine Precision Functions --- */
real   zerom(void);                            /* real machine precision */
float  fzerom(void);                           /* float machine precision */
double dzerom(void);                           /* double machine precision */

void   crossVtoM(real *a, real *b, MAT C, int dim);
void   crossMtoM(MAT A, MAT B, MAT C, int dima, int dimb, int dimc);
void   axis(int a, AXIS A);
void   traslat(VECTOR u, real h, MAT4 m);
void   traslat2(int a, real h, MAT4 m);
void   traslat24(int a, real h, POINT p, MAT4 m);
real   psedot(MAT4 A, MAT4 B);

void runtime_error(char *name, char *msg, int ex);
#endif