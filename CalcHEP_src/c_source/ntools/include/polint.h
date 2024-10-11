#ifndef __POLINT__
#define __POLINT__

extern double polint3(double x, int dim, const  double *xa, const double *ya);
extern double polint1(double x, int n,  double *xa, double *ya);
extern double polint1Exp(double x, int n,  double *xa, double *ya);
extern double polintN(double x, int n, const  double *xa, const double *ya);
extern void   polintDiff(int n, const  double *xa, const double *ya, double * dxdy);

typedef struct { int dim; double *x; double *y;}   polintStr;
extern double polint_arg(double x, void * arg);
#endif
