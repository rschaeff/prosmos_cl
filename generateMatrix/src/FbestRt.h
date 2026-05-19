
/* this is file is the header file for bestrt.c*/
#ifndef _BESTRT_H
#define _BESTRT_H
double* generatedIdealHelix(int npoints);
double measure_distance ( double ix, double iy, double iz, double jx, double jy, double jz );
double angle(double mix,double miy,double miz,double ix,double iy,double iz,double mjx,double mjy,double mjz,double jx,double jy,double jz);
double angle_between_lines(double* xpoint1,double* xpoint2,double* ypoint1,double* ypoint2);
double **bestRT(double **a, double **b, int n );

#endif

