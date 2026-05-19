
/* this is file is the header file for rotfit.c*/
#ifndef _ROTFIT_H
#define _ROTFIT_H
#include "external.h"
#define NUMBER 6000
void rotFit(double *data,int npoints, int dim,VECTOR* retData);
void pointProjToVect(double* data,int npoints,VECTOR vect,point* pP);

#endif


