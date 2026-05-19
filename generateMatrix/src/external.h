#ifndef EXTERNAL_H
#define EXTERNAL_H
struct point
{
        double   xcoord;
        double   ycoord;
        double   zcoord;
};
typedef struct point point;
struct _vector
{
        point initPoint;
        point endPoint;
};
typedef struct _vector myvector;
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1):(dmaxarg2))
static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1):(dminarg2))
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1):(imaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1):(iminarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define fabslmy(a) ((a) < 0.0 ? -a : a)
#define SIGNL(a,b) ((b) >= 0.0 ? fabslmy(a) : -fabslmy(a))
#ifndef M_PI
#define M_PI    3.1415926
#endif
#define  DISTANCE_DEFAULT    11.0 // 8 A this parament need to be adjusted
#define  OVERLAP_DEFAULT     2.5  // 2 A, adjusted to real world
#define  ANGLE_DEFAULT       85.0
#define  SENDSIZE 7000
#define  FILENUMBER 500
static double **dmatrix(long nrl, long nrh, long ncl, long nch);
static void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch);
static void eigensystem_symm(double **a, long n, double *eval, double 
**evec, double *helper);
static void tred2(double **a, int n, double d[], double e[]);
static void tqli(double d[], double e[], int n, double **z);
static double dpythag(double a, double b);
double measure_distance ( double ix, double iy, double iz, double jx, double jy, double jz );
double angle ( double mix, double miy, double miz, double ix, double iy, double iz,
               double mjx, double mjy, double mjz, double jx, double jy, double jz )
{
        double dot,ma,mb;
        double  cos_theta,theta;
  

// done by a dot product 
 //  * a.b = |a||b|cos_theta 
        dot = ( ix - mix ) * ( jx - mjx ) + ( iy - miy ) * ( jy - mjy ) 
                + ( iz - miz ) * ( jz - mjz );

        ma = measure_distance ( mix, miy, miz, ix, iy, iz );
        mb = measure_distance ( mjx, mjy, mjz, jx, jy, jz );

        cos_theta = dot / ( ma * mb ); 
  // cos_theta can only take values from -1 to 1. remove floating
   // point errors beyond 4th place of decimal 
        if ( ( cos_theta > 1.01 )|| ( cos_theta < -1.01 ) )
        {       
                fprintf(stderr, "ERROR: the cos_theta is over the range from 1 to -1 \n");
        
        }
        if ( cos_theta > 1 )cos_theta = 1;
        if ( cos_theta < -1 )
                cos_theta = -1;
  // convert to degrees and return value. return is 0 to 180 
        theta = fabs ( acos ( cos_theta ) / M_PI * 180 );
        return theta;
}
double **bestRT(double **a, double **b, int n );

double angle_between_lines(double* xpoint1,double* xpoint2,double* ypoint1,double* ypoint2);

void pointProjToVect(double* data,int npoints,myvector vect,point* pP);
void rotFit(double *data, int npoints, int dim,myvector* retVect)
{
        double **matrix1;
        double *data2;
        double *data3;
        int temPoint;
        double **a;
        double **b;
        double **c;
        double firstP[3],massCenter[3],lastP[3], vertualP[3];
        double retData[6];
        double  theta1, theta2;
        // store the angel between line
        double length1=0,length2,t1=0,t2=0;
        double  ca[3];
        int i=0, j=0, k=0;
        if(data == NULL )
        {
                printf("No Data inputed \n");
                return;
        }
        if( retData == NULL)
        {
                printf(" No valid return pointer \n");
                return;
        }

        a = dmatrix(1,(npoints-1),1,dim);
        b = dmatrix(1,(npoints-1),1,dim);
        c = dmatrix(1,(npoints-1),1,dim);

        // calculate the mass center of helix or strand

        for(i= 0; i<3; i++)
        {
                ca[i]=0;
                for(j=0;j<npoints;j++)
                {
                        ca[i] =ca[i]+data[dim*j+i];
                }
                ca[i] = ca[i]/npoints;
        }


        data2 = data+1*dim;
        data3 = data2+1*dim;
        temPoint = npoints;
        for(i=1;i<temPoint;i++)
                for(j = 1; j<= dim; j++)
                {
                        a[i][j] = *(data+(i-1)*dim +(j-1));
                        b[i][j] = *(data + i*dim + (j-1));
                        if(i<(temPoint-1))
                        {
                                c[i][j] = *(data+(i+1)*dim +(j-1));
                        }
                        else
                        {
                                c[i][j] = 0.0;
                        }
                }

        firstP[0] = a[1][1];
        firstP[1] = a[1][2];
        firstP[2] = a[1][3];
        lastP[0] = data[(npoints-1)*dim];
        lastP[1] = data[(npoints-1)*dim+1];
        lastP[2] = data[(npoints-1)*dim+2];

        matrix1 = bestRT(a,b,temPoint-1);
        //get a point of rotational axis
        //  matrix2 = bestRT(b, c,(temPoint-2)); get another point of rotational axis
        for(j=0;j<6;j++)
        {

                if(j<3)
                        retData[j] = matrix1[5][j+1];// get the unite rotational axis vector
                else
                        retData[j] = matrix1[6][j-2];// get the center of mass of helix or strand
        }
         massCenter[0] = ca[0];
         massCenter[1] = ca[1];
         massCenter[2] = ca[2];
        //massCenter[0] = retData[3]; Center of Mass of Helix X coordinate 
        //massCenter[1] = retData[4]; Center of Mass of Helix Y coordinate
        //massCenter[2] = retData[5];Center of Mass of Helix Z coordinate

        //calculate another point on the rotational axis based on center of mass of helix, the formula is                //  x= x0 + tX, y = y0 + tY, z = z0 + tZ     t is any real number  

         if((retData[0] != 0.0) || (retData[1] != 0.0) || (retData[2] != 0.0))
         {

                 for(j=0;j<3;j++)
                 {
                        vertualP[j] =massCenter[j] + 5* retData[j]; // t= 1
                 }

                 theta1 = angle(massCenter[0],massCenter[1],massCenter[2],firstP[0],firstP[1],firstP[2],
                        massCenter[0],massCenter[1],massCenter[2],vertualP[0],vertualP[1],vertualP[2]);
                 theta2 = angle(massCenter[0],massCenter[1],massCenter[2],lastP[0],lastP[1],lastP[2],
                        massCenter[0],massCenter[1],massCenter[2],vertualP[0],vertualP[1],vertualP[2]);
// I will figure out the projected points of the first point and last point, using them to stand vector

                 length1 =  measure_distance(firstP[0],firstP[1],firstP[2],massCenter[0],massCenter[1],massCenter[2]);
// the length between the first point of helix and mass center
             length2 =  measure_distance(lastP[0],lastP[1],lastP[2],massCenter[0],massCenter[1],massCenter[2]);
// the length between the last point of helix and mass center
             t1 =  (length1 * (cos(theta1/180*M_PI)));
// t1 is negative because of first points is in the opposite direction of rotational vector
             t2 =  (length2 * (cos(theta2/180*M_PI)));
             for(j=0; j<3; j++)
                 {

                     retData[j] =massCenter[j] + t1* matrix1[5][j+1];
                     retData[3+j] = massCenter[j] + t2* matrix1[5][j+1];
                 }

             retVect->initPoint.xcoord = retData[0];
             retVect->initPoint.ycoord =retData[1];
             retVect->initPoint.zcoord =retData[2];
             retVect->endPoint.xcoord =retData[3];
             retVect->endPoint.ycoord =retData[4];
             retVect->endPoint.zcoord =retData[5];
         }
         else
         {
             retVect->initPoint.xcoord = firstP[0];
             retVect->initPoint.ycoord = firstP[1];
             retVect->initPoint.zcoord = firstP[2];
             retVect->endPoint.xcoord = lastP[0];
             retVect->endPoint.ycoord = lastP[1];
             retVect->endPoint.zcoord = lastP[2];
         }

        free_dmatrix(a,1,(npoints-1),1,dim);
        free_dmatrix(b,1,(npoints-1),1,dim);
        free_dmatrix(c,1,(npoints-1),1,dim);
        free_dmatrix(matrix1,1,6,1,3);

        return;
       // two point is in the same rotational axis, so decide the vector
}

		
    

double **bestRT(double **a, double **b, int n )
{
	double **m, **s, **evec;
	double helper[5];
    double eval[5];
    double x[5],ca[4],cb[4];
    double c;
    int i,j,t;
        
    m= dmatrix(1,6,1,3);
	s= dmatrix(1,6,1,6);
	evec = dmatrix(1,6,1,6);
	/* calculate the center of mass coordinates for both sets */
	for(j=1;j<=3;j++)
	{
		ca[j]=0.0;cb[j]=0.0;
        for(i=1;i<=n;i++)
		{
			ca[j]=ca[j]+a[i][j];
			cb[j]=cb[j]+b[i][j];
		}
		ca[j]=ca[j]/n;
		cb[j]=cb[j]/n;
	}
	/* calculate the magic matrix s */
    for(i=1;i<=3;i++)
		for(j=1;j<=3;j++)
		{
			m[i][j]=0.0;
			for(t=1;t<=n;t++)
				m[i][j]=m[i][j]+(a[t][i]-ca[i])*(b[t][j]-cb[j]);
		}
    s[1][1]=m[1][1] - m[2][2] - m[3][3];
    s[1][2]=s[2][1]=m[1][2] + m[2][1];
    s[1][3]=s[3][1]=m[1][3] + m[3][1];
    s[1][4]=s[4][1]=m[2][3] - m[3][2];
    s[2][1]=s[1][2]=m[1][2] + m[2][1];
    s[2][2]=-m[1][1] + m[2][2] - m[3][3];
    s[2][3]=s[3][2]=m[2][3] + m[3][2];
    s[2][4]=s[4][2]=-m[1][3] + m[3][1];
    s[3][1]=s[1][3]=m[1][3] + m[3][1];
    s[3][2]=s[2][3]=m[2][3] + m[3][2];
    s[3][3]=-m[1][1] - m[2][2] + m[3][3];    
    s[3][4]=s[4][3]= m[1][2] - m[2][1];
    s[4][1]=s[1][4]= m[2][3] - m[3][2];
    s[4][2]=s[2][4]=-m[1][3] + m[3][1];
    s[4][3]=s[3][4]=m[1][2] - m[2][1];
    s[4][4]=m[1][1] + m[2][2] + m[3][3];
	/* solve eigensystem for the matrix s*/
	eigensystem_symm(s, 4, eval, evec, helper );
	/* find maximum eigenvalue */
	t=1;
	for(i=2;i<=4;i++)
	{
		if(eval[t]<eval[i])t=i;
	}

 /* take the eigenvector corresponding to the largest eigenvalue */

    for(i=1;i<=4;++i)x[i]=evec[i][t];

 /* calculate the rotation matrix */

    m[1][1]=x[1]*x[1] - x[2]*x[2] - x[3]*x[3] + x[4]*x[4];
    m[1][2]=2*x[1]*x[2] - 2*x[3]*x[4];
    m[1][3]=2*x[1]*x[3] + 2*x[2]*x[4];
    m[2][1]=2*x[1]*x[2] + 2*x[3]*x[4];
    m[2][2]=-x[1]*x[1] + x[2]*x[2] -  x[3]*x[3] + x[4]*x[4];
    m[2][3]=2*x[2]*x[3] - 2*x[1]*x[4];
    m[3][1]=2*x[1]*x[3] - 2*x[2]*x[4];
    m[3][2]=2*x[2]*x[3] + 2*x[1]*x[4];
    m[3][3]=-x[1]*x[1] - x[2]*x[2] +  x[3]*x[3] + x[4]*x[4];

 /* calculate the rotation axis */

    c=(m[1][1]+m[2][2]+m[3][3]-1.0)/2.0;
	c = c*c;
    if((double)c>=1){
        fprintf(stderr,"Matrix is not orthogonal in bestR():FATAL\n\n");
	exit(0);}
    c=2.0*sqrt(1-c);
    m[5][1]=(m[3][2]-m[2][3])/c;
    m[5][2]=(m[1][3]-m[3][1])/c;
    m[5][3]=(m[2][1]-m[1][2])/c;

 /* calculate the translation vector */

    for(j=1;j<=3;j++)
	{
		m[4][j]=cb[j];
		for(i=1;i<=3;i++)
			m[4][j]=m[4][j]-m[j][i]*ca[i];
	}
	/* store the center of mass coordinates  for data **a */
	for(j=1; j<=3;j++)
	{
		m[6][j] = ca[j];
	}

        free_dmatrix(s,1,6,1,6);
	free_dmatrix(evec,1,6,1,6);
	return m;
}

static double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	m=(double **)malloc((size_t)((nrow+1)*sizeof(double*)));
	if (!m) 
	{
		fprintf(stderr,"allocation failure 1 in dmat()\n");
		exit(0);
	}
	m += 1;
	m -= nrl;
	m[nrl]=(double *)malloc((size_t)((nrow*ncol+1)*sizeof(double)));
	if (!m[nrl]) 
	{
		fprintf(stderr,"allocation failure 2 in dmat()\n");
		exit(0);
	}
	m[nrl] += 1;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return m;
}


static void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch)
{
	free((char *) (m[nrl]+ncl-1));
	free((char *) (m + nrl -1));
}
static void eigensystem_symm(double **a, long n, double *eval, double **evec, double *helper)
/* solves eigensystem for symmetric real matrix a[1.n][1.n]
stores eigenvalues in eval[1..n] and correspoding vectors as columns in evec[1..n][1..n]
helper[1..n] is a storage vector to save time for its allocation
*/
{
	int i,j;
	for(i=1;i<=n;++i)
	{
		for(j=1;j<=n;++j)
		{
			evec[i][j]=a[i][j];
		}
	}
tred2(evec,n,eval,helper);
tqli(eval,helper,n,evec);
}


static void tred2(double **a, int n, double d[], double e[])
/* reduces a real symmetric matrix a[1.n][1.n] to tridiagonal form
a is replaced by orthoginal matrix, necessary further
d[1..n] - diagonal elements of the tridiagonal matrix;
e[1..n] - off-diagonal elements with e[1]=0
*/

{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

static void tqli(double d[], double e[], int n, double **z)
/*
determines eigenvalues and eigenvectors of real symmetric tridiagonal matrix:
on input d[1..n] dagonal elements of tridiagonal matrix; on output - eigenvalues;
e[1..n]
the k-th column of z returns eigenvector corresponding to d[k];
on input z is matrix, outputed by tred2 
*/

{
	int m,l,iter,i,k;
	double  s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break; 
			/*	if (fabs(e[m]) <1.0e-14) break;  */

			}
			if (m != l) {
				if (iter++ == 30) fprintf(stdout,"Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=dpythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=dpythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=1;k<=n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}


static double dpythag(double a, double b)
/* computes (a^2+b^2)^0.5 without underflow or overflow */
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

double angle_between_lines(double* xpoint1,double* xpoint2,double* ypoint1,double* ypoint2)
{
	double mix,  miy,  miz;
	double ix,  iy, iz;
	double mjx, mjy, mjz;
	double jx, jy,  jz; 
	double cos_theta, theta;
    double dot, ma, mb;
	mix = xpoint1[0];
	miy = xpoint1[1];
	miz = xpoint1[2];
	ix = xpoint2[0];
	iy = xpoint2[1];
	iz = xpoint2[2];
	mjx = ypoint1[0];
	mjy = ypoint1[1];
	mjz = ypoint1[2];
	jx = ypoint2[0];
	jy = ypoint2[1];
	jz = ypoint2[2];
	
  /* done by a dot product
   * a.b = |a||b|cos_theta */
	dot = ( ix - mix ) * ( jx - mjx ) + ( iy - miy ) * ( jy - mjy )
		+ ( iz - miz ) * ( jz - mjz );
	ma = measure_distance ( mix, miy, miz, ix, iy, iz );
	mb = measure_distance ( mjx, mjy, mjz, jx, jy, jz );
	cos_theta = dot / ( ma * mb );
        /* cos_theta can only take values from -1 to 1. remove floating
        * point errors beyond 4th place of decimal */
	if ( ( cos_theta > 1.01 )|| ( cos_theta < -1.01 ) )
	{
		fprintf(stderr, "ERROR: the cos_theta is over the range from 1 to -1 \n");
/*		return 400;*/
		
	}
	if ( cos_theta > 1 )cos_theta = 1;
	if ( cos_theta < -1 )
		cos_theta = -1;
  /* convert to degrees and return value. return is 0 to 180 */
	theta = fabs ( acos ( cos_theta ) / M_PI * 180 );
  /* angle can only be between 0 to 90 */
	if ( theta > 90 )theta = 180 - theta;
  /* check if theta falls within allowed range */
  	return theta;
}
double measure_distance ( double ix, double iy, double iz, double jx, double jy, double jz )
{
	return sqrt (((jx-ix)*(jx-ix))+((jy-iy)*(jy-iy))+((jz-iz)*(jz-iz)));
}
    


double* generatedIdealHelix(int npoints)/* this function is about to generate the ideal helix along z-coord, the diameter is about 4 and up 2 for every circle*/
{
	int i;
	double * retP = (double*)malloc(sizeof(double)*npoints*3);
    if( retP == NULL) return retP;

	for(i=0;i<npoints;i++)
	{
		retP[3*i] =2*cos(M_PI/2*i);
		retP[3*i+1] = 2*sin(M_PI/2*i);
		retP[3*i+2] = i/2;
	}
	return retP;
}

void pointProjToVect(double* data,int npoints,myvector vect,point* pP)
{
        int i;
        double length, theta,t1;
        double tempP[3];
        double firstP[3],endP[3];
        point dirV;
        firstP[0] = vect.initPoint.xcoord;
        firstP[1] = vect.initPoint.ycoord;
        firstP[2] = vect.initPoint.zcoord;
        endP[0] = vect.endPoint.xcoord;
        endP[1] = vect.endPoint.ycoord;
        endP[2] = vect.endPoint.zcoord;
        dirV.xcoord = vect.endPoint.xcoord - vect.initPoint.xcoord;
        dirV.ycoord = vect.endPoint.ycoord - vect.initPoint.ycoord;
        dirV.zcoord = vect.endPoint.zcoord - vect.initPoint.zcoord;
        length = sqrt((dirV.xcoord*dirV.xcoord)+(dirV.ycoord*dirV.ycoord)+(dirV.zcoord*dirV.zcoord));
        dirV.xcoord = dirV.xcoord/length;
        dirV.ycoord = dirV.ycoord/length;
        dirV.zcoord = dirV.zcoord/length;
        for(i = 0; i<npoints; i++)
        {
                tempP[0] = data[3*i];
                tempP[1] = data[3*i+1];
                tempP[2] = data[3*i+2];
                theta = angle(firstP[0],firstP[1],firstP[2],tempP[0],tempP[1],tempP[2],
                        firstP[0],firstP[1],firstP[2],endP[0],endP[1],endP[2]);
                if(theta >180.001)
                        return;
                length =  measure_distance(firstP[0],firstP[1],firstP[2],tempP[0],tempP[1],tempP[2]);
                t1 =  (length * (cos(theta/180*M_PI)));
                pP[i].xcoord = firstP[0] + t1 *dirV.xcoord;
                pP[i].ycoord = firstP[1] + t1 *dirV.ycoord;
                pP[i].zcoord = firstP[2] + t1 *dirV.zcoord;
        }
        return ;
}
double overlapBetweenVector (point vector1point1, point vector1point2, point vector2point1, point vector2point2)
{
        point Mpoint1,Mpoint2;
        point Npoint1,Npoint2;
        point dirvector1,dirvector2;
        double length, t1,t2,t3,t4,temp;
        double overlap1=0,overlap2=0,retdata=0;
        /*calcute the first set of two midpoints of vector1 to vector2*/
        Mpoint1.xcoord = (vector1point1.xcoord+vector2point1.xcoord)/2;
        Mpoint1.ycoord = (vector1point1.ycoord+vector2point1.ycoord)/2;
        Mpoint1.zcoord = (vector1point1.zcoord+vector2point1.zcoord)/2;

        Mpoint2.xcoord = (vector1point2.xcoord+vector2point2.xcoord)/2;
        Mpoint2.ycoord = (vector1point2.ycoord+vector2point2.ycoord)/2;
        Mpoint2.zcoord = (vector1point2.zcoord+vector2point2.zcoord)/2;
        length = sqrt((Mpoint2.xcoord-Mpoint1.xcoord)*(Mpoint2.xcoord-Mpoint1.xcoord)+
                (Mpoint2.ycoord-Mpoint1.ycoord)*(Mpoint2.ycoord-Mpoint1.ycoord)+
                (Mpoint2.zcoord-Mpoint1.zcoord)*(Mpoint2.zcoord-Mpoint1.zcoord));

        dirvector1.xcoord = (Mpoint2.xcoord-Mpoint1.xcoord)/length;
        dirvector1.ycoord = (Mpoint2.ycoord-Mpoint1.ycoord)/length;
        dirvector1.zcoord = (Mpoint2.zcoord-Mpoint1.zcoord)/length;

        t1 = (vector1point1.xcoord-Mpoint1.xcoord)*dirvector1.xcoord +
             (vector1point1.ycoord-Mpoint1.ycoord)*dirvector1.ycoord +
             (vector1point1.zcoord-Mpoint1.zcoord)*dirvector1.zcoord;

        t2 = (vector1point2.xcoord-Mpoint1.xcoord)*dirvector1.xcoord +
             (vector1point2.ycoord-Mpoint1.ycoord)*dirvector1.ycoord +
             (vector1point2.zcoord-Mpoint1.zcoord)*dirvector1.zcoord;

        t3 = (vector2point1.xcoord-Mpoint1.xcoord)*dirvector1.xcoord +
             (vector2point1.ycoord-Mpoint1.ycoord)*dirvector1.ycoord +
             (vector2point1.zcoord-Mpoint1.zcoord)*dirvector1.zcoord;

        t4 = (vector2point2.xcoord-Mpoint1.xcoord)*dirvector1.xcoord +
             (vector2point2.ycoord-Mpoint1.ycoord)*dirvector1.ycoord +
             (vector2point2.zcoord-Mpoint1.zcoord)*dirvector1.zcoord;

        if(t1 >t2)
        {
                temp = t1;
                t1= t2;
                t2 = temp;
        }
        if( t1 == t2 ) overlap1 = 0; /* just one point, no overlap*/
        else
        {
                if((t3 >=t1) & (t3 <=t2))/* t3 is between t1 and t2*/
                {
                        if(t4<=t1)/* the ovelap is between t1-t3*/
                                overlap1 = t3-t1;
                        else if(t4>=t2)/* the overlap is between t3-t2*/
                                overlap1 = t2-t3;
                        else       /* the overlap is between t3-t4*/
                        {
                                overlap1 = t4-t3;
                                if(overlap1 < 0)
                                        overlap1 = t3-t4;
                        }
                }
                else if(t3 < t1)
                {

                        if(t4<=t1)/* the ovelap is 0*/
                                overlap1 = 0;
                        else if(t4>=t2)/* the overlap is between t1-t2*/
                                overlap1 = t2-t1;
                        else       /* the overlap is between t1-t4*/
                                overlap1 = t4-t1;
                }
                else /* t3>t2*/
                {
                        if(t4<=t1)/* the ovelap is t1-t2*/
                                overlap1 = t2-t1;
                        else if(t4>=t2)/* the overlap is 0*/
                                overlap1 = 0;
                        else       /* the overlap is between t2-t4*/
                                overlap1 = (t2-t4);
                }
        }

        /*calcute the second set of two midpoints of vector1 to vector2*/
        Npoint1.xcoord = (vector1point1.xcoord+vector2point2.xcoord)/2;
        Npoint1.ycoord = (vector1point1.ycoord+vector2point2.ycoord)/2;
        Npoint1.zcoord = (vector1point1.zcoord+vector2point2.zcoord)/2;

        Npoint2.xcoord = (vector1point2.xcoord+vector2point1.xcoord)/2;
        Npoint2.ycoord = (vector1point2.ycoord+vector2point1.ycoord)/2;
        Npoint2.zcoord = (vector1point2.zcoord+vector2point1.zcoord)/2;
        length = sqrt((Npoint2.xcoord-Npoint1.xcoord)*(Npoint2.xcoord-Npoint1.xcoord)+
                (Npoint2.ycoord-Npoint1.ycoord)*(Npoint2.ycoord-Npoint1.ycoord)+
                (Npoint2.zcoord-Npoint1.zcoord)*(Npoint2.zcoord-Npoint1.zcoord));

        dirvector2.xcoord = (Npoint2.xcoord-Npoint1.xcoord)/length;
        dirvector2.ycoord = (Npoint2.ycoord-Npoint1.ycoord)/length;
        dirvector2.zcoord = (Npoint2.zcoord-Npoint1.zcoord)/length;
        t1 = (vector1point1.xcoord-Npoint1.xcoord)*dirvector2.xcoord +
                (vector1point1.ycoord-Npoint1.ycoord)*dirvector2.ycoord +
                (vector1point1.zcoord-Npoint1.zcoord)*dirvector2.zcoord;
        t2 = (vector1point2.xcoord-Npoint1.xcoord)*dirvector2.xcoord +
                (vector1point2.ycoord-Npoint1.ycoord)*dirvector2.ycoord +
                (vector1point2.zcoord-Npoint1.zcoord)*dirvector2.zcoord;
        t3 = (vector2point1.xcoord-Npoint1.xcoord)*dirvector2.xcoord +
                (vector2point1.ycoord-Npoint1.ycoord)*dirvector2.ycoord +
                (vector2point1.zcoord-Npoint1.zcoord)*dirvector2.zcoord;
        t4 = (vector2point2.xcoord-Npoint1.xcoord)*dirvector2.xcoord +
                (vector2point2.ycoord-Npoint1.ycoord)*dirvector2.ycoord +
                (vector2point2.zcoord-Npoint1.zcoord)*dirvector2.zcoord;
        if(t1 >t2)
        {
                temp = t1;
                t1= t2;
                t2 = temp;
        }
        if( t1 == t2 ) overlap2 = 0; /* just one point, no overlap*/
        else
        {
                if((t3 >=t1) & (t3 <=t2))/* t3 is between t1 and t2*/
                {
                        if(t4<=t1)/* the ovelap is between t1-t3*/
                                overlap2 = t3-t1;
                        else if(t4>=t2)/* the overlap is between t3-t2*/
                                overlap2 = t2-t3;
                        else       /* the overlap is between t3-t4*/
                        {
                                overlap2 = t4-t3;
                                if(overlap2 < 0)
                                        overlap2 = t3-t4;
                        }
                }
                else if(t3 < t1)
                {

                        if(t4<=t1)/* the ovelap is 0*/
                                overlap2 = 0;
                        else if(t4>=t2)/* the overlap is between t1-t2*/
                                overlap2 = t2-t1;
                        else       /* the overlap is between t1-t4*/
                                overlap2 = t4-t1;
                }
                else /* t3>t2*/
                {
                        if(t4<=t1)/* the ovelap is t1-t2*/
                                overlap2 = t2-t1;
                        else if(t4>=t2)/* the overlap is 0*/
                                overlap2 = 0;
                        else       /* the overlap is between t2-t4*/
                                overlap2 = (t2-t4);
                }
        }
        if(overlap1>=overlap2) retdata = overlap1;
        else
                retdata = overlap2;

        return  retdata;
}	
double distanceBetweenTwoPoint(point mpoint,point npoint)
{
        return sqrt((mpoint.xcoord-npoint.xcoord)*(mpoint.xcoord-npoint.xcoord) + (mpoint.ycoord-npoint.ycoord)*(
mpoint.ycoord-npoint.ycoord) + (mpoint.zcoord-npoint.zcoord)*(mpoint.zcoord-npoint.zcoord));
}


double distanceBetweenTwoElm(point* kElm,int k, point* lElm, int l)
{
        double prevdis=1000.0,currdis = 1000.0;
        int i,j;
        if((kElm == NULL) || (lElm == NULL) || (k==0) || (l == 0))
        {
                fprintf(stdout, " NOT VALID DATA \n");
                fprintf(stdout,"kElm = %d, lElm = %d, k = %d, l= %d\n",kElm, lElm, k,l);
                return -1;
        }
        for(i=0;i<k;i++)
        {
                for(j=0;j<l;j++)
                {
                        currdis = distanceBetweenTwoPoint(kElm[i],lElm[j]);
                        if(currdis <prevdis)
                                prevdis = currdis;
                }
        }
        return prevdis;
}
#endif

