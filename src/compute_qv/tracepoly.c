/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

/* 
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 *
 * 2.7 2003/11/06 18:18:49
 *
 * Produce a polynomial fit to a trace.
 * The underlying numerical code is from Numerical Recipes
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>

#include "Btk_qv.h"
#include "tracepoly.h"



/*******************************************************************************
 * fudge factors, it should be OK to tweak these 
 ******************************************************************************/
const double
    DEFAULT_CONST_POLY  = 500.0,
    CUTOFF_FACTOR       = 0.5,      /* range=(0,1) */
    CUTOFF_NUM_SAMPLES  = 100.0,
    IIR_CONST           = 0.9,      /* Infinite Impulse Response; range=(0,1) */
    NORMALIZATION_SLACK = 16.0,  /* NOTICE */
    MAX_POLY_LOG_Y      = 7.3777589,/* log(1600) */
    MIN_POLY_Y          = 1.0,
    POINTS_PER_COEF     = 5.0;

#define USE_SVD 0       /* svdfit is slower but more robust */


/*******************************************************************************
 * various debug stuff 
 ******************************************************************************/
#define POLY_DBG 0       /* debug for fitting code */

/* Global variables to get useful debug info to the routines below */
#if TELL_ME_IF_I_EXCEED_THE_LIMITS      /* tracepoly.h */
char* Filename;
int ReadNum, BadReadNum, BaseIndx, ReportNum;
#define IFR if(0) /* if( BaseIndx==2 ) */
#else
#define IFR if(0) /* if( ReadNum==BadReadNum && BaseIndx == 1 ) */
#endif

#define POLY_EVEN 0
#define POLY_ALL  2
#define POLY_FLAG POLY_EVEN     /* don't change this */


/*******************************************************************************
 * typedefs
 ******************************************************************************/
typedef struct {
    int poly_type;
    double xmin, xmax;
} FuncArgs;

typedef void MyFunc( FuncArgs*, double, double[], int );


#define MYMIN(A,B)      ( (A)<(B) ? (A) : (B) )
#define MYMAX(A,B)      ( (A)>(B) ? (A) : (B) )

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


/*******************************************************************************
 * The next 3 routines allow evaluation of Legendre Polynomials 
 ******************************************************************************/

/* This is "fleg()" from Numerical Recipes (NR) 
 * Slightly modified by adding argument of type FuncArgs and using it to do
 * normalization of argument "x".
 */
static void 
func_legendre_1based_sub( FuncArgs* fa, double x, double pl[], int nl )
{
    int j;
    double twox,f2,f1,d;
    double delta = fa->xmax - fa->xmin;

    x = (delta==0.0 ? 0.0 : (2.0*x-fa->xmin-fa->xmax)/delta );

    pl[1]=1.0;
    if( nl==1 ) return; /* added by Scott Tyler to handle constant poly */
    pl[2]=x;
    if (nl > 2) {
        twox=2.0*x;
        f2=x;
        d=1.0;
        for (j=3;j<=nl;j++) {
            f1=d++;
            f2 += twox;
            pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
        }
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

/* modification to "fleg()" to allow only even powers of x if desired */
static void 
func_legendre_1based( FuncArgs* fa, double x, double pl1[], int nl )
{
    int n, i;
    double  *pl0 = pl1+1, ptemp0[2*NumCoef-1];

    switch( fa->poly_type ) {
    case POLY_EVEN:
        n = 2*nl-1;
        func_legendre_1based_sub( fa, x, ptemp0-1, n );
        for( i=0; i<nl; i++ ) { pl0[i] = ptemp0[i*2]; }
        break;
    case POLY_ALL:
        func_legendre_1based_sub( fa, x, pl1, nl );
        break;
    default:
        fprintf( stderr, "illegal poly_type=%d\n", fa->poly_type );
        exit(-1);
    }
}

/* interface to "fleg()" to allow 0-based arrays instead of NR type arrays */
static void 
func_legendre_0based( FuncArgs* fa, double x, double pl[], int nl )
{
    func_legendre_1based( fa, x, pl-1, nl );
}


       /*************************************************
	*						*
	*		NR Utilities			*
	*						*
	*************************************************/

#define NR_END 1
#define FREE_ARG char*

static void 
nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error in tracepoly.c ...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

static double*
vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

#if !USE_SVD
static int *
ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
#endif

static double **
matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

static void 
free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

#if !USE_SVD
static void 
free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
#endif

static void 
free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

#if USE_SVD /*================================================================*/

/*******************************************************************************
 * The following routines are mostly from NR for doing least squares fit
 * using the SVD.
 ******************************************************************************/

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* NR routine svbksb() */
static void 
svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
    int jj,j,i;
    double s,*tmp;

    tmp=vector(1,n);
    for (j=1;j<=n;j++) {
	s=0.0;
	if (w[j]) {
	    for (i=1;i<=m;i++) s += u[i][j]*b[i];
	    s /= w[j];
	}
	tmp[j]=s;
    }
    for (j=1;j<=n;j++) {
	s=0.0;
	for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
	x[j]=s;
    }
    free_vector(tmp,1,n);
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

/* NR routine pythag() */
static double 
pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

/* NR routine svdcmp() */
static void 
svdcmp(double **a, int m, int n, double w[], double **v)
{
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

    rv1=vector(1,n);
    g=scale=anorm=0.0;
    for (i=1;i<=n;i++) {
	l=i+1;
	rv1[i]=scale*g;
	g=s=scale=0.0;
	if (i <= m) {
	    for (k=i;k<=m;k++) scale += fabs(a[k][i]);
	    if (scale) {
		for (k=i;k<=m;k++) {
		    a[k][i] /= scale;
		    s += a[k][i]*a[k][i];
		}
		f=a[i][i];
		g = -SIGN(sqrt(s),f);
		h=f*g-s;
		a[i][i]=f-g;
		for (j=l;j<=n;j++) {
		    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
		    f=s/h;
		    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
		}
		for (k=i;k<=m;k++) a[k][i] *= scale;
	    }
	}
	w[i]=scale *g;
	g=s=scale=0.0;
	if (i <= m && i != n) {
	    for (k=l;k<=n;k++) scale += fabs(a[i][k]);
	    if (scale) {
		for (k=l;k<=n;k++) {
		    a[i][k] /= scale;
		    s += a[i][k]*a[i][k];
		}
		f=a[i][l];
		g = -SIGN(sqrt(s),f);
		h=f*g-s;
		a[i][l]=f-g;
		for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
		for (j=l;j<=m;j++) {
		    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
		    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
		}
		for (k=l;k<=n;k++) a[i][k] *= scale;
	    }
	}
	anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) {
	if (i < n) {
	    if (g) {
		for (j=l;j<=n;j++)
		    v[j][i]=(a[i][j]/a[i][l])/g;
		for (j=l;j<=n;j++) {
		    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
		    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
		}
	    }
	    for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
	}
	v[i][i]=1.0;
	g=rv1[i];
	l=i;
    }
    for (i=IMIN(m,n);i>=1;i--) {
	l=i+1;
	g=w[i];
	for (j=l;j<=n;j++) a[i][j]=0.0;
	if (g) {
	    g=1.0/g;
	    for (j=l;j<=n;j++) {
		for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
		f=(s/a[i][i])*g;
		for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	    }
	    for (j=i;j<=m;j++) a[j][i] *= g;
	} else for (j=i;j<=m;j++) a[j][i]=0.0;
	++a[i][i];
    }
    for (k=n;k>=1;k--) {
	for (its=1;its<=30;its++) {
	    flag=1;
	    for (l=k;l>=1;l--) {
		nm=l-1;
		if ((double)(fabs(rv1[l])+anorm) == anorm) {
		    flag=0;
		    break;
		}
		if ((double)(fabs(w[nm])+anorm) == anorm) break;
	    }
	    if (flag) {
		c=0.0;
		s=1.0;
		for (i=l;i<=k;i++) {
		    f=s*rv1[i];
		    rv1[i]=c*rv1[i];
		    if ((double)(fabs(f)+anorm) == anorm) break;
		    g=w[i];
		    h=pythag(f,g);
		    w[i]=h;
		    h=1.0/h;
		    c=g*h;
		    s = -f*h;
		    for (j=1;j<=m;j++) {
			y=a[j][nm];
			z=a[j][i];
			a[j][nm]=y*c+z*s;
			a[j][i]=z*c-y*s;
		    }
		}
	    }
	    z=w[k];
	    if (l == k) {
		if (z < 0.0) {
		    w[k] = -z;
		    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
		}
		break;
	    }
	    if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
	    x=w[l];
	    nm=k-1;
	    y=w[nm];
	    g=rv1[nm];
	    h=rv1[k];
	    f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	    g=pythag(f,1.0);
	    f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	    c=s=1.0;
	    for (j=l;j<=nm;j++) {
		i=j+1;
		g=rv1[i];
		y=w[i];
		h=s*g;
		g=c*g;
		z=pythag(f,h);
		rv1[j]=z;
		c=f/z;
		s=h/z;
		f=x*c+g*s;
		g = g*c-x*s;
		h=y*s;
		y *= c;
		for (jj=1;jj<=n;jj++) {
		    x=v[jj][j];
		    z=v[jj][i];
		    v[jj][j]=x*c+z*s;
		    v[jj][i]=z*c-x*s;
		}
		z=pythag(f,h);
		w[j]=z;
		if (z) {
		    z=1.0/z;
		    c=f*z;
		    s=h*z;
		}
		f=c*g+s*y;
		x=c*y-s*g;
		for (jj=1;jj<=m;jj++) {
		    y=a[jj][j];
		    z=a[jj][i];
		    a[jj][j]=y*c+z*s;
		    a[jj][i]=z*c-y*s;
		}
	    }
	    rv1[l]=0.0;
	    rv1[k]=f;
	    w[k]=x;
	}
    }
    free_vector(rv1,1,n);
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


#define TOL 1.0e-05

/* NR routine svdfit() 
 * Slightly modified by adding argument of type FuncArgs.
 */
static void
svdfit( FuncArgs* fa, 
        double x[], double y[], double sig[], int ndata, double a[], int ma,
	double **u, double **v, double w[], double *chisq,
        MyFunc* funcs )
{
    int j,i;
    double wmax,tmp,thresh,sum,*b,*afunc;

    b=vector(1,ndata);
    afunc=vector(1,ma);
    for (i=1;i<=ndata;i++) {
	(*funcs)(fa,x[i],afunc,ma);
	tmp=1.0/sig[i];
	for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
	b[i]=y[i]*tmp;
    }
    svdcmp(u,ndata,ma,w,v);
    wmax=0.0;
    for (j=1;j<=ma;j++)
	if (w[j] > wmax) wmax=w[j];
    thresh=TOL*wmax;
    for (j=1;j<=ma;j++)
	if (w[j] < thresh) w[j]=0.0;
    svbksb(u,w,v,ndata,ma,b,a);
    *chisq=0.0;
    for (i=1;i<=ndata;i++) {
	(*funcs)(fa,x[i],afunc,ma);
	for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
	*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
    }
    free_vector(afunc,1,ma);
    free_vector(b,1,ndata);
}
#undef TOL
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

#if 0
/* NR routine svdvar() */
static void
svdvar(double **v, int ma, double w[], double **cvm)
{
    int k,j,i;
    double sum,*wti;

    wti=vector(1,ma);
    for (i=1;i<=ma;i++) {
	wti[i]=0.0;
	if (w[i]) wti[i]=1.0/(w[i]*w[i]);
    }
    for (i=1;i<=ma;i++) {
	for (j=1;j<=i;j++) {
	    for (sum=0.0,k=1;k<=ma;k++) sum += v[i][k]*v[j][k]*wti[k];
	    cvm[j][i]=cvm[i][j]=sum;
	}
    }
    free_vector(wti,1,ma);
}
#endif
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


/*******************************************************************************
 * This is the code for fitting a polynomial to a trace of a read.
 * It is derived from Curtis' polyfit code.
 * Description: Find a legendre polynomial fit to data in x[] and y[].
 *
 * Inputs:   x, y    arrays of data
 *           sig     array of uncertainties of y's (NULL ===> sigs all the same)
 *           n_data  length of both data arrays
 *           na      order of polynomial + 1 (e.g., 4 for a cubic fit)
 * Outputs:  a       Legendre polynomial coefficients; 
 *                     a[n] is coefficient of poly of order x^n
 * Return:   void
 * Comments: If no uncertainties are available pass a NULL pointer for sig.
 *           Uses NR singular value decomposition fitting routine.
 ******************************************************************************/
static void 
legendre_fit( FuncArgs* fa, 
              double x[], double y[], double sig[], int n_data, 
              double a[], int na)
{
    int     i;
    double   **u, **v, *w;
    double  chisq;
    int     free_sig = 0;

    if ( POLY_DBG > 2 ) fprintf(stderr, "  legendre_fit() ...\n");

    /* Allocate memory required by NR (NR's unusual index convention) */
    u = matrix(1, n_data, 1, na);
    v = matrix(1, na, 1, na);
    w = (double *) malloc(na*sizeof(double));

    if ( sig == NULL ) {        /* Create dummy uncertainties */
        sig = (double *) malloc(n_data*sizeof(double));
        free_sig = 1;
        for ( i=0; i < n_data; ++i ) sig[i] = 1.0;
    }

    svdfit( fa, x-1, y-1, sig-1, n_data, a-1, na, u, v, w-1,
            &chisq, func_legendre_1based);
#if POLY_DBG > 1
    if( ReadNum==BadReadNum && BaseIndx==0 ) {
        fprintf(stderr, "    singular values = (");
        for ( i=0; i < na; ++i )
            fprintf(stderr, "%12.3e%s", w[i], i==na-1 ? ")\n" : ",");

    }
#endif

    free_matrix(u, 1, n_data, 1, na);
    free_matrix(v, 1, na, 1, na);
    FREE(w);
    if ( free_sig ) FREE(sig);
}

#else /* USE_SVD =============================================================*/

/*******************************************************************************
 * The following routines are mostly from NR for doing least squares fit
 * NOT using the SVD.
 ******************************************************************************/

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

/* NR routine "gaussj()" */
static void 
gaussj(double **a, int n, double **b, int m)
{
    int *indxc,*indxr,*ipiv;
    int i,icol=0,irow=0,j,k,l,ll;
    double big,dum,pivinv,temp;

    indxc=ivector(1,n);
    indxr=ivector(1,n);
    ipiv=ivector(1,n);
    for (j=1;j<=n;j++) ipiv[j]=0;
    for (i=1;i<=n;i++) {
	big=0.0;
	for (j=1;j<=n;j++)
	    if (ipiv[j] != 1)
		for (k=1;k<=n;k++) {
		    if (ipiv[k] == 0) {
			if (fabs(a[j][k]) >= big) {
			    big=fabs(a[j][k]);
			    irow=j;
			    icol=k;
			}
		    } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
		}
		++(ipiv[icol]);
		if (irow != icol) {
		    for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
		    for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
		    if (ll != icol) {
			dum=a[ll][icol];
			a[ll][icol]=0.0;
			for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
			for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
		    }
    }
    for (l=n;l>=1;l--) {
	if (indxr[l] != indxc[l])
	    for (k=1;k<=n;k++)
		SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n);
    free_ivector(indxc,1,n);
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

/* NR routine "covsrt()" */
static void 
covsrt(double **covar, int ma, int ia[], int mfit)
{
    int i,j,k;
    double swap;

    for (i=mfit+1;i<=ma;i++)
	for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit;
    for (j=ma;j>=1;j--) {
	if (ia[j]) {
	    for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
	    for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
	    k--;
	}
    }
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

/* NR routine "lfit()" */
static void 
lfit(FuncArgs* fa, 
     double x[], double y[], double sig[], int ndat, double a[], int ia[],
     int ma, double **covar, double *chisq, 
     MyFunc* funcs )
{
    int i,j,k,l,m,mfit=0;
    double ym,wt,sum,sig2i,**beta,*afunc;

    beta=matrix(1,ma,1,1);
    afunc=vector(1,ma);
    for (j=1;j<=ma;j++)
	if (ia[j]) mfit++;
    if (mfit == 0) nrerror("lfit: no parameters to be fitted");
    for (j=1;j<=mfit;j++) {
	for (k=1;k<=mfit;k++) covar[j][k]=0.0;
	beta[j][1]=0.0;
    }
    for (i=1;i<=ndat;i++) {
        (*funcs)(fa,x[i],afunc,ma);
	ym=y[i];
	if (mfit < ma) {
	    for (j=1;j<=ma;j++)
		if (!ia[j]) ym -= a[j]*afunc[j];
	}
	sig2i=1.0/SQR(sig[i]);
	for (j=0,l=1;l<=ma;l++) {
	    if (ia[l]) {
		wt=afunc[l]*sig2i;
		for (j++,k=0,m=1;m<=l;m++)
		    if (ia[m]) covar[j][++k] += wt*afunc[m];
		beta[j][1] += ym*wt;
	    }
	}
    }
    for (j=2;j<=mfit;j++)
	for (k=1;k<j;k++)
	    covar[k][j]=covar[j][k];
    gaussj(covar,mfit,beta,1);
    for (j=0,l=1;l<=ma;l++)
	if (ia[l]) a[l]=beta[++j][1];
    *chisq=0.0;
    for (i=1;i<=ndat;i++) {
	(*funcs)(fa,x[i],afunc,ma);
	for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
	*chisq += SQR((y[i]-sum)/sig[i]);
    }
    covsrt(covar,ma,ia,mfit);
    free_vector(afunc,1,ma);
    free_matrix(beta,1,ma,1,1);
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

/*******************************************************************************
 * This is the code for fitting a polynomial to a trace of a read.
 * It is derived from Curtis' polyfit code.
 * Description: Find a legendre polynomial fit to data in x[] and y[].
 *
 * Inputs:   x, y    arrays of data
 *           sig     array of uncertainties of y's (NULL ===> sigs all the same)
 *           n_data  length of both data arrays
 *           na      order of polynomial + 1 (e.g., 4 for a cubic fit)
 * Outputs:  a       Legendre polynomial coefficients; 
 *                     a[n] is coefficient of poly of order x^n
 * Return:   void
 * Comments: If no uncertainties are available pass a NULL pointer for sig.
 *           Does NOT use NR singular value decomposition fitting routine.
 ******************************************************************************/
static void 
legendre_fit(FuncArgs* fa, 
             double x[], double y[], double sig[], int n_data, 
             double a[], int na )
{
    int     i;
    double  **v;
    int     *ia;
    double  chisq;
    int     free_sig = 0;

    if ( POLY_DBG > 2 ) fprintf(stderr, "  legendre_fit() ...\n");

    /* Allocate memory required by NR (NR's unusual index convention) */
    v = matrix(1, na, 1, na);

    if ( sig == NULL ) {        /* Create dummy uncertainties */
        sig = (double *) malloc(n_data*sizeof(double));
        free_sig = 1;
        for ( i=0; i < n_data; ++i ) sig[i] = 1.0;
    }

    ia = (int *) malloc(na*sizeof(int));
    for ( i=0; i < na; ++i ) ia[i] = 1;
    lfit(fa, x-1, y-1, sig-1, n_data, a-1, ia-1, na, v, &chisq, 
         func_legendre_1based);

    FREE(ia);
    free_matrix(v, 1, na, 1, na);
    if ( free_sig ) FREE(sig);
}

#endif /* USE_SVD ============================================================*/


/*******************************************************************************
 * End:   strictly numerical code.
 * Begin: application specific (trace tuner) code.
 ******************************************************************************/


/*******************************************************************************
 * Methods for class TracePoly
 ******************************************************************************/

/* print TracePoly to file */
void 
tracePrint( TracePoly*tp, FILE* fp )
{
    int i;
    fprintf( fp, "coef(%d) = [", tp->num_coef_ );
    for( i=0; i<tp->num_coef_; i++ ) {
        fprintf( fp, " %g", tp->coef_[i] );
    }
    fprintf( fp, "]  " );
    fprintf( fp, "minx=%g max(x,y)=(%g,%g), cutoff(x,y)=(%g,%g)\n",
             tp->x_min_, tp->x_max_, tp->y_max_, tp->x_cutoff_, tp->y_cutoff_ );
}

/* evaluate the TracePoly polynomial */
double 
traceEvalPoly( TracePoly* tp, double x )
{
    int j;
    double pp[NumCoef], y;
    FuncArgs fa;
    if( tp->num_coef_ <= 0 ) {
        return DEFAULT_CONST_POLY;
    }

    fa.poly_type = POLY_FLAG;
    fa.xmin = tp->x_min_;
    fa.xmax = 2.0*tp->x_max_-tp->x_min_;
    func_legendre_0based( &fa, x, pp, tp->num_coef_ );
    
    y = 0.0;
    for( j=0; j<tp->num_coef_; j++ ) { y += tp->coef_[j]*pp[j]; }
 
    y = MYMIN( y, MAX_POLY_LOG_Y );
    y = exp(y);
    y = MYMAX( y, MIN_POLY_Y );
    return y;
}

/* evaluate the TracePoly polynomial with some extra sanity checking */
double 
traceEvalFunc( TracePoly* tp, double x )
{ 
    if( x < tp->x_cutoff_ ) {
        return tp->y_cutoff_;
    } else {
        double y = traceEvalPoly( tp, x );
        if( tp->x_cutoff_ < -1.0 ) { return y; } /* was set to -2.0 elsewhere*/
        else { return MYMIN( y, tp->y_cutoff_ ); }
    }
}

/* This routine assumes: (1) poly is concave up. (2) ymin occurs at xmax.
 * It approximates the largest x value "x_cutoff" (<xmax) such that
 *    poly(x_cutoff) >= y_cutoff
 */
static double
traceCalcCutoff( TracePoly* tp, double y_cutoff, double x_max )
{
    double x, y, delta = x_max/(CUTOFF_NUM_SAMPLES-1);
    
    assert( x_max > 0 );
    assert( x_max < 1000000 );
    for( x=x_max; x>(-delta/2.0); x-=delta ) {
        y = traceEvalPoly( tp, x );
        if( y>y_cutoff ) { return x; }
    }
    return -2.0;        /* < -1.0 (effectively -inf) */
}

/* Function: traceGetCoef
 * This routine computes the polynomial coefficients for the trace.
 * The coefficients are stored in the TracePoly* variable.
 * Everything else is an input.
 * read_num is for debug output only. (any integer will do)
 * The function returns a boolean indicating whether or
 * not it was successful.
 */
static int
traceGetCoef( double xx[MaxNumBases], double yy[MaxNumBases], 
             int num_xy, double ignore_if_less_than,
             TracePoly* tp, int read_num  )
{
    int i, num=0, nc;
    double xxxx[2*MaxNumBases], yyyy[2*MaxNumBases], sigma[2*MaxNumBases];

    for( i=0; i<2*MaxNumBases; i++ ) { sigma[i] = 1.0; }

    for( i=0; i<num_xy; i++ ) {
        double x=xx[i], y, log_y;
        if( x<ignore_if_less_than ) {  continue;  }
        y=yy[i];
        log_y=log(y);        
        xxxx[num] = x;
        yyyy[num] = log_y;
        num++; /* increment AFTER using */
    }
    /* "POINTS_PER_COEF" is somewhat arbitrary.
     * The objective is to avoid over-fitting of the polynomial.
     * (I.e., you want many more points to fit than the 
     * degree of the poly.)
     */
    nc = INT_DBL(num/POINTS_PER_COEF);
    nc = MYMAX(nc,1);
    nc = MYMIN(nc,NumCoef);
    tp->num_coef_ = nc;
    
#if 0
    IFR {
        int i;
        for( i=0; i<num; i++ ) {
            printf( "%4.0f %8.2f\n", xxxx[i], exp(yyyy[i]) );
        }
    }
#endif
    /* Now "flip" the data about tp->x_max_ to force symmetry */
    for( i=0; i<num; i++ ) {
        int i1 = num-1-i, i2=2*num-1-2*i;
        xxxx[i2]   = 2.0*tp->x_max_ - xxxx[i1];
        xxxx[i2-1] = xxxx[i1];
        yyyy[i2]   = yyyy[i1];
        yyyy[i2-1] = yyyy[i1];
    }
    num *= 2;
   
    if( num==0 ) {
        tp->num_coef_ = 1; /* flag ===> poly doesn't exist */
        tp->coef_[0]  = DEFAULT_CONST_POLY;
        tp->x_cutoff_ = -2.0;   /* < -1.0 */
        tp->y_cutoff_ = DEFAULT_CONST_POLY;
        return 1; /* return 0; */  /* not ok */                
    } else {
        FuncArgs fa;
        fa.poly_type = POLY_FLAG;
        fa.xmin = tp->x_min_;
        fa.xmax = 2.0*tp->x_max_-tp->x_min_;
        legendre_fit( &fa, xxxx, yyyy, sigma, num, 
                 tp->coef_, tp->num_coef_ );
    }
    return num;
}


/*******************************************************************************
 * Methods for class ReadInfo
 ******************************************************************************/

/* print ReadInfo to file */
void
readPrint( ReadInfo* read_info, FILE* fp )
{
    static char label[] = { 'A', 'C', 'G', 'T' };
    int i;
    for( i=0; i<4; i++ ) {
        fprintf( fp, "%c: ", label[i] );
        tracePrint( &(read_info->read_polys[i]), fp );
    }
}

/* evaluate all the polys of ReadInfo and compute average */
double
readEvalAveFunc( ReadInfo* read_info, int x )
{
    int base_indx, num=0;
    double y=0.0;
    for( base_indx=0; base_indx<4; base_indx++ ) {
        TracePoly* tp = &(read_info->read_polys[base_indx]);
        if( tp->num_coef_ > 0 ) {
            y += traceEvalFunc( &(read_info->read_polys[base_indx]), x );
            num++;
        }
    }
    if( num ) { return y/num; }
    return DEFAULT_CONST_POLY;
}

    
/* Function readNormFactor
 *   Returns a number by which to multiply a given peak height.
 *   Peak positions are assumed to be in the same units as
 *   those used in the input to "readGetCoefficients()".
 */
double
readNormFactor( ReadInfo* read_info,
                int base_code,      /* 0,1,2,3 ==> A,C,G,T */
                int x )           /* Position of peak */
{
    double y_ref, y, norm;
    const double        /* somewhat arbitrary limits on norm (fudge factors) */
        norm_max = NORMALIZATION_SLACK,
        norm_min = 1.0/norm_max;
    TracePoly* tp;

    if (read_info == NULL) 
        return (1.);

    tp = &(read_info->read_polys[base_code]);

    /* y_ref = 1000.0; */
    /* y_ref = traceEvalPoly( &(rp[4]), x ); */
    y_ref = readEvalAveFunc( read_info, x );

    y = traceEvalFunc( tp, x );

    norm = y_ref/y;

    if( norm < norm_min ) {
#if TELL_ME_IF_I_EXCEED_THE_LIMITS
        if( Filename && ReportNum==0 ) { 
            fprintf(stderr,"Read=%d File=%s",ReadNum,Filename); 
            fprintf( stderr, " norm(%d)=%g < %g\n", x, norm, norm_min );
            Filename=NULL;
            ReportNum++;
        }
#endif
        norm = norm_min;
    }
    if( norm > norm_max ) {
#if TELL_ME_IF_I_EXCEED_THE_LIMITS
        if( Filename && ReportNum==0 ) { 
            fprintf(stderr,"Read=%d File=%s",ReadNum,Filename); 
            fprintf( stderr, " norm(%d)=%g > %g\n", x, norm, norm_max );
            Filename=NULL;
            ReportNum++;
        }
#endif
        norm = norm_max;
    }

    return norm;
}


/* For each trace of the read, compute:
 *    x_min_, x_max_, y_min_, 
 * and set default values for:
 *    x_cutoff_, y_cutoff_, num_coef_, coef_
 */
static void
readCalcXYMaxAndSetDefaults(  
    ReadInfo* read_info,  double xx[4][MaxNumBases], double yy[4][MaxNumBases], 
    int num_xy[4] )
{
    int base_code, i;
    read_info->x_min_ =  DBL_MAX;
    read_info->x_max_ = -DBL_MAX;
    for( base_code=0; base_code<4; base_code++ ) {
        TracePoly* tp = &(read_info->read_polys[base_code]);
        tp->x_min_ = DBL_MAX;
        tp->x_max_ = tp->y_max_ = -DBL_MAX;
        for( i=0; i<num_xy[base_code]; i++ ) {
            tp->x_min_ = MYMIN( xx[base_code][i], tp->x_min_ );
            tp->x_max_ = MYMAX( xx[base_code][i], tp->x_max_ );
            tp->y_max_ = MYMAX( yy[base_code][i], tp->y_max_ );            
        }
        read_info->x_min_ = MYMIN( tp->x_min_, read_info->x_min_ );
        read_info->x_max_ = MYMAX( tp->x_max_, read_info->x_max_ );
    }
    /* Make x_min & x_max the same for all traces */
    for( base_code=0; base_code<4; base_code++ ) {
        TracePoly* tp = &(read_info->read_polys[base_code]);    
        int n;
        tp->x_min_ = read_info->x_min_;
        tp->x_max_ = read_info->x_max_;
        tp->x_cutoff_ = -2.0; /* < -1.0 */
        tp->y_cutoff_ = DEFAULT_CONST_POLY;
        tp->num_coef_ = 0;
        for( n=0; n<NumCoef; n++ ) { tp->coef_[n] = 0.0; }
    }     
}


/* IIR == Infinite Impulse Response */
/* Run an IIR filter backwards & evaluate the result at n=0 */
static double
est_y_at_xmin(  double xx[MaxNumBases], double yy[MaxNumBases], int num )
{
    int n;
    double ave, mean;
    if( num==0 ) { return DEFAULT_CONST_POLY; }
    ave = mean = yy[num-1];
    for( n=num-2; n>=0; n-- ) { 
        ave = (IIR_CONST*ave)+(1.0-IIR_CONST)*yy[n];
        mean += yy[n];
    }
    mean /= num;
    return MYMAX(ave,mean);
}

/* Function: readGetCoefficients
 * This routine computes the polynomials for the read.
 * (One poly per trace)
 * The coefficients are stored in the ReadInfo* variable.
 * Everything else is an input.
 * read_num is for debug output only. (any integer will do)
 * The function returns a boolean indicating whether or
 * not it was successful.
 */
int
readGetCoefficients( ReadInfo* read_info, 
                     double xx_a[4][MaxNumBases], double yy_a[4][MaxNumBases],
                     int num_xy_a[4], int read_num )
{
    int b_indx;

    readCalcXYMaxAndSetDefaults( read_info, xx_a, yy_a, num_xy_a );
    
    for( b_indx=0; b_indx<4; b_indx++ ) {
        TracePoly save, *tp = &(read_info->read_polys[b_indx]);
        int num_xy = num_xy_a[b_indx], ok;
        int one_pass_only = 0;
        double
            *xx = xx_a[b_indx], 
            *yy = yy_a[b_indx],
            y_at_xmin = est_y_at_xmin( xx, yy, num_xy ),
            y_at_xmax,
            y_delta, ignore_if_less_than=0;
#if TELL_ME_IF_I_EXCEED_THE_LIMITS
        BaseIndx = b_indx;      /* set global for debug use */
#endif
        IFR { 
            fprintf( stderr, "before traceGetCoef()\n" );
            tracePrint( tp, stderr );
        }
        ok = traceGetCoef( xx, yy, num_xy, ignore_if_less_than, tp, read_num );
        if( !ok ) { return 0; }
        IFR { 
            fprintf( stderr, "after traceGetCoef()\n" );
            tracePrint( tp, stderr );
        }

        /* I am assuming a quadratic polynomial & testing for concave down */
        if( (tp->num_coef_>=2) && ( tp->coef_[1]<0.0 ) ) {
            tp->coef_[1] = 0.0;
            tp->num_coef_ = 1;
        }

        IFR { 
            fprintf( stderr, "before traceEvalPoly()\n" );
            tracePrint( tp, stderr );
        }

        y_at_xmax = traceEvalPoly( tp, tp->x_max_ );

        IFR { 
            fprintf( stderr, "after traceEvalPoly()\n" );
            tracePrint( tp, stderr );
        }


        y_delta = y_at_xmin - y_at_xmax;

        IFR {
            fprintf( stderr, 
                   "\nA: base=%d, num=%d y(xmin), y(xmax), ydelta= %g %g %g\n", 
                   b_indx, ok, y_at_xmin, y_at_xmax, y_delta );
            tracePrint( tp, stderr );
        }
        
        if( one_pass_only || y_delta <= 0.0 || tp->num_coef_==0 ) { 
            /* Just use the current poly */
            tp->y_cutoff_ = y_at_xmin;
            tp->x_cutoff_ = traceCalcCutoff( tp, tp->y_cutoff_, tp->x_max_ );
            IFR {
                fprintf( stderr, 
                  "\nB1: base=%d, num=%d y(xmin), y(xmax), ydelta= %g %g %g\n", 
                         b_indx, ok, y_at_xmin, y_at_xmax, y_delta );
                tracePrint( tp, stderr );
            }
            continue;
        } else {
            /* Try to fix up the poly */
            double 
                y_cutoff = y_at_xmax + CUTOFF_FACTOR*y_delta,
                x_cutoff = traceCalcCutoff( tp, y_cutoff, tp->x_max_ );
            save = *tp;
            IFR {
                fprintf( stderr, "cutoff(x,y)=(%g,%g)\n", x_cutoff, y_cutoff );
            }
            ok = traceGetCoef( xx, yy, num_xy, x_cutoff, tp, read_num );
            if( !ok || ((tp->num_coef_>=2) && (tp->coef_[1]<=0.0)) ) {
                IFR { fprintf( stderr, "restoring y\n" ); }
                *tp = save;        
            } else {
                tp->y_cutoff_ = y_at_xmin;
                tp->x_cutoff_ = traceCalcCutoff(tp, tp->y_cutoff_, tp->x_max_ );
            }
        }
     
        IFR {
            fprintf( stderr, 
                  "\nB2: base=%d num=%d y(xmin), y(xmax), ydelta= %g %g %g\n", 
                     b_indx, ok, y_at_xmin, y_at_xmax, y_delta );
            tracePrint( tp, stderr );
        }
    }
    return 1;   /* ok */
}

