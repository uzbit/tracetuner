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
 * 2.9 2003/11/06 18:18:49
 */

// Taken from Numerical Recipes in C, 2nd ed.
// NB: Need permission to include this in any commercial product

#include "nr.h"
// #include "nrutil.h"
#include <math.h>

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
	if (j > i) {
	    SWAP(data[j],data[i]);
	    SWAP(data[j+1],data[i+1]);
	}
	m=n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax=2;
    while (n > mmax) {
	istep=mmax << 1;
	theta=isign*(6.28318530717959/mmax);
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0;
	wi=0.0;
	for (m=1;m<mmax;m+=2) {
	    for (i=m;i<=n;i+=istep) {
		j=i+mmax;
		tempr=(float)(wr*data[j]-wi*data[j+1]);
		tempi=(float)(wr*data[j+1]+wi*data[j]);
		data[j]=data[i]-tempr;
		data[j+1]=data[i+1]-tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr=(wtemp=wr)*wpr-wi*wpi+wr;
	    wi=wi*wpr+wtemp*wpi+wi;
	}
	mmax=istep;
    }
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


void realft(float data[], unsigned long n, int isign)
{
    void four1(float data[], unsigned long nn, int isign);
    unsigned long i,i1,i2,i3,i4,np3;
    float c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;

    theta=3.141592653589793/(double) (n>>1);
    if (isign == 1) {
	c2 = -0.5;
	four1(data,n>>1,1);
    } else {
	c2=0.5;
	theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2;i<=(n>>2);i++) {
	i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
	h1r=c1*(data[i1]+data[i3]);
	h1i=c1*(data[i2]-data[i4]);
	h2r = -c2*(data[i2]+data[i4]);
	h2i=c2*(data[i1]-data[i3]);
	data[i1]=(float)(h1r+wr*h2r-wi*h2i);
	data[i2]=(float)(h1i+wr*h2i+wi*h2r);
	data[i3]=(float)(h1r-wr*h2r+wi*h2i);
	data[i4]=(float)(-h1i+wr*h2i+wi*h2r);
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
	data[1] = (h1r=data[1])+data[2];
	data[2] = h1r-data[2];
    } else {
	data[1]=c1*((h1r=data[1])+data[2]);
	data[2]=c1*(h1r-data[2]);
	four1(data,n>>1,-1);
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

void correl(float data1[], float data2[], unsigned long n, float ans[])
{
        void realft(float data[], unsigned long n, int isign);
        void twofft(float data1[], float data2[], float fft1[], float fft2[],
                unsigned long n);
        unsigned long no2,i;
        float dum,*fft;

        fft=vector(1,n<<1);
        twofft(data1,data2,fft,ans,n);
        no2=n>>1;
        for (i=2;i<=n+2;i+=2) {
                ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
                ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
        }
        ans[2]=ans[n+1];
        realft(ans,n,-1);
        free_vector(fft,1,n<<1);
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */

void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[])
{
    int jj,j,i;
    float s,*tmp;

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

void twofft(float data1[], float data2[], float fft1[], float fft2[],
        unsigned long n)
{
        void four1(float data[], unsigned long nn, int isign);
        unsigned long nn3,nn2,jj,j;
        float rep,rem,aip,aim;

        nn3=1+(nn2=2+n+n);
        for (j=1,jj=2;j<=n;j++,jj+=2) {
                fft1[jj-1]=data1[j];
                fft1[jj]=data2[j];
        }
        four1(fft1,n,1);
        fft2[1]=fft1[2];
        fft1[2]=fft2[2]=0.0;
        for (j=3;j<=n+1;j+=2) {
                rep=0.5*(fft1[j]+fft1[nn2-j]);
                rem=0.5*(fft1[j]-fft1[nn2-j]);
                aip=0.5*(fft1[j+1]+fft1[nn3-j]);
                aim=0.5*(fft1[j+1]-fft1[nn3-j]);
                fft1[j]=rep;
                fft1[j+1]=aim;
                fft1[nn2-j]=rep;
                fft1[nn3-j] = -aim;
                fft2[j]=aip;
                fft2[j+1] = -rem;
                fft2[nn2-j]=aip;
                fft2[nn3-j]=rem;
        }
}

/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


void svdcmp(float **a, int m, int n, float w[], float **v)
{
    float pythag(float a, float b);
    int flag, i, its, j, jj, k, l=0, nm=0;
    float anorm, c, f, g, h, s, scale, x, y, z, *rv1;

    rv1=vector(1,n);
    g=scale=anorm=0.0;
    for (i=1;i<=n;i++) {
	l=i+1;
	rv1[i]=scale*g;
	g=s=scale=0.0;
	if (i <= m) {
	    for (k=i;k<=m;k++) scale += (float)fabs(a[k][i]);
	    if (scale) {
		for (k=i;k<=m;k++) {
		    a[k][i] /= scale;
		    s += a[k][i]*a[k][i];
		}
		f=a[i][i];
		g = -(float)SIGN(sqrt(s),f);
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
	    for (k=l;k<=n;k++) scale += (float)fabs(a[i][k]);
	    if (scale) {
		for (k=l;k<=n;k++) {
		    a[i][k] /= scale;
		    s += a[i][k]*a[i][k];
		}
		f=a[i][l];
		g = -(float)SIGN(sqrt(s),f);
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
	anorm=FMAX(anorm,((float)fabs(w[i])+(float)fabs(rv1[i])));
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
	    g=(float)(1.0/g);
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
		if ((float)(fabs(rv1[l])+anorm) == anorm) {
		    flag=0;
		    break;
		}
		if ((float)(fabs(w[nm])+anorm) == anorm) break;
	    }
	    if (flag) {
		c=0.0;
		s=1.0;
		for (i=l;i<=k;i++) {
		    f=s*rv1[i];
		    rv1[i]=c*rv1[i];
		    if ((float)(fabs(f)+anorm) == anorm) break;
		    g=w[i];
		    h=pythag(f,g);
		    w[i]=h;
		    h=(float)(1.0/h);
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
	    f=(float)(((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y));
	    g=pythag(f,1.0);
	    f=(float)(((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x);
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
		    z=(float)(1.0/z);
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

void svdfit(float x[], float y[], float sig[], int ndata, float a[], int ma,
	float **u, float **v, float w[], float *chisq,
	void (*funcs)(float, float [], int))
{
    void svbksb(float **u, float w[], float **v, int m, int n, float b[],
		float x[]);
    void svdcmp(float **a, int m, int n, float w[], float **v);
    int j,i;
    float wmax,tmp,thresh,sum,*b,*afunc;

    b=vector(1,ndata);
    afunc=vector(1,ma);
    for (i=1;i<=ndata;i++) {
	(*funcs)(x[i],afunc,ma);
	tmp=(float)(1.0/sig[i]);
	for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
	b[i]=y[i]*tmp;
    }
    svdcmp(u,ndata,ma,w,v);
    wmax=0.0;
    for (j=1;j<=ma;j++)
	if (w[j] > wmax) wmax=w[j];
    thresh=(float)(TOL*wmax);
    for (j=1;j<=ma;j++)
	if (w[j] < thresh) w[j]=0.0;
    svbksb(u,w,v,ndata,ma,b,a);
    *chisq=0.0;
    for (i=1;i<=ndata;i++) {
	(*funcs)(x[i],afunc,ma);
	for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
	*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
    }
    free_vector(afunc,1,ma);
    free_vector(b,1,ndata);
}
#undef TOL
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


void svdvar(float **v, int ma, float w[], float **cvm)
{
    int k,j,i;
    float sum,*wti;

    wti=vector(1,ma);
    for (i=1;i<=ma;i++) {
	wti[i]=0.0;
	if (w[i]) wti[i]=(float)(1.0/(w[i]*w[i]));
    }
    for (i=1;i<=ma;i++) {
	for (j=1;j<=i;j++) {
	    for (sum=0.0,k=1;k<=ma;k++) sum += v[i][k]*v[j][k]*wti[k];
	    cvm[j][i]=cvm[i][j]=sum;
	}
    }
    free_vector(wti,1,ma);
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


float pythag(float a, float b)
{
	float absa,absb;
	absa=(float)fabs(a);
	absb=(float)fabs(b);
	if (absa > absb) return (float)(absa*sqrt(1.0+SQR(absb/absa)));
	else return (float)(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(float **a, int n, float **b, int m)
{
    int *indxc, *indxr, *ipiv;
    int i, icol=0, irow=0, j, k, l, ll;
    float big, dum, pivinv, temp;

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
			    big=(float)fabs(a[j][k]);
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
		pivinv=(float)(1.0/a[icol][icol]);
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

void covsrt(float **covar, int ma, int ia[], int mfit)
{
    int i,j,k;
    float swap;

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


void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
	int ma, float **covar, float *chisq, void (*funcs)(float, float [], int))
{
    void covsrt(float **covar, int ma, int ia[], int mfit);
    void gaussj(float **a, int n, float **b, int m);
    int i,j,k,l,m,mfit=0;
    float ym,wt,sum,sig2i,**beta,*afunc;

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
	(*funcs)(x[i],afunc,ma);
	ym=y[i];
	if (mfit < ma) {
	    for (j=1;j<=ma;j++)
		if (!ia[j]) ym -= a[j]*afunc[j];
	}
	sig2i=(float)(1.0/SQR(sig[i]));
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
	(*funcs)(x[i],afunc,ma);
	for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
	*chisq += (float)(SQR((y[i]-sum)/sig[i]));
    }
    covsrt(covar,ma,ia,mfit);
    free_vector(afunc,1,ma);
    free_matrix(beta,1,ma,1,1);
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */



void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda)
{
	void covsrt(float **covar, int ma, int ia[], int mfit);
	void gaussj(float **a, int n, float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
		int ia[], int ma, float **alpha, float beta[], float *chisq,
		void (*funcs)(float, float [], float *, float [], int));
	int j,k,l,m;
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=vector(1,ma);
		beta=vector(1,ma);
		da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=(float)0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=(float)(alpha[j][j]*(1.0+(*alamda)));
			oneda[j][1]=beta[j];
		}
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= (float)0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **alpha, float beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int))
{
	int i,j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy,*dyda;

	dyda=vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=(float)(1.0/(sig[i]*sig[i]));
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

int qv_select(unsigned long k, unsigned long n, int arr[])
   /**********************************************************************
    * Returns k-th smallest element of arr[1..n].
    * Comments: Rearranges data in arr[]!  k-th smallest element will be
    *		moved to arr[k]; all smaller elements, to arr[1..k-1]; and
    *		all larger elements, to arr[k+1..n].
    *		Replaced floats with ints (Curtis Gehman, 2001.02.05)
    */

{
	unsigned long i,ir,j,l,mid;
	int a,temp;

	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software <%8+)k[16. */




       /*************************************************
	*						*
	*		NR Utilities			*
	*						*
	*************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error in nr.c ...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
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

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

/* Function: tridiag1 (from Numerical Recipes)
 * Purpose: Solve a tri-diagonal linear system of eqs.
 * Inputs:  a = lower-diag
 *          b = diagonal
 *          c = upper-diag
 *          r = right hand side of linear system
 *          n = size of arrays
 * Output:  u = unknown array
 */
int 
tridiag1( const double a[], const double b[], 
          const double c[], const double r[],
          double u[], int n )
{
    int j;
    double bet,*gam, *gam_temp;
    gam_temp = (double*) malloc( n * sizeof(double) );
    if( !gam_temp ) { return 0; }
    gam = gam_temp-1;
    if (b[1] == 0.0) { free(gam_temp); return 0; }
    u[1]=r[1]/(bet=b[1]);
    for (j=2;j<=n;j++) {
        gam[j]=c[j-1]/bet;
        bet=b[j]-a[j]*gam[j];
        if (bet == 0.0)	{ free(gam_temp); return 0; }
        u[j]=(r[j]-a[j]*u[j-1])/bet;
    }
    for (j=(n-1);j>=1;j--) {
        u[j] -= gam[j+1]*u[j+1];
    }
    free(gam_temp);
    return 1;
}

/* Function: solve_sym_tridiag
 * Purpose: Solve a symmetric tri-diagonal linear system of eqs.
 * Inputs:  diag = diagonal terms
 *          off_diag = off diagonal terms
 *          rhs = right hand side of linear system of eqs.
 *          n = size of arrays
 * Output:  u = unknown array
 */
int
solve_sym_tridiag( const double diag[], const double off_diag[],
                   const double rhs[],
                   double u[], int n )
{
    return tridiag1( off_diag-2, diag-1, off_diag-1, rhs-1, u-1, n );
}


