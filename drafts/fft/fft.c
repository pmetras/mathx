#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double *data, int n, const int isign) {
/* Replaces data[0..2*n-1] by its discrete Fourier transform, if isign is input as 1; or replaces
data[0..2*n-1] by n times its inverse discrete Fourier transform, if isign is input as -1. data
is a complex array of length n stored as a real array of length 2*n. n must be an integer power
of 2. */
	int nn,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) { printf("n must be power of 2 in four1"); exit(1); }
	nn = n << 1;
	j = 1;
	for (i=1;i<nn;i+=2) {
		if (j > i) {
			SWAP(data[j-1],data[i-1]);
			SWAP(data[j],data[i]);
		}
		m=n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}


void realft(double *data, int n, const int isign) {
/* Calculates the Fourier transform of a set of n real-valued data points. Replaces these data
(which are stored in array data[0..n-1]) by the positive frequency half of their complex Fourier
transform. The real-valued ﬁrst and last components of the complex transform are returned
as elements data[0] and data[1], respectively. n must be a power of 2. This routine also
calculates the inverse transform of a complex data array if it is the transform of real data.
(Result in this case must be multiplied by 2/n.) */
	int i,i1,i2,i3,i4;
	double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp;
	double theta=3.141592653589793238/((double) (n>>1));
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	for (i=1;i<(n>>2);i++) {
		i2=1+(i1=i+i);
		i4=1+(i3=n-i1);
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r= -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[0] = (h1r=data[0])+data[1];
		data[1] = h1r-data[1];
	} else {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(data,n,-1);
	}
}


int main(int argc, char** argv) {
	int size = 8;
	double test[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
	double *data;
	data = (double *) malloc(sizeof(test));
	int i;
	for (i = 0; i < size; i++) {
		data[i] = test[i];
	}
	for (i = 0; i < size; i++) {
		printf("orig %d --> %f\n", i, data[i]);
	}
	printf("\n");
	//four1(data, size/2, 1);
	realft(data, size, 1);
	for (i = 0; i < size; i++) {
		printf("fft %d --> %f\n", i, data[i]);
	}
}

