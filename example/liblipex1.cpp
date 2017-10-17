/**************************************************************************

    begin                : April 19 2005
	version				 : 2.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au
	
 	An example of how to use liblip package with class  interface,
	using additional upper and lower bounds

	this program:
	1. randomly generates data
	2. derives a class from SLipInt, which uses upper/lower bounds
	3. computes the value of the interpolant subject to monotonicity

 *                                                                         *
 *  © Gleb Beliakov, 2005												   *
 *                                                                         *
 * This program is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU General Public License as published by the   *
 * Free Software Foundation; either version 2 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * This program is distributed in the hope that it will be useful, but     *
 * WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * General Public License for more details.                                *
 *                                                                         *
 * You should have received a copy of the GNU General Public License       *
 * along with this program; if not, write to the Free Software Foundation, *
 * Inc., 59 Temple Place Suite 330, Boston, MA 02111-1307 USA.             *
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#define MATH_UTILS_H
#include "../src/liblip.h"
#include "../src/liblipc.h"

int dim=4;
int	npts=100;

// test function, here just a product x1*x2*.....
double fun2(double* dat)
{	int j;
	double s=1;
	for(j=0;j<dim;j++) 	s*=dat[j];  
	return s;
}

// generates a random number (real) between x and y
double myrand(double x, double y)
{	 
	if(x>y)	     return y + ((x-y)*rand())/(RAND_MAX);
	else if(x<y) return x + ((y-x)*rand())/(RAND_MAX);
	return x;
}

/* This class illustrates using upper and lower bounds and monotone apprxomation
  we require the approximation f(x) to satisfy Lip * x[0]/20 <=f(x)<= max(x[0],x[1],...)
*/
class SLipIntDerived:public SLipInt
{
public:
	virtual double  ExtraUpperBound(int dim, double* x, double * param) {
		int i; double t=0;
		for(i=0;i<dim;i++) if(x[i]>t) t=x[i]; // compute max(x)
		return t;
	};
	virtual double  ExtraLowerBound(int dim, double* x, double * param) {
		// here we use the Lipschitz constant passed in param
		return x[0]* (*param)/20.0;
	};

};

void main(int argc, char *argv[])
{	
	int  j,i;
	double w;
	int k2,K2=100;
	double w1, w2, err, err2; // compute the error of approximation

	double *x, *XData, *YData;
// arrays to store the data
	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));
	double *TData=(double*)malloc(npts*sizeof(double));
	int *Cons=(int*)malloc(dim*sizeof(int));

// generate data randomly
	for(i=0;i<npts;i++) {
		for(j=0;j<dim;j++) {
			x[j]=myrand(1.5,0);
			XData[i*dim + j]=x[j];
		}
		YData[i]=fun2(x);
	}


	for(j=0;j<dim;j++)  Cons[j]=0; // require monotone increasing f

	Cons[dim-1]=-1;
// simple Lipschitz interpolant section

	printf("Using SLipInt class with bounds...\n");

	double LC=2.0;
	SLipIntDerived SLi;

// if we change the next line to =0, bounds will not be used, and 
// and the approximation will probably fail to be between the bounds.
// the errors will appear in the printout.	
	SLi.UseOtherBounds=1;

	SLi.SmoothLipschitzCons(dim,npts,Cons,XData,YData,TData,LC);

// checking bounds
	for(k2=0;k2<K2;k2++) {
		for(j=0;j<dim;j++) 	x[j]=myrand(1.5,0);
		w=SLi.ValueCons(dim,npts,Cons,x,XData, TData, LC);
		w1=SLi.ExtraLowerBound(dim,x,&LC);
		w2=SLi.ExtraUpperBound(dim,x,&LC);
		// w is supposed to be in [w1,w2] 
		if(w<w1 || w>w2)
			printf("bounds not satisfied %f <= %f <=%f \n", w1,w,w2);
	}

	printf("test completed\n");
}

