/**************************************************************************

    begin                : April 19 2005
	version				 : 2.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au
	
 	An example of how to use liblip package with class and procedural interface

	this program:
	1. randomly generates data
	2. Builds the interpolant
	3. computes the value of the interpolant and compates it with the test 
	   data (model function)
	4. reports preprocessing and evaluation time and the accuracy of 
	   approximation
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

#include "../include/liblip.h"
#include "../include/liblipc.h"

int dim=3;
int	npts=100;

// test function, here just a product of sin(2x)sin(2y),...
double fun2(double* dat)
{	int j;
	double s=1;
	for(j=0;j<dim;j++) 	s*=sin(2*dat[j]);
	return s;
}

// generates a random number (real) between x and y
double myrand(double x, double y)
{	 
	if(x>y)	     return y + ((x-y)*rand())/(RAND_MAX);
	else if(x<y) return x + ((y-x)*rand())/(RAND_MAX);
	return x;
}

int  main(int argc, char *argv[])
{	
	int  j,i;
	double w;
	int k2,K2=100;
	double w1, err, err2; // compute the error of approximation

	double *x, *XData, *YData;
// arrays to store the data
	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));
	double *TData=(double*)malloc(npts*sizeof(double));

// generate data randomly
	for(i=0;i<npts;i++) {
		for(j=0;j<dim;j++) {
			x[j]=myrand(3.0,0);
			XData[i*dim + j]=x[j];
		}
		YData[i]=fun2(x);
	}

// simple Lipschitz interpolant section

	printf("Using SLipInt class ...\n");

	SLipInt SLi;

	SLi.ComputeLocalLipschitz(dim,npts,XData,YData);
//	printf("Local Lipschitz penalty %f\n",SLi.m_minvalue);

	err2=err=0;
	for(k2=0;k2<K2;k2++) {
		for(j=0;j<dim;j++) 	x[j]=myrand(3.0,0);
		w=SLi.Value(dim,npts,x,XData, YData);
		w1=fun2(x);				// the true function 
		w=fabs(w-w1);				// compute the error 
		if(err<w) err=w;
		err2+=w*w;
	}
	err2=sqrt(err2/K2);  // average error RMSE
	printf("Interpolation max error %f\n",err);
	printf("average error %f\n",err2);

	printf("Using procedural interface ...\n");
	double ratio=0.5;
	int type=0;
   LipIntComputeLipschitzSplit(&dim,&npts,  XData,YData,TData,&ratio,&type,0,0,0);

   LipIntComputeLocalLipschitz(&dim,&npts,  XData,YData);

   w=LipIntValueLocal(&dim,&npts, x, XData, TData);

   double LC=10;
	for(k2=0;k2<10;k2++) {
		for(j=0;j<dim;j++) 	x[j]=myrand(3.0,0);
		w=LipIntValue(&dim,&npts,x,XData, YData, &LC,0);
		w=LipIntInfValue(&dim,&npts,x,XData, YData, &LC,0);
	}

	printf("Using STC methods ...\n");

	STCSetLipschitz(&LC);
	STCBuildLipInterpolant(&dim, &npts, XData, YData);

	for(k2=0;k2<K2;k2++) {
		for(j=0;j<dim;j++) 	x[j]=myrand(3.0,0);
		w=STCValue(x);
	}

	STCFreeMemory();


	printf("Using STCInterpolant class ...\n");
	// interpolant object
	STCInterpolant LipInt;
// set Lipschitz const
	LipInt.SetData(dim,npts,XData,YData,0);  
	//   assumes all the data are distinct. 
	//	If this needs to be tested, use  LipInt.GetData(dim,npts,&(XData[0][0]),YData.begin(), 1); command, with the last
	//  parameter =1. This may be slow. 

// if necessary, compute Lipschitz constant (slow)
//	LipConst=LipInt.DetermineLipschitz();
	LipInt.SetConstants(5);
	LipInt.Construct(); 
	err2=err=0;

	for(k2=0;k2<K2;k2++) {
		for(j=0;j<dim;j++) 	x[j]=myrand(3.0,0); // randomly choose a test point
		w=LipInt.Value(dim,x);	// evaluate the interpolant
		w1=fun2(x);				// the true function 
		w=fabs(w-w1);				// compute the error 
		if(err<w) err=w;
		err2+=w*w;
	}


	err2=sqrt(err2/K2);  // average error RMSE
	printf("max interpolation error %f\n",err);
	printf("average error %f\n",err2);

	LipInt.FreeMemory();

return 0;
}

