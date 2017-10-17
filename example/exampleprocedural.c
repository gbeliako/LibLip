/**************************************************************************

    begin                : April 19 2004
    version		 : 1.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au
	
 	An example of how to use lint package with procedural interface

	this program:
	1. randomly generates data
	2. Builds the interpolant
	3. computes the value of the interpolant and compates it with the test 
	   data (model function)
	4. reports preprocessing and evaluation time and the accuracy of 
	   approximation
 *                                                                         *
 *   Gleb Beliakov, 2004												   *
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

int dim=3;
int	npts=1500;

// test function, here just a product of sin(2x)sin(2y),...
double fun2(double* dat)
{
//	return 0;
	int j;
	double s=1;
	for(j=0;j<dim;j++)
		s*=sin(2*dat[j]);
	return s;
}

// generates a random number (real) between x and y
double myrand(double x, double y)
{	 
	if(x>y)
		return y + ((x-y)*rand())/(RAND_MAX);
	else if(x<y)
		return x + ((y-x)*rand())/(RAND_MAX);

	return x;
}

int main(int argc, char *argv[])
{
	int  j,i;
	double w;
	int k2,K2=100;

	double *x, *XData, *YData;
// arrays to store the data
	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));

// generate data randomly
	for(i=0;i<npts;i++) {
		for(j=0;j<dim;j++) {
			x[j]=myrand(3.0,0);
			XData[i*dim + j]=x[j];
		}
		YData[i]=fun2(x);
	}

	double LipConst=10;

	STCSetLipschitz(&LipConst);
	STCBuildLipInterpolant(&dim, &npts, XData, YData);

	for(k2=0;k2<K2;k2++) {
		for(j=0;j<dim;j++) 	x[j]=myrand(3.0,0);
		w=STCValue(x);
	}
     printf("completed.\n"); 
	return 0;
}

