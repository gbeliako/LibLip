/**************************************************************************

    begin                : April 26 2004
	version				 : 1.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au
	
 	An example of how to use lipax package 

	this program:
	1. randomly generates data
	2. Smoothens the data gived the desired Lipschitz constant
 *                                                                         *
 *  © Gleb Beliakov, 2004												   *
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

#include <cstdlib>
#include <iostream>

#include <fstream>


using namespace std;

#include "../src/slipint.h"

int dim=1;
int	npts=150;  

double sqr_(double a) {return a*a; }
// generates a random number (real) between x and y
double myrand(double x, double y)
{	 
	if(x>y)
		return y + ((x-y)*rand())/(RAND_MAX);
	else if(x<y)
		return x + ((y-x)*rand())/(RAND_MAX);

	return x;
}



// test function, here just a product of sin(2x)sin(2y),...
double fun2(double* dat)
{
//	return 0;
	int j;
	double s=1,s2=1,r;
	for(j=0;j<dim;j++) {
		s*=sin(2*dat[j]);
		r=sin(20*dat[j]);
		s2*=r*r*r;
	}
	return s+0.1*s2;
}

double fun3(double* dat)
{
//	return 0;
	int j;
	double s=0;
	for(j=0;j<dim;j++) {
		s+=dat[j];
	}
	return s;
}

double fun4(double* dat)
{
//	return 0;
	int j;
	double s=0;
	for(j=0;j<dim;j++) {
		s+=dat[j]*dat[j];
	}
	return s;
}


double *x, *XData, *YData, *TData, *SData, noise;
SLipInt LipAx;
int DimTest, NPTSTest;
double  *XDataTest, *YDataTest;
int* Cons;
double *W, *Region;

void ReadTrainFile(char* filename)
{
	int i,j;
	ifstream myFile; 

	myFile.open( filename,ios::in );

	if(!myFile.is_open()) {
		std::cout << "Data file cannot be opened."<<endl;
		return;
	};

	myFile >> dim >> npts ;

	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));
	TData=(double*)malloc(npts*sizeof(double));
	SData=(double*)malloc(npts*sizeof(double));
	 W=(double*)malloc(npts*sizeof(double));
	Cons=(int*)malloc(dim*sizeof(int));
	 Region=(double*)malloc(dim*sizeof(double));

	for(i=0;i<npts;i++) {
		for(j=0;j<dim;j++) {
			myFile >> XData[i*dim+j];
		}
		myFile >> YData[i]; 
	}
	myFile.close();
}
void ReadTestSample(char* filename)
{
	int i,j;
	ifstream myFile; 

	myFile.open( filename,ios::in );

	if(!myFile.is_open()) {
		std::cout << "Data file cannot be opened."<<endl;
		return;
		;
	};

	myFile >> DimTest >> NPTSTest ;
	
	if(dim != DimTest) 
	{
		myFile.close();;
		return;
	}

	XDataTest=(double*)malloc(dim*NPTSTest*sizeof(double));
	YDataTest=(double*)malloc(NPTSTest*sizeof(double));

	for(i=0;i<NPTSTest;i++) {
		for(j=0;j<dim;j++) {
			myFile >> XDataTest[i*dim+j];
		}
		myFile >> YDataTest[i]; 
	}
	myFile.close();
}


void RunTest(double& RMSE, double& MaxErr)
{
	RMSE=MaxErr=0;

	int i;
	double r;

	for(i=0;i<NPTSTest;i++) {
		r=LipAx.ValueLocal(dim,npts,&(XDataTest[i*dim]), XData, YData);

		RMSE += sqr_(r-YDataTest[i]);
		if(MaxErr<fabs(r-YDataTest[i])) MaxErr =fabs(r-YDataTest[i]);
	}

	RMSE/= NPTSTest;

	RMSE= sqrt(RMSE);
}

void Runplot(char* f, double a, double b, int res)
{
	int i,j;
	ofstream myFile; 

	myFile.open( f,ios::out );

	if(!myFile.is_open()) {
		std::cout << "Results file cannot be created."<<endl;
		return ;
	};

	double h=(b-a)/(res+1.);
	double xx,yy;

	myFile<<dim<<" "<<res+1<<endl;

	xx=a;
	for(i=0;i<=res;i++) 
	{
		yy=LipAx.ValueCons(dim,npts,Cons,&xx, XData, YData);
//		yy=LipAx.ValueLocal(dim,npts,&xx, XData, SData);
//		yy=LipAx.Value(dim,npts,&xx, XData, SData);
//		yy=LipAx.ValueLocalCons(dim,npts,Cons,&xx, XData, SData);
		myFile << xx << " "<< yy<<endl;
		xx+=h;
	}

	myFile.close();

}

void Runplot2(char* f, double a, double b, double c, double d, int res)
{
	int i,j;
	ofstream myFile; 

	myFile.open( f,ios::out );

	if(!myFile.is_open()) {
		std::cout << "Results file cannot be created."<<endl;
		return ;
	};

	double hx=(b-a)/(res+1.);
	double hy=(d-c)/(res+1.);
	double xx,yy,zz;

	myFile<<dim<<" "<<(res+1)*(res+1)<<endl;

	xx=a; yy=c;
	double r[2];
	for(i=0;i<=res;i++) 
	{
		yy=c;
		for(j=0;j<=res;j++) 
		{
			r[0]=xx;
			r[1]=yy;
//			zz=LipAx.ValueLocal(dim,npts,r, XData, YData);
//			zz=LipAx.ValueLocalCons(dim,npts,Cons,r, XData, YData);
			zz=LipAx.ValueCons(dim,npts,Cons,r, XData, YData);
			myFile << fabs( xx) << " "<< yy<<" "<<zz<<endl;
			yy+=hy;
		}
		xx+=hx;
	}

	myFile.close();

}


void PrintTriangulation(char* f)
{
	int i,j,iind,m;
	ofstream myFile; 

	myFile.open( f,ios::out );

	if(!myFile.is_open()) {
		std::cout << "Triangulation file cannot be created."<<endl;
		return ;
	};

	for(i=0;i<npts;i++) { // for every datum (row) // parallelisation can be done here
		iind=i*dim;
		myFile << i <<" "<< XData[i]<<" "<<endl;
		for(m=LipAx.pneighbors[i];m<LipAx.pneighbors[i+1];m++) {
			j=LipAx.neighbors[m];
			myFile <<j<<" "<<XData[j] << " ";
		}
		myFile<<LipAx.LocalLipConstants[i]<<endl;
	}
	myFile.close();
}

__cdecl main(int argc, char *argv[])
{
	ofstream myFile; 

//	ReadTrainFile("..//aluminiumdata.txt");
//	ReadTrainFile("..//train2_2.txt");
	ReadTrainFile("..//trainyager.txt");
//	ReadTrainFile("..//train15.txt");
//	ReadTestSample("..//test15.txt");

	for(int i=0;i<dim;i++) Cons[i]=1;

//	npts=200;
	//ComputeLocalLipschitz
//	LipAx.ComputeLocalLipschitzSmooth( dim,  npts, XData,  YData,  TData,  0.5,
//			1,  Cons);
	LipAx.ComputeLipschitz(dim,npts,XData, YData);
//	LipAx.ComputeLocalLipschitz(dim,npts,XData, YData);
//	LipAx.ComputeLocalLipschitzCons(dim,npts,1,Cons,XData, YData);
//	LipAx.ComputeLipschitzSplit(dim,npts,XData, YData, SData,0.5, 0 , Cons, NULL, NULL);
//	LipAx.ComputeLocalLipschitzSmooth(dim,npts,XData, YData, SData,0.5, 0 , Cons, NULL, NULL);
//	LipAx.ComputeLocalLipschitzPenalty(dim,npts,XData, YData, SData,0.05, 0 , Cons, NULL, NULL);
//	LipAx.ComputeLocalLipschitzSmooth(dim,npts,XData, YData, SData,0.06, 0 , NULL, NULL, NULL);
//	LipAx.ComputeLocalLipschitz(dim,npts,XData, YData);

//	PrintTriangulation("..//triang.txt");

;
//	Runplot2("..//plotalum1.txt", -2.3, 0, -0.07,1.13, 100);
//	Runplot2("..//plotalum1.txt", 0, 9, 0,9, 100);
	Runplot2("..//plotalum1.txt", 0, 1, 0,1, 100);
//	Runplot("..//plotalum.txt", -3, 0, 3000);

	x[0]=2.17;

//	LipAx.ValueLocal(dim,npts,x, XData, SData);


double RMSE, MaxErr;

//	RunTest(RMSE,  MaxErr);
//	cout<<"RMSE "<<RMSE<<" MaxErr "<<MaxErr<<endl;

	return 0;
}

/*__cdecl main(int argc, char *argv[])
{	
	int  j,i;
	double a=1.0;


	srand(13);
	double *x, *XData, *YData, *TData, *SData, noise;
// arrays to store the data
	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));
	TData=(double*)malloc(npts*sizeof(double));
	SData=(double*)malloc(npts*sizeof(double));
	double* W=(double*)malloc(npts*sizeof(double));
	int* Cons=(int*)malloc(dim*sizeof(int));
	double* Region=(double*)malloc(dim*sizeof(double));

// generate data randomly
	for(i=0;i<npts;i++) {
		for(j=0;j<dim;j++) {
			x[j]=myrand(3.0,0);
			XData[i*dim + j]=x[j];
		}
		TData[i]=fun4(x) ;
		noise=myrand(0.8, -0.8);
		YData[i]= TData[i] + 0.4*noise; // + error
//		cout<<noise<<endl;

		W[i]=1./dim;
	}
	for(i=0;i<dim;i++) Cons[i]=0;
	Cons[0]=1;
 	Cons[1]=1;

	for(i=0;i<dim;i++) Region[i]=2;


	SLipInt LipAx;
	double LipConst=2.0/1.2; 
//	LipConst=4.0;

	LipAx.ComputeLipschitz(dim,npts,XData, YData);
	cout<<LipAx.MaxLipConst<<endl;

//	int		ComputeLipschitzSplit(int dim, int npts, double* XData, double* YData, double* TData,  double ratio=0.5,
//			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);//	LipAx.ComputeLocalLipschitz(dim,  npts,  XData,  YData);
//	LipAx.ComputeLipschitzCV(dim,npts,XData, YData,SData,1 , Cons, NULL, W);
	LipAx.ComputeLocalLipschitzSmooth(dim,npts,XData, YData, SData,0.7, 0 , NULL, NULL, NULL);
//	LipAx.ComputeLipschitzSplit(dim,npts,XData, YData,SData,0.5, 1 , Cons, NULL, W);
//	LipAx.ComputeLipschitzSplit(dim,npts,XData, YData,SData,0.5);
//	LipAx.ComputeLipschitzCV(dim,npts,XData, YData,SData);
	cout<<LipAx.MaxLipConst<<endl;
	LipConst=LipAx.MaxLipConst;
//	LipConst=2;

//	LipAx.SmoothLipschitzWConsLeftRegion( dim,  npts,  Cons, XData,  YData,  SData,  LipConst,W,Region);

//	for(i=0;i<npts;i++) {
//		cout << YData[i]<<" "<<SData[i]<<endl;
//	}

	for(i=0;i<15;i++) {
		for(j=0;j<dim;j++) 
			x[j]=myrand(3.0,0);

		cout<<fun3(x)<<" "<<LipAx.Value(dim,npts,x,XData,SData,LipConst)<<endl;
	}
//	LipAx.ValueLocal(dim,npts,x,XData,YData);

//	LipAx.SmoothLipschitzInfW( dim,  npts,  XData,  YData,  SData,  LipConst,W);
//	LipAx.SmoothLipschitzSimp( dim,  npts,  XData,  YData,  SData,  LipConst);
//	LipAx.SmoothLipschitz2WCons( dim,  npts,  Cons, XData,  YData,  SData,  LipConst,W);
//	LipAx.SmoothLipschitzInfWCons( dim,  npts,  Cons, XData,  YData,  SData,  LipConst,W);
 
	cout<<LipAx.m_number_constraints<<" "<<LipAx.m_minvalue<<endl;


	double err=0, err1=0, err2=0;
	for(i=0;i<npts;i++) {

		err += fabs(YData[i] - SData[i]);
		err1+= fabs(TData[i] - SData[i]);
		err2+= fabs(TData[i] - YData[i]);
	}
	cout <<"ineq "<< npts*(npts-1)<<endl;
	cout << npts<< " err y-g "<<err/npts<<" err f-g "<<err1/npts<<endl;
	cout << a <<" err y-f "<<err2/npts << endl;



	return 0;
}

*/