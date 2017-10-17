

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


#include "../src/slipint.h"

	SLipClass LipClass;
	double LC;

	double *x, *XData, *YData;
	double  *XDataTr, *YDataTr;
	double  *XDataTe, *YDataTe;
	double   *YDataSmooth;

	int dim,npts,TestNPTS;
int* Cons;
double *W, *Region;

// generates a random number (real) between x and y
double myrand(double x, double y)
{	 
	if(x>y)
		return y + ((x-y)*rand())/(RAND_MAX);
	else if(x<y)
		return x + ((y-x)*rand())/(RAND_MAX);

	return x;
}

void ReadData(int &dim, int npts, int Ntrain, double * XData,double * YData, double * XDataTr, double* YDataTr,double * XDataTe, double* YDataTe)
{
	srand(10);


	int	Ntest  = npts - Ntrain ;
	double P = double(Ntrain)/npts;

	double p;
	int itr=0, ite=0,i,j;
// generate data randomly
	for(i=0;i<npts;i++) {
		if(ite<Ntest && itr<Ntrain) 
			//p=0;
		   p=myrand(1.0,0);
		else if(ite<Ntest) p = 1; else p=0;

		if(p<=P) { // train
			for(j=0;j< dim;j++) 
				XDataTr[itr*dim +j ]=XData[i*dim + j];	
			YDataTr[itr]=YData[i];
			itr++;

		} else { //test
			for(j=0;j< dim;j++) 
				XDataTe[ite*dim +j ]=XData[i*dim + j];	
			YDataTe[ite]=YData[i];
			ite++;
		}
//		cout<<itr<<" "<<ite<<endl;
	}

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
//		yy=LipClass.ValueClass(dim,npts,&xx, XData, YDataSmooth,LC);
//		yy=LipClass.Value(dim,npts,&xx, XData, YDataSmooth,LC);
//		yy=LipClass.ValueCons(dim,npts,Cons,&xx, XData, YDataSmooth,LC);
//		yy=LipClass.ValueConsLeftRegion(dim,npts,Cons,&xx, XData, YDataSmooth,LC,Region);
		yy=LipClass.ValueConsRightRegion(dim,npts,Cons,&xx, XData, YDataSmooth,LC,Region);
//		yy=LipAx.Value(dim,npts,&xx, XData, SData);
//		yy=LipAx.ValueLocalCons(dim,npts,Cons,&xx, XData, SData);
		myFile << xx << " "<< yy<<endl;
		xx+=h;
	}

	myFile.close();

}


void main(int argc, char *argv[])
{	
	double w;
	int k2,K2=100;


	srand(11);

	ifstream myFile; 
//	myFile.open( "cancer683.txt" );
//	myFile.open( "..//classtrain.txt" );
//	myFile.open( "..//concentric.dat" );
//	myFile.open( "..//cancer683.txt" );
	myFile.open( "..//wdbcclass.txt" );

	int  j,i;
	myFile >> dim >> npts;
//	npts=20000;
	int Cnt1=0,Cnt2=0;

// arrays to store the data
	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));
	YDataSmooth=(double*)malloc(npts*sizeof(double));

	for(i=0;i<npts;i++) {
		for(j=0;j< dim;j++) {
			myFile >> XData[i*dim + j];
//			XData[i*dim + j]=-XData[i*dim + j];
		}
		myFile >> YData[i];
		if(YData[i]==1) Cnt1++; else Cnt2++;
//		if(YData[i]==4) YData[i]=1; else YData[i]=-1;
	}
	myFile.close();


// randomly split into train and test
	int Ntrain, Ntest;
	Ntrain = npts;
	Ntest  = npts - Ntrain ;

	XDataTr=(double*)malloc(dim*Ntrain*sizeof(double));
	YDataTr=(double*)malloc(Ntrain*sizeof(double));

	XDataTe=(double*)malloc(dim*Ntest*sizeof(double));
	YDataTe=(double*)malloc(Ntest*sizeof(double));

//	ReadData(dim, npts, Ntrain, XData,YData,  XDataTr, YDataTr, XDataTe, YDataTe);


//	free(XData); free(YData); free(XDataTe); free(YDataTe);

// if necessary, compute Lipschitz constant (slow)
	double* Lip=(double*)malloc(npts*sizeof(double));
	for(i=0;i<Ntrain;i++) {
			Lip[i]=1.152;
	}

//	double beta=0.5225/Ntrain;
	double beta=0.4225/10;

	 W=(double*)malloc(npts*sizeof(double));
	Cons=(int*)malloc(dim*sizeof(int));
	 Region=(double*)malloc(dim*sizeof(double));

	double w1,w2, err, err2; // compute the error of approximation
	err2=err=0;
	int counterr=0, ii;

	for( i=0;i<dim;i++) Cons[i]=0;

	Region[0]=-3;


	LipClass.SmoothingParam=0.00;
	LC=1.99;

	LipClass.ComputeScaling(dim,npts,XData,YData);

	for(i=0;i<npts;i++) {
		for(j=0;j< dim;j++) {
			XData[i*dim + j] = XData[i*dim + j]*LipClass.Scaling[j];
		}
	}

//	PrepareLipschitzSplit(0.9);
	TestNPTS=569;
	LipClass.ComputeLipschitz(dim,TestNPTS,XData,YData);
	LC=LipClass.MaxLipConst;
	cout<<LC <<endl;

	int r;
	for(i=0;i<dim;i++){
		Cons[i]=1;
		r=LipClass.VerifyMonotonicity(dim,npts,Cons,XData, YData,LC);
		if(r<1) {
			Cons[i]=-1;
		}
		r=LipClass.VerifyMonotonicity(dim,npts,Cons,XData, YData,LC);
		cout<<r<<endl;
	}
	Cons[dim-1]=1;
	cout<<LipClass.VerifyMonotonicity(dim,npts,Cons,XData, YData, LC)<<endl;
	Cons[dim-1]=-1;
	cout<<LipClass.VerifyMonotonicity(dim,npts,Cons,XData, YData,LC)<<endl;


	LipClass.ComputeLipschitz(dim,TestNPTS,XData,YData);
	LC=LipClass.MaxLipConst;
	cout<<LC <<endl;
	LC /=1;
//	LC= 10;

	double confusion1a,confusion2a,confusion3a,confusion4a;
	double confusion1,confusion2,confusion3,confusion4;
	confusion1=confusion2=confusion3=confusion4=0;
	confusion1a=confusion2a=confusion3a=confusion4a=0;
	int totclass1=0, totclass2=0;

for(k2=0;k2<1; k2++) {
	LipClass.NPTS=npts;
	LipClass.PrepareLipschitzSplit(0.90);

	Cnt1=Cnt2=0;
	for(i=0;i<LipClass.Indexsize;i++) {
			if(YData[LipClass.Index[i]]==1) Cnt1++; else Cnt2++;
	}

	for(i=0;i<LipClass.Indexsize;i++) {
		if(YData[LipClass.Index[i]]==1)
			W[LipClass.Index[i]]=1.0/Cnt1; else W[LipClass.Index[i]]=1.0/Cnt2;
	}


//	LipClass.SmoothLipschitz2Classinternal(dim,LipClass.Indexsize, XData,  YData,  YDataSmooth, 0,1,0, &LC,  W,NULL,0,NULL, LipClass.Index);



//	LipClass.SmoothLipschitzWClass(dim,TestNPTS,XData,YData,YDataSmooth,&LC,W);

//	LipClass.SmoothLipschitzConsClass(dim,npts,Cons,XData,YData,YDataSmooth,&LC);
//	LipClass.SmoothLipschitzConsLeftRegionClass(dim,npts,Cons,XData,YData,YDataSmooth,&LC, Region);
//	LipClass.SmoothLipschitzConsRightRegionClass(dim,npts,Cons,XData,YData,YDataSmooth,&LC, Region);
//	SmoothClassifier( dim,  Ntrain,  beta,  XDataTr,  YDataTr,  Lip);	
//	SmoothClassifierLimited( dim,  Ntrain, 10, beta,  XDataTr,  YDataTr,  Lip);	

//	Runplot("..//plotclass.txt",-10,8,1000);


	Ntest=LipClass.IndexsizeComp;
//	Ntest=npts-TestNPTS;

	for(i=0;i<Ntest;i++) {
		ii=LipClass.IndexComp[i];
		for(j=0;j<dim;j++) x[j]=XData[ii *dim + j ];

//		w=LipClass.Value(dim,LipClass.Indexsize,x, XData, YDataSmooth,LC, LipClass.Index);
		w=LipClass.Value(dim,LipClass.Indexsize,x, XData, YData,LC, LipClass.Index);
		if(w>=0) w=1; else w=-1;
//		w=LipClass.ValueClass(dim,TestNPTS,x, XData, YData,200);

		if(YData[ii]==1) {
			if(w==1) confusion1++; else confusion2++;
			totclass1++;
		} else { // -1
			if(w==1) confusion3++; else confusion4++;
			totclass2++;
		}
		
	}
	cout << (double)confusion1 /totclass1*100 << " "<<(double)confusion2/totclass1*100<<endl;
	cout << (double)confusion3/totclass2*100 << " "<<(double)confusion4/totclass2*100<<endl;
} // k2
	cout << (double)confusion1 /totclass1*100 << " "<<(double)confusion2/totclass1*100<<endl;
	cout << (double)confusion3/totclass2*100 << " "<<(double)confusion4/totclass2*100<<endl;
	//	cout <<"Lip "<<LipConst<<endl;

}

/* abalone
void main(int argc, char *argv[])
{	
	double w;
	int k2,K2=1000;
	double *x, *XData, *YData;
	double  *XDataTr, *YDataTr;
	double  *XDataTe, *YDataTe;

	ifstream myFile; 
//	myFile.open( "cancer683.txt" );
	myFile.open( "abalone.txt" );

	int dim,npts, j,i;
	myFile >> dim >> npts;


// arrays to store the data
	x=(double*)malloc((dim+1)*sizeof(double));
	XData=(double*)malloc(dim*npts*sizeof(double));
	YData=(double*)malloc(npts*sizeof(double));
	

	for(i=0;i<npts;i++) {
		for(j=0;j< dim;j++) {
			myFile >> XData[i*dim + j];

			if(XData[i*dim + j]==1)XData[i*dim + j]=2;
			else if(XData[i*dim + j]==1)XData[i*dim + j]=1;
		}
		myFile >> YData[i];
		
	}
	myFile.close();


// randomly split into train and test
	int Ntrain, Ntest;
	Ntrain = 3133;
	Ntest  = npts - Ntrain ;
	double P = double(Ntrain)/npts;

	XDataTr=(double*)malloc(dim*Ntrain*sizeof(double));
	YDataTr=(double*)malloc(Ntrain*sizeof(double));

	XDataTe=(double*)malloc(dim*Ntest*sizeof(double));
	YDataTe=(double*)malloc(Ntest*sizeof(double));

	srand(11);

	double p;
	int itr=0, ite=0;
// generate data randomly
	for(i=0;i<npts;i++) {
		if(ite<Ntest && itr<Ntrain) 
			p=0;
		  // p=myrand(1.0,0);
		else if(ite<Ntest) p = 1; else p=0;

		if(p<=P) { // train
			for(j=0;j< dim;j++) 
				XDataTr[itr*dim +j ]=XData[i*dim + j];	
			YDataTr[itr]=YData[i];
			itr++;

		} else { //test
			for(j=0;j< dim;j++) 
				XDataTe[ite*dim +j ]=XData[i*dim + j];	
			YDataTe[ite]=YData[i];
			ite++;
		}
		cout<<itr<<" "<<ite<<endl;
	}

/*	for(i=0;i<Ntest;i++) {
		if(YDataTe[i]<9) YDataTe[i]=1;
		else if(YDataTe[i]<10.99) YDataTe[i]=2;
		else YDataTe[i]=3;
	}

/*	for(i=0;i<Ntrain;i++) {
		if(YDataTr[i]<9) YDataTr[i]=1;
		else if(YDataTr[i]<10.99) YDataTr[i]=2;
		else YDataTr[i]=3;
	}

	// interpolant object
	STCInterpolant LipInt;
// set Lipschitz const
	LipInt.SetData(dim,Ntrain,XDataTr,YDataTr,1);  
	//   assumes all the data are distinct. 
	//	If this needs to be tested, use  LipInt.GetData(dim,npts,&(XData[0][0]),YData.begin(), 1); command, with the last
	//  parameter =1. This may be slow. 

// if necessary, compute Lipschitz constant (slow)
	double LipConst=1200;
//	LipConst=LipInt.DetermineLipschitz();

	LipInt.SetConstants(LipConst);

	LipInt.ConstructExplicit(); 

	double w1,w2, err, err2; // compute the error of approximation
	err2=err=0;
	int counterr=0;

	for(i=0;i<Ntest;i++) {
		for(j=0;j<dim;j++) x[j]=XDataTe[i *dim +j ];
		w=LipInt.ValueExplicit(dim,x);
		w2=w;
/*		if(w2<9) w=1;
		else if(w2<10.99) w=2;
		else w=3;
/*
		if(w2<1.5) w=1;
		else if(w2<2.5) w=2;
		else w=3;


	//	cout << w << " " << YDataTe[i]<<endl;
		w1=fabs(w-YDataTe[i]);
		if(w1>=0.5) {counterr++;
		cout<<x[0]<<" "<<w2<<" " <<w<<" "<<YDataTe[i]<<endl;
		}

		err2+=w1*w1;
	}

	err2=sqrt(err2/K2);  // average error RMSE
	cout << "average error RMSE "<<err2 <<endl;
	err=double(counterr)/Ntest;
	cout << "errors: "<<counterr<< " "<<err<<":"<<1-err<<" "<<Ntest<<endl;
	cout <<"Lip "<<LipConst<<endl;

}

*/