#include <iostream>
#include <cmath>
#include <math.h>
#include "molecule.h"

//setting the dimensions of all vectors
const int n=3;

//declaring necessary vectors and variables
//vector< vector<double> > a{{156.154,0,52.8556},{0,199.371,0},{52.8556,0,54.4595}};
Matrix rot(n,n);

double theta=0;
double t=0,c=0,s=0,tau=0;
int sweep_count=0;

//function to reset the rotation matrix Ppq for every p,q
void initialize_rot()
{
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<n;k++)
		{
			if(i==k)
				rot(i,k)=1;
			else
				rot(i,k)=0;
		}
	}
}

//signum function			
double signum(double num)
{
	if(num<0.0)
		return -1.0;
  else 
		return 1.0;
}

//function to update the sum of the off-diagonal elements for convergence
double update_s(Matrix vec)
{
	double sum=0;
	for(int p=0;p<n;p++)
	{
		for(int q=0;q<n;q++)
		{
			if(q!=p)
				sum+=abs(vec(p,q));
		}
	}
	return sum;
	//cout<<sum<<'\n';
}

//function to do the jacobi transformation
void jacobi(Matrix &veca, Matrix &vecb)
{
	double s0=update_s(veca);
	while(update_s(veca) >pow(10,-20))
	{
		sweep_count++;
		for(int p=0;p<n;p++)
		{
			for(int q=p+1;q<n;q++)
			{
				initialize_rot();
			//all variables created here
				theta=(veca(q,q)-veca(p,p))/(2.0*veca(p,q));
				t = signum(theta)/ (abs(theta) + pow(pow(theta,2.0) + 1.0,0.5));
				c=1.0/sqrt(pow(t,2.0)+1.0);
				s=c*t;
				tau=s/(1.0+c);
			//updating diagonal elements of the vector to be diagonalized
				veca(p,p)=veca(p,p)-t*veca(p,q);
				veca(q,q)=veca(q,q)+t*veca(p,q);
			//updating non-diagonal elements of the vector
				for(int r=0;r<n;r++)
				{
					if(r!=p && r!=q)
					{
						double temp=veca(r,p);
						double temp1=veca(r,q);
						veca(r,p)=temp-(s*(temp1+(tau*temp)));
						veca(p,r)=veca(r,p);
						veca(r,q)=temp1+(s*(temp-(tau*temp1)));
						veca(q,r)=veca(r,q);
					}
				}
				veca(p,q)=0.0;
				veca(q,p)=0.0;
			//creating the rotation matrix for p,q
				rot(p,p)=c;
				rot(q,q)=c;
				rot(p,q)=s;
				rot(q,p)=-s;
			//multiplying to get the eigenvector matrix
				vecb=vecb*rot;
				//print_matrix(vecb);
				//cout<<"theta: "<<theta<<"\tt: "<<t<<"\tc: "<<c<<"\ts: "<<s<<'\t'<<'\n';
			  cout<<'\n'<<veca<<'\n';
				//cout<<'\n';
				//print_matrix(rot);
				//cout<<'\n';
			}
		}
	}
	cout<<"\nJacobi done.\nNumber of sweeps: "<<sweep_count<<'\n';
}

