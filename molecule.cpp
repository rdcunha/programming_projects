#include "molecule.h"
#include "masses.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>

void Molecule::read_in()
{
  ifstream input(iname);
  //if(!input) error("Can't open input file ",iname);
  //read number of atoms into no_atoms
  input>>no_atoms;
  cout<<no_atoms<<'\n';
  z_vals.resize(no_atoms);
  geom.resize(3,no_atoms);
  //read zvals and coords into single multidimensional array
  while(input)
  {
    for(int j=0; j<no_atoms;j++)
    {
      input>>z_vals[j];
      for(int i=0;i<3;i++)
      {
        input>>geom(i,j);
      }
    }
  }
}

void Molecule::print_geom()
{
  for(int i=0; i<no_atoms;i++)
  {
    printf("%d\t%6.6f\t%6.6f\t%6.6f\t\n",z_vals[i],geom(0,i),geom(1,i),geom(2,i));
  }
}

void Molecule::translate(double x, double y, double z)
{
  for(int i=0; i<no_atoms;i++)
  {
    geom(0,i)+=x;
    geom(1,i)+=y;
    geom(2,i)+=z;
  }
}

double Molecule::bond(int a1,int a2)
{
  double length=0;
  for(int j=0;j<3;j++)
  {
    length+=pow(geom(j,a1)-geom(j,a2),2.0);
  }
  length=sqrt(length);
  return length;
}

double Molecule::angle(int a1, int a2, int a3)
{
  double ang=0;
  for(int j=0;j<3;j++)
  {
    ang+=((geom(j,a1)-geom(j,a2))/bond(a2,a1))*((geom(j,a3)-geom(j,a2))/bond(a3,a2));
  }
  ang=acos(ang)*(180.0/acos(-1.0));
  return ang;
}

double Molecule::oop(int a1, int a2, int a3, int a4)
{
  double oops=0;
  vector <vector<double> > coords(3,vector<double>(3));
  for(int j=0;j<3;j++)
  {
    coords[j][0]=(geom(j,a2)-geom(j,a3))/bond(a2,a3);
    coords[j][1]=(geom(j,a4)-geom(j,a3))/bond(a4,a3);
    coords[j][2]=(geom(j,a1)-geom(j,a3))/bond(a1,a3);
  }
  for(int i=0;i<3;i++)
  {
    oops+=coords[i][2]*(coords[(i+1)%3][0]*coords[(i+2)%3][1]-coords[(i+2)%3][0]*coords[(i+1)%3][1]);
    //std::cout<<'\n'<<i<<" 2.("<<(i+1)%3<<" 0 x "<<(i+2)%3<<" 1 - "<<(i+2)%3<<" 0 x "<<(i+1)%3<<" 1 )\n";
  }
  oops*=1.0/(sin(angle(a2,a3,a4)*acos(-1.0)/180.0));
  if(oops < -1.0)
    oops=asin(-1.0)*180.0/acos(-1.0);
  else if(oops > 1.0)
    oops=asin(1.0)*180.0/acos(-1.0);
  else
    oops=asin(oops)*180.0/acos(-1.0);
  return oops;
}

double Molecule::torsion(int a1, int a2, int a3, int a4)
{
  double tors=0;
  vector <vector<double> > coords(3,vector<double>(4));
  for(int j=0;j<3;j++)
  {
    coords[j][0]=(geom(j,a2)-geom(j,a1))/bond(a2,a1);
    coords[j][1]=(geom(j,a3)-geom(j,a2))/bond(a3,a2);
    coords[j][2]=(geom(j,a4)-geom(j,a3))/bond(a4,a3);
  }
  for(int i=0;i<3;i++)
  {
    tors+=(coords[(i+1)%3][0]*coords[(i+2)%3][1]-coords[(i+2)%3][0]*coords[(i+1)%3][1])*(coords[(i+1)%3][1]*coords[(i+2)%3][2]-coords[(i+2)%3][1]*coords[(i+1)%3][2]);
  }
  tors*=1.0/sin(angle(a1,a2,a3)*acos(-1.0)/180.0);
  tors*=1.0/sin(angle(a2,a3,a4)*acos(-1.0)/180.0);
  if(tors < -1.0)
    tors=acos(-1.0)*180.0/acos(-1.0);
  else if(tors > 1.0)
    tors=acos(1.0)*180.0/acos(-1.0);
  else
    tors=acos(tors)*180.0/acos(-1.0);
  return tors;
}

void Molecule::com()
{
  cofm.resize(3);
  double sum_mass=0;
  for(int i=0;i<no_atoms;i++)
  {
    for(int j=0;j<3;j++)
    {
      cofm[j]+=mass[z_vals[i]]*geom(j,i);
    }
    sum_mass+=mass[z_vals[i]];
  }
  for(int j=0;j<3;j++)
  {
    cofm[j]/=sum_mass;
  }
}

void Molecule::inertia()
{
  I.resize(3,3);
  for(int i=0;i<no_atoms;i++)
  {
    I(0,0)+=mass[z_vals[i]]*((geom(1,i)*geom(1,i))+(geom(2,i)*geom(2,i)));
    I(1,1)+=mass[z_vals[i]]*(pow(geom(0,i),2.0)+pow(geom(2,i),2.0));
    I(2,2)+=mass[z_vals[i]]*(pow(geom(1,i),2.0)+pow(geom(0,i),2.0));
  }
  for(int j=0;j<3;j++)
  {
    for(int k=j+1;k<3;k++)
    {
      for(int i=0;i<no_atoms;i++)
      {
        I(j,k)+=(mass[z_vals[i]]*geom(j,i)*geom(k,i));
      }
      I(k,j)=I(j,k);
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(3);
  es.compute(I);
  evals=es.eigenvalues();
}

int Molecule::rotor_type()
{
  if(evals(0)==evals(1))
  {
    if(evals(1)==evals(2))
    {
      return 1;
    }
    else
    {
      return 2;
    }
  }
  else if(evals(1)==evals(2))
  {
    return 3;
  }
  else if(evals(0)<pow(10,-4))
  {
    return 4;
  }
  else
    return 5;
}


Molecule::Molecule(){}
Molecule::~Molecule(){}
