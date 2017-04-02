#include <iostream>
#include <cmath>
#include "molecule.h"

using namespace std;

int main()
{
  Molecule test;
  test.iname="geom.dat";
  test.read_in();
  //checking if the geometry was read correctly
  test.print_geom();

  //creating and printing bond lengths
  cout<<"\nBond lengths:\n";
  for(int i=0;i<test.no_atoms;i++)
  {
    for(int j=i+1;j<test.no_atoms;j++)
    {
      cout<<i<<'\t'<<j<<'\t'<<test.bond(i,j)<<'\n';
    }
  }

  //creating and printing bond angles
  cout<<"\nBond angles:\n";
  for(int i=0;i<test.no_atoms;i++)
  {
    for(int j=i+1;j<test.no_atoms;j++)
    {
      for(int k=j+1;k<test.no_atoms;k++)
      {
        if(test.bond(i,j)<4.0 && test.bond(j,k)<4.0)
          cout<<i<<'\t'<<j<<'\t'<<k<<'\t'<<test.angle(i,j,k)<<'\n';
      }
    }
  }
  
  //creating and printing out of plane angles
  cout<<"\nOut-of-plane angles:\n";
  for(int i=0;i<test.no_atoms;i++)
  {
    for(int j=0;j<test.no_atoms;j++)
    {
      for(int k=0;k<test.no_atoms;k++)
      {
        for(int l=k+1;l<test.no_atoms;l++)
        {
          if(i!=j && j!=k && k!=l && i!=l && i!=k && j!=l)
          {
            if(test.bond(j,k)<4.0 && test.bond(i,k)<4.0 && test.bond(k,l)<4.0)
            {
              cout<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\t'<<test.oop(i,j,k,l)<<'\n';
            }
          }
        }
      }
    }
  }

  //creating and printing torsional angles
  cout<<"\nTorsional angles:\n";
  for(int i=0;i<test.no_atoms;i++)
  {
    for(int j=0;j<i;j++)
    {
      for(int k=0;k<j;k++)
      {
        for(int l=0;l<k;l++)
        {
            if(test.bond(i,j)<4.0 && test.bond(j,k)<4.0 && test.bond(k,l)<4.0)
            {
              cout<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\t'<<test.torsion(i,j,k,l)<<'\n';
            }
        }
      }
    }
  }

  //finding the coordinates of the center of mass
  test.com();
  cout<<"\nCenter of mass coordinates:\n"<<test.cofm[0]<<" "<<test.cofm[1]<<" "<<test.cofm[2]<<'\n';
  //translating to center-of-mass coordinates
  test.translate(-test.cofm[0],-test.cofm[1],-test.cofm[2]);
  test.print_geom();

  //finding the inertia tensor and the principal moments of inertia
  test.inertia();
  cout<<"\nMoment of inertia tensor:\n"<<test.I<<'\n';
  cout<<"\nPrincipal moments of inertia (amu*bohr^2):\nIa = "<<test.evals(0)<<"\nIb = "<<test.evals(1)<<"\nIc = "<<test.evals(2)<<'\n'; 
  cout<<"\nPrincipal moments of inertia (amu*Ang^2):\n";
  for(int i=0;i<3;i++){cout<<test.evals(i)*(pow(0.52917721067,2)) <<'\n';}
  cout<<"\nPrincipal moments of inertia (g*cm^2):\n";
  for(int i=0;i<3;i++){cout<<test.evals(i)*(1.66053904*(pow(10,-40))) <<'\n';}
  test.evecs.resize(3,3);
  jacobi(test.I,test.evecs);

  //checking the kind of rotor
  switch(test.rotor_type())
  {
    case 1:
      cout<<"The molecule is a spherical top.\n";
    case 2:
      cout<<"The molecule is an oblate symmetric top.\n";
    case 3:
      cout<<"The molecule is a prolate symmetric top.\n";
    case 4:
      cout<<"The molecule is linear.\n";
    default:
      cout<<"The molecule is an asymmetric top.\n";
  }
  
  //finding the rotational constants in cm-1
  string con="ABC";
  double factor=6.626070040*pow(10,-34)/(pow(acos(-1),2)*8.0*29979245800);
  for(int i=0;i<3;i++)
  {
    cout<<"Rotational constant "<<con[i]<<"(cm-1) = "<<factor/(test.evals(i)*pow(0.52917721067,2)*(1.66053904*(pow(10,-47))))<<'\n';
  }
  cout<<'\n';
  //For Rotational constants in MHz, multiplying by c*10^-6
  factor*=29979.2458;
  for(int i=0;i<3;i++)
  {
    cout<<"Rotational constant "<<con[i]<<"(MHz) = "<<factor/(test.evals(i)*pow(0.52917721067,2)*(1.66053904*(pow(10,-47))))<<'\n';
  }

  return 0;
  }
