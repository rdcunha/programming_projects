#include "hf.h"

void Hfock::read_in(string ifile,Matrix &mat)
{
  double a1=0,a2=0;
  ifstream input(ifile);
  mat.resize(no_parts,no_parts);
  //cout<<mat<<'\n';
   while(input)
   {
     for(int i=0;i<no_parts;i++)
     {
       for(int j=0;j<i+1;j++)
       {
         input>>a1>>a2>>mat(i,j);
         mat(j,i)=mat(i,j);
       }
     }
   }
   input.close();
  //cout<<'\n'<<mat<<'\n';
}

int index(int a, int b)
{
  int val=0;
  if(a>b)
  {
    val=(a*(a+1.0)/2.0 + b);
  }
  else
  {
    val=(b*(b+1.0)/2.0 + a);
  }
  return val;
}

void Hfock::read_2e(string ifile,vector <double> &mat)
{
  int p=1,q=1,r=1,s=1;
  ifstream input(ifile);
  mat.resize(26106);
  while(input)
  {
    input>>p>>q>>r>>s>>mat[index(index(p,q),index(r,s))];
    //cout<<'\n'<<p<<'\t'<<q<<'\t'<<r<<'\t'<<s<<'\t'<<mat[index(index(p,q),index(r,s))];
  }
  input.close();
}

void Hfock::orthogonalize(Matrix &mat,Matrix &mat1)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(7);
  es.compute(overlap);
  Matrix eval;
  eval=es.eigenvalues();
  Matrix evals;
  evals.resize(7,7);
  for(int i=0;i<7;i++)
  {
    evals(i,i)= 1.0/std::sqrt(eval(i));
  }
  //cout<<'\n'<<evals<<'\n';
  Matrix evecs;
  evecs=es.eigenvectors();
  mat1=evecs*evals*evecs.transpose();
}

void Hfock::make_density(Matrix &mat, Matrix &mat1, Matrix &mat2, Matrix &mat3)
{
  mat2=mat.transpose()*mat1*mat;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(7);
  es.compute(mat2);
  Matrix eval;
  Matrix evecs;
  eval=es.eigenvalues();
  evecs=es.eigenvectors();
  evecs=mat*evecs;
  //cout<<'\n'<<evecs<<'\n';
  mat3.resize(no_parts,no_parts);
  for(int i=0;i<no_parts;i++)
  {
    for(int j=0;j<no_parts;j++)
    {
      mat3(i,j)=0;
      for(int m=0;m<no_occ;m++)
      {
        mat3(i,j)+=evecs(i,m)*evecs(j,m);
      }
    }
  }
}

double Hfock::scf_energy(Matrix &mat, Matrix &mat1)
{
  double energy=0;
  for(int i=0;i<no_parts;i++)
  {
    for(int j=0;j<no_parts;j++)
    {
      energy+= mat(i,j)*(core_ham(i,j)+mat1(i,j));
    }
  }
  energy+=enuc;
  return energy;
}

void Hfock::make_fock(Matrix &mat, Matrix &mat1)
{
  mat1=core_ham;
  for(int p=1;p<no_parts+1;p++)
  {
    for(int q=1;q<no_parts+1;q++)
    {
      for(int r=1;r<no_parts+1;r++)
      {
        for(int s=1;s<no_parts+1;s++)
        {
          mat1(p-1,q-1)+=mat(r-1,s-1)*(2.0*eri[index(index(p,q),index(r,s))] - eri[index(index(p,s),index(r,q))]);
        }
      }
    }
  }
}

bool Hfock::converge_density()
{
  old_density.resize(no_parts,no_parts);
  double sum=0;
  for(int i=0;i<no_parts;i++)
  {
    for(int j=0;j<no_parts;j++)
    {
      sum+=pow(density(i,j)-old_density(i,j),2.0);
    }
  }
  sum=std::sqrt(sum);
  cout<<"RMS(D) "<<sum<<std::endl;
  if(sum>pow(10,-8) )
    return true;
  else
    return false;
}

