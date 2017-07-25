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
    val=(a*(a+1)/2 + b);
  }
  else
  {
    val=(b*(b+1)/2 + a);
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
    cout<<'\n'<<p<<'\t'<<q<<'\t'<<r<<'\t'<<s<<'\t'<<mat[index(index(p,q),index(r,s))];
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

void Hfock::mp2(Matrix &mat)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(7);
  es.compute(mat);
  Matrix eigvals=es.eigenvalues();
  cout<<"Eigenvalues of final Fock matrix: "<<std::endl<<eigvals<<std::endl;
  Matrix eigvecs=es.eigenvectors();
  eigvecs=ortho*eigvecs;
  vector <double> eri_trans(406);
  double sum=0,sum1=0,sum2=0,sum3=0;
  //cout<<"Transformed 2-electron integrals: "<<std::endl;
  int ijkl;
//commented code is the working noddy algorithm
//uncommented code is the smarter algorithm
/*  for(int i=1;i<no_parts+1;i++)
  {
    for(int j=1;j<=i;j++)
    {
      for(int k=1;k<=i;k++)
      {
        for(int l=1;l<= (i==k ? j : k);l++,ijkl++)
        {
          sum3=0;
          //sum=0;
          for(int p=1;p<no_parts+1;p++)
          {
            sum2=0;
            for(int q=1;q<no_parts+1;q++)
            {
              sum1=0;
              for(int r=1;r<no_parts+1;r++)
              {
                sum=0;
                for(int s=1;s<no_parts+1;s++)
                {
                  sum+=eigvecs(s-1,l-1)*eri[index(index(p,q),index(r,s))];
                }
                sum1+=eigvecs(r-1,k-1)*sum;
              }
              sum2+=eigvecs(q-1,j-1)*sum1;
            }
            sum3+=eigvecs(p-1,i-1)*sum2;
          }
          eri_trans[ijkl]=sum3;
          //eri_trans[ijkl]=sum;
          if(std::abs(eri_trans[ijkl])<pow(10.0,-14))
          {
            eri_trans[ijkl]=0.0;
          }
          cout<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\t'<<eri_trans[ijkl]<<std::endl;
        }
      }
    }
  }*/
 //cout<<"Eigenvals: "<<eigvals<<std::endl;
  Matrix tmp1;
  Matrix tmp;
  tmp1.setZero(28,28);
  tmp.resize(7,7);
for(int p=1,pq=0;p<no_parts+1;p++)
{
  for(int q=1;q<=p;q++,pq++)
  {
    for(int r=1;r<no_parts+1;r++)
    {
      for(int s=1;s<=r;s++)
      {
        //cout<<index(index(p,q),index(r,s))<<std::endl;
        tmp(r-1,s-1)=tmp(s-1,r-1)=eri[index(index(p,q),index(r,s))];
      }
    }
    tmp=eigvecs.transpose()*tmp*eigvecs;
    for(int r=0,rs=0;r<no_parts;r++)
    {
      for(int s=0;s<=r;s++,rs++)
      {
        tmp1(rs,pq)=tmp(r,s);
      }
    }
    tmp.setZero(no_parts,no_parts);
  }
}

for(int r=0,rs=0;r<no_parts;r++)
{
  for(int s=0;s<=r;s++,rs++)
  {
    for(int p=0,pq=0;p<no_parts;p++)
    {
      for(int q=0;q<=p;q++,pq++)
      {
        tmp(p,q)=tmp(q,p)=tmp1(rs,pq);
      }
    }
    tmp=eigvecs.transpose()*tmp*eigvecs;
    for(int p=0,pq=0;p<no_parts;p++)
    {
      for(int q=0;q<=p;q++,pq++)
      {
        eri_trans[index(rs,pq)]=tmp(p,q);
        //cout<<eri_trans[index(rs,pq)]<<std::endl;
      }
    }
  }
}

/*cout<<"Temporary matrices: "<<std::endl;
cout<<tmp<<std::endl;
cout<<tmp1<<std::endl;*/

  double sume=0;
  for(int i=0;i<5;i++)
  {
    for(int a=5;a<7;a++)
    {
      for(int j=0;j<5;j++)
      {
        for(int b=5;b<7;b++)
        {
          //cout<<"eri_trans ijab: "<<eri_trans[index(index(i,a),index(j,b))]<<std::endl;
          sume += (eri_trans[index(index(i,a),index(j,b))]*(2.0*eri_trans[index(index(i,a),index(j,b))] - eri_trans[index(index(i,b),index(j,a))]))/(eigvals(i)+eigvals(j)-eigvals(a)-eigvals(b));
        }
      }
    }
  }           
  cout<<"MP2 energy: "<<sume<<std::endl;
}
