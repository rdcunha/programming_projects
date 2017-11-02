/*
Code for reading in various energy values and 2-electron integrals and iteratively calculating the HF energy using the SCF method.
Also includes a function to compute the MP2 correction to the HF energy.

Requires an input file *.cpp creating an instance of the class Hfock and calling its member function do_hf with user-set arguments.

Outputs the final Fock matrix, final orbital energies and the total SCF energy for each iteration until convergence is reached.

*/

#include "cc.h"

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
    //cout<<'\n'<<p<<'\t'<<q<<'\t'<<r<<'\t'<<s<<'\t'<<mat[index(index(p,q),index(r,s))];
  }
  input.close();
}

void Hfock::orthogonalize(Matrix &mat,Matrix &mat1)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(no_parts);
  es.compute(overlap);
  Matrix eval;
  eval=es.eigenvalues();
  Matrix evals;
  evals.resize(no_parts,no_parts);
  for(int i=0;i<no_parts;i++)
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
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(no_parts);
  es.compute(mat2);
  Matrix eval;
  Matrix evecs;
  eval=es.eigenvalues();
  evecs=es.eigenvectors();
  evecs=mat*evecs;
  //cout<<'\n'<<evecs<<'\n';
  mat3.setZero(no_parts,no_parts);
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
  if(sum>pow(10,-11) )
    return true;
  else
    return false;
}

void Hfock::mp2(Matrix &mat)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(no_parts);
  es.compute(mat);
  Matrix eigvals=es.eigenvalues();
  cout<<"Eigenvalues of final Fock matrix: "<<std::endl<<eigvals<<std::endl;
  Matrix eigvecs=es.eigenvectors();
  eigvecs=ortho*eigvecs;
  eri_trans.resize(406);
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

//                                                                   //
//                                                                   //
//Smarter algorithm starts here; done according to Hint 2, Project 4 //
//                                                                   //
//                                                                   //
  Matrix tmp1;
  Matrix tmp;
  tmp1.setZero(2*no_parts*no_occ,2*no_parts*no_occ);
  tmp.resize(no_parts,no_parts);
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
//        cout<<eri_trans[index(rs,pq)]<<std::endl;
      }
    }
  }
}

/*cout<<"Temporary matrices: "<<std::endl;
cout<<tmp<<std::endl;
cout<<tmp1<<std::endl;*/

//Computing the MP2 Energy

  double sume=0;
  for(int i=0;i<no_occ;i++)
  {
    for(int a=no_occ;a<no_parts;a++)
    {
      for(int j=0;j<no_occ;j++)
      {
        for(int b=no_occ;b<no_parts;b++)
        {
          //cout<<"eri_trans ijab: "<<eri_trans[index(index(i,a),index(j,b))]<<std::endl;
          sume += (eri_trans[index(index(i,a),index(j,b))]*(2.0*eri_trans[index(index(i,a),index(j,b))] - eri_trans[index(index(i,b),index(j,a))]))/(eigvals(i)+eigvals(j)-eigvals(a)-eigvals(b));
        }
      }
    }
  }           
  emp2 = sume;
  printf("MP2 energy: %8.12f \n", sume);
}

void Hfock::do_hf(int np,int no, string enucfile, string overlapfile, string kinfile,string nuc_attfile, string erifile)
{
  no_parts = np;
  no_occ = no;
  no_vir = no_parts - no_occ;
  ifstream input(enucfile);
  input>>enuc;
  cout<<"Nuclear energy: "<<enuc<<std::endl;

//reading in the integrals
  read_in(overlapfile,overlap);
  read_in(kinfile,kinetic);
  read_in(nuc_attfile,nuc_att);

//forming the core hamiltonian
  core_ham.resize(no_parts,no_parts);
  core_ham = kinetic + nuc_att;

//reading in the two electron integrals
  read_2e(erifile,eri);

//building the orthogonalization matrix
  orthogonalize(overlap,ortho);

//transforming the fock matrix and getting the initial density matrix
  make_density(ortho,core_ham,fock,density);

//computing the initial SCF energy
  newe=scf_energy(density,core_ham);

//building the new fock matrix
  make_fock( density, fock);
  cout<<"\nInitial Fock matrix: \n"<< fock<<'\n';
  make_density( ortho, fock, newfock, density);
  olde= newe;
  newe= scf_energy( density, fock);
  cout<<"\nSecond SCF energy: "<< newe<<'\n';

//SCF procedure
  while(std::abs(olde-newe)>std::pow(10,-12) && converge_density())
  {    
     make_fock(density,fock);
     old_density= density;
     make_density(ortho,fock,newfock,density);
     olde=newe;
     newe=scf_energy( density, fock);
    printf("Total SCF energy: %8.12f \n", newe);
  }
  
//obtaining the MO coefficients and energies from the final Fock matrix
   make_fock( density, fock);
   make_density( ortho, fock, newfock, density);

//Calculating the MP2 correction to the energy using the MO coefficients and energies
   mp2(newfock);

}
