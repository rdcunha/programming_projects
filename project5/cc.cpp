/*
Code for computing the Coupled Cluster Singles and Doubles energy using an SCF converged HF reference wavefunction.

Requires an input file *.cpp which creates an instance of the class Hfock and calls its member function do_cc with user-set arguments.

Outputs the Fock matrix in the spin orbital basis as well as the CCSD energy at each iteration until convergence.
*/

#include "cc.h"

//new index function
int indexing(int a, int b)
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

//transforming the 2-electron MO-basis ints into the spin orbital basis
void Hfock::transform_ints(vector<double> vec, Eigen::Tensor<double,4> &ten)
{
  int ik,jl,ikjl,il,jk,iljk;
  for(int i=0;i<(2*no_parts);i++)
  {
    for(int j=0;j<(2*no_parts);j++)
    {
      //cout<<"ij = "<<ij<<std::endl;
      for(int k=0;k<(2*no_parts);k++)
      {
        for(int l=0;l<(2*no_parts);l++)
        {
          ik=indexing(i/2,k/2);
          jl=indexing(j/2,l/2);
          //cout<<"kl = "<<kl<<std::endl;
          il=indexing(i/2,l/2);
          jk=indexing(j/2,k/2);

          ikjl=indexing(ik,jl);
          iljk=indexing(il,jk);
          //cout<<"ijkl = "<<ijkl<<std::endl;
          ten(i,j,k,l)=vec[ikjl]*(i%2==k%2)*(j%2==l%2) - vec[iljk]*(i%2==l%2)*(j%2==k%2);
        }
      }
    }
  }
  //int i=13,j=13,k=13,l=13;
  //cout<<"ERI_TRANS value at 13,13,13,13 = "<<vec[indexing(indexing(7,7),indexing(7,7))];
}

//creating the new spin-orbital basis Fock matrix
Matrix Hfock::make_fock_new(Eigen::Tensor<double,4> &ten)
{
  Matrix mat;
  mat.setZero(14,14);
  Matrix core_ham_mo(7,7);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(7);
  es.compute(newfock);
  Matrix eigvecs = es.eigenvectors();
  eigvecs = ortho*eigvecs;
  core_ham_mo = eigvecs.transpose()*core_ham*eigvecs;

  for(int p=0;p<(2*no_parts);p++)
  {
    for(int q=0;q<(2*no_parts);q++)
    {
      mat(p,q)=core_ham_mo(p/2,q/2)*(p%2 == q%2);
      for(int m=0;m<(2*no_occ);m++)
      {
        mat(p,q)+=ten(p,m,q,m);
      }
      if(std::abs(mat(p,q)) < std::pow(10,-7))
        {
        mat(p,q) =0.0;
        }
    }
  }
  cout<<"Fock matrix in spin orbital basis:"<<std::endl<<mat<<std::endl;

/*  fock = eigvecs.transpose()*fock*eigvecs;
  cout<<"Fock matrix cheating in spin orbital basis:"<<std::endl<<fock<<std::endl;*/
  return mat;
  
}

//creating initial t-amplitudes and the MP2 energy as a check
Eigen::Tensor<double, 4> Hfock::initial_amps(Matrix &mat, Eigen::Tensor<double,4> &ten)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(14);
  es.compute(mat);
  Matrix eigvals=es.eigenvalues();
  //cout<<"eigenvalues of the fock matrix:"<<std::endl;
  //cout<<eigvals;

  Eigen::Tensor<double,4> t(10,10,4,4);
  double Emp2=0;
  //cout<<"t(i,j,a,b)'s"<<std::endl;
  for(int i=0;i<(2*no_occ);i++)
  {
    for(int j=0;j<(2*no_occ);j++)
    {
      for(int a=0;a<4;a++)
      {
        for(int b=0;b<4;b++)
        {
          t(i,j,a,b)=(ten(i,j,a+(2*no_occ),b+(2*no_occ)))/(eigvals(i)+eigvals(j)-eigvals((a+(2*no_occ)))-eigvals((b+(2*no_occ))));
        }
      }
    }
  }
  for(int i=0;i<2*no_occ;i++)
  {
    for(int j=0;j<2*no_occ;j++)
    {
      for(int a=0;a<4;a++)
      {
        for(int b=0;b<4;b++)
        {
          Emp2+=t(i,j,a,b)*(ten(i,j,a+(2*no_occ),b+(2*no_occ)));
          //cout<<"MP2 Energy: "<<Emp2<<std::endl;
        }
      }
    }
  }
  Emp2*=0.25;
  //cout<<"MP2 Energy using spin orbitals: "<<Emp2<<std::endl;
  printf("MP2 energy using spin orbitals: %8.12f \n", Emp2);
  return t;
}

std::pair<Eigen::Tensor<double,4>, Eigen::Tensor<double,4> > Hfock::make_tau(Matrix &mat, Eigen::Tensor<double,4> &ten)
{
    Eigen::Tensor<double, 4> tau;
    Eigen::Tensor<double, 4> tau_t;
    tau = ten;
    tau_t = ten;
    for(int i=0; i < 2*no_occ; i++)
    {
        for(int j=0; j < 2*no_occ; j++)
        {
            for(int a=0; a < 2*(no_parts - no_occ); a++)
            {
                for(int b=0; b < 2*(no_parts - no_occ); b++)
                {
                    tau(i,j,a,b) += mat(i,a)*mat(j,b);
                    tau(i,j,a,b) -= mat(i,b)*mat(j,a);
                    tau_t(i,j,a,b) += 0.5*mat(i,a)*mat(j,b);
                    tau_t(i,j,a,b) -= 0.5*mat(i,b)*mat(j,a);
                }
            }
        }
    }
/*    cout<<"t_ijab:"<<std::endl;
    cout<<ten<<std::endl;
    cout<<"Tau:"<<std::endl;
    cout<<tau<<std::endl;
    cout<<"Tau_t:"<<std::endl;
    cout<<tau_t<<std::endl;*/
    return std::make_pair(tau,tau_t);
}

double Hfock::make_F_ae(int a, int e)
{
    double f_ae;
    double temp1 = 0;
    double temp2 = 0;
    double temp3 = 0;
    for(int m=0; m < 2*no_occ; m++)
    {
        temp1 += fock(m,(2*no_occ)+e) * t_ia(m,a);
        for(int f=0; f < 2*(no_vir); f++)
        {
            temp2 += t_ia(m,f) * spin_eri(m,(2*no_occ)+a,(2*no_occ)+f,(2*no_occ)+e);
            for(int n=0; n < 2*no_occ; n++)
            {
                temp3 += tau_t(m,n,a,f) * spin_eri(m,n,(2*no_occ)+e,(2*no_occ)+f);
            }
        }
    }
    f_ae = (1- (a==e)) * fock((2*no_occ)+a,(2*no_occ)+e) - 0.5 * temp1 + temp2 - 0.5*temp3;
    return f_ae;
}

double Hfock::make_F_mi(int m, int i)
{
    double f_mi;
    double temp1 = 0;
    double temp2 = 0;
    double temp3 = 0;
    for(int e=0; e < 2*(no_vir); e++)
    {
        temp1 += fock(m,(2*no_occ)+e) * t_ia(i,e);
        for(int n=0; n < 2*no_occ; n++)
        {
            temp2 += t_ia(n,e) * spin_eri(m,n,i,(2*no_occ)+e);
            for(int f=0; f < 2*(no_vir); f++)
            {
                temp3 += tau_t(i,n,e,f) * spin_eri(m,n,(2*no_occ)+e,(2*no_occ)+f);
            }
        }
    }
    f_mi = (1- (m==i)) * fock(m,i) + 0.5 * temp1 + temp2 + 0.5*temp3;
    return f_mi;
}

double Hfock::make_F_me(int m, int e)
{
    double f_me = fock(m,(2*no_occ)+e);
    for(int n=0; n < 2*no_occ; n++)
    {
        for(int f=0; f < 2*(no_vir); f++)
        {
            f_me += t_ia(n,f) * spin_eri(m,n,(2*no_occ)+e,(2*no_occ)+f);
        }
    }
    return f_me;
}

double Hfock::make_W_mnij(int m, int n, int i, int j)
{
    double w_mnij;
    double temp1 = 0;
    double temp2 = 0;
    for(int e=0; e < 2*(no_vir); e++)
    {
        temp1 += t_ia(j,e) * spin_eri(m,n,i,(2*no_occ)+e);
        temp1 -= t_ia(i,e) * spin_eri(m,n,j,(2*no_occ)+e);
        for(int f=0; f < 2*(no_vir); f++)
        {
            temp2 += tau(i,j,e,f) * spin_eri(m,n,(2*no_occ)+e,(2*no_occ)+f);
        }
    }
    w_mnij = spin_eri(m,n,i,j) + temp1 + 0.25 * temp2;
    return w_mnij;
}

double Hfock::make_W_abef(int a, int b, int e, int f)
{
    double w_abef;
    double temp1 = 0;
    double temp2 = 0;
    for(int m=0; m < 2*no_occ; m++)
    {
        temp1 += t_ia(m,b) * spin_eri((2*no_occ)+a,m,(2*no_occ)+e,(2*no_occ)+f);
        temp1 -= t_ia(m,a) * spin_eri((2*no_occ)+b,m,(2*no_occ)+e,(2*no_occ)+f);
        for(int n=0; n < 2*no_occ; n++)
        {
            temp2 += tau(m,n,a,b) * spin_eri(m,n,(2*no_occ)+e,(2*no_occ)+f);
        }
    }
    w_abef = spin_eri((2*no_occ)+a,(2*no_occ)+b,(2*no_occ)+e,(2*no_occ)+f) - temp1 + 0.25 * temp2;
    return w_abef;
}

double Hfock::make_W_mbej(int m, int b, int e, int j)
{
    double w_mbej;
    double temp1 = 0;
    double temp2 = 0;
    double temp3 = 0;
    for(int f=0; f < 2*(no_vir); f++)
    {
        temp1 += t_ia(j,f) * spin_eri(m,(2*no_occ)+b,(2*no_occ)+e,(2*no_occ)+f);
    }

    for(int n=0; n < 2*no_occ; n++)
    {
        temp2 += t_ia(n,b) * spin_eri(m,n,(2*no_occ)+e,j);
        for(int f=0; f < 2*(no_vir); f++)
        {
            temp3 += (0.5 * t_ijab(j,n,f,b) + t_ia(j,f) * t_ia(n,b) ) * spin_eri(m,n,(2*no_occ)+e,(2*no_occ)+f);
        }
    }
    w_mbej = spin_eri(m,(2*no_occ)+b,(2*no_occ)+e,j) + temp1 - temp2 - temp3;
    return w_mbej;
}

Matrix Hfock::update_t_ia()
{
    Matrix up_tia(2*no_occ,2*(no_vir));
    double temp1=0, temp2=0, temp3=0, temp4=0, temp5=0, temp6=0;
    for(int i=0; i < 2*no_occ; i++)
    {
        for(int a=0; a < 2*(no_vir); a++)
        {
            temp1=0;
            for(int e=0; e < 2*(no_vir); e++)
            {
                temp1 += t_ia(i,e) * make_F_ae(a,e);
            }
            temp2=0;
            temp3=0;
            temp4=0;
            temp5=0;
            for(int m=0; m < 2*(no_occ); m++)
            {
                temp2 += t_ia(m,a) * make_F_mi(m,i);
                for(int e=0; e < 2*(no_vir); e++)
                {
                    temp3 += t_ijab(i,m,a,e) * make_F_me(m,e);
                    for(int f=0; f < 2*(no_vir); f++)
                    {
                        temp4 += t_ijab(i,m,e,f) * spin_eri(m,(2*no_occ)+a,(2*no_occ)+e,(2*no_occ)+f);
                    }
                    for(int n=0; n < 2*no_occ; n++)
                    {
                        temp5 += t_ijab(m,n,a,e) * spin_eri(n,m,(2*no_occ)+e,i);
                    }
                }
            }
            temp6=0;
            for(int n=0; n< 2*no_occ; n++)
            {
                for(int f=0; f < 2*(no_vir); f++)
                {
                    temp6 += t_ia(n,f) * spin_eri(n,(2*no_occ)+a,i,(2*no_occ)+f);
                }
            }
            up_tia(i,a) = fock(i,(2*no_occ) + a) + temp1 - temp2 + temp3 - temp6 - 0.5 * temp4 - 0.5 * temp5;
            up_tia(i,a) /= (fock(i,i) - fock((2*no_occ)+a,(2*no_occ)+a) );
        }
    }
    return up_tia;
}

Eigen::Tensor<double, 4> Hfock::update_t_ijab()
{
    Eigen::Tensor<double, 4> up_tijab(2*no_occ, 2*no_occ,2*(no_vir),2*(no_vir));
    double temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    for(int i=0; i < 2*no_occ; i++)
    {
        for(int j=0; j < 2*no_occ; j++)
        {
            for(int a=0; a < 2*(no_parts - no_occ); a++)
            {
                for(int b=0; b < 2*(no_parts - no_occ); b++)
                {
                    temp1=0;
                    temp4=0;
                    temp6=0;
                    for(int e=0; e < 2*(no_parts - no_occ); e++)
                    {
                        double sumb = 0;
                        double suma = 0;
                        for(int m=0; m < 2*no_occ; m++)
                        {
                            sumb+= t_ia(m,b) * make_F_me(m,e);
                            suma+= t_ia(m,a) * make_F_me(m,e);
                        }
                        temp1+= t_ijab(i,j,a,e) * (make_F_ae(b,e) - 0.5 * sumb);
                        temp1-= t_ijab(i,j,b,e) * (make_F_ae(a,e) - 0.5 * suma);
                        for(int f=0; f < 2*(no_parts - no_occ); f++)
                        {
                            temp4 += tau(i,j,e,f) * make_W_abef(a,b,e,f);
                        }
                        temp6 += t_ia(i,e) * spin_eri((2*no_occ)+a,(2*no_occ)+b,(2*no_occ)+e,j);       
                        temp6 -= t_ia(j,e) * spin_eri((2*no_occ)+a,(2*no_occ)+b,(2*no_occ)+e,i);       
                    }
                    temp2=0;
                    temp3=0;
                    temp5=0;
                    temp7=0;
                    for(int m=0; m < 2*no_occ; m++)
                    {
                        //cout<<"works in m loop "<<m<<std::endl;
                        double sumj = 0;
                        double sumi = 0;
                        for(int e=0; e < 2*(no_vir); e++)
                        {
                         //   cout<<"works in e loop"<<std::endl;
                            sumj+= t_ia(j,e) * make_F_me(m,e);
                            sumi+= t_ia(i,e) * make_F_me(m,e);
                        }
                        //cout<<"works here"<<m<<std::endl;
                        temp2+= t_ijab(i,m,a,b) * (make_F_mi(m,j) + 0.5 * sumj);
                        temp2-= t_ijab(j,m,a,b) * (make_F_mi(m,i) + 0.5 * sumi);
                        for(int n=0; n < 2*no_occ; n++)
                        {
                            temp3 += tau(m,n,a,b) * make_W_mnij(m,n,i,j);
                        }
                        for(int e=0; e < 2*(no_vir); e++)
                        {
                            temp5 += t_ijab(i,m,a,e) * make_W_mbej(m,b,e,j) - t_ia(i,e) * t_ia(m,a) * spin_eri(m,(2*no_occ)+b,(2*no_occ)+e,j);
                            temp5 -= t_ijab(i,m,b,e) * make_W_mbej(m,a,e,j) - t_ia(i,e) * t_ia(m,b) * spin_eri(m,(2*no_occ)+a,(2*no_occ)+e,j);
                            temp5 -= t_ijab(j,m,a,e) * make_W_mbej(m,b,e,i) - t_ia(j,e) * t_ia(m,a) * spin_eri(m,(2*no_occ)+b,(2*no_occ)+e,i);
                            temp5 += t_ijab(j,m,b,e) * make_W_mbej(m,a,e,i) - t_ia(j,e) * t_ia(m,b) * spin_eri(m,(2*no_occ)+a,(2*no_occ)+e,i);
                        }
                        temp7 += t_ia(m,a) * spin_eri(m,(2*no_occ)+b,i,j);
                        temp7 -= t_ia(m,b) * spin_eri(m,(2*no_occ)+a,i,j);
                    }
                    up_tijab(i,j,a,b) = spin_eri(i,j,(2*no_occ)+a, (2*no_occ)+b) + temp1 - temp2 + 0.5 * temp3 + 0.5 * temp4 + +temp5 + temp6 - temp7;
                    up_tijab(i,j,a,b) /= fock(i,i) + fock(j,j) - fock((2*no_occ)+a,(2*no_occ)+a) - fock((2*no_occ)+b,(2*no_occ)+b);
                }
            }
        }
    }
    return up_tijab;
}

double Hfock::cc_energy(Matrix &t_ia, Eigen::Tensor<double,4> &t_ijab)
{
    double e_cc = 0;
    double temp1=0, temp2=0, temp3=0;
    for(int i=0; i < 2*no_occ; i++)
    {
        for(int a=0; a < 2*(no_parts - no_occ); a++)
        {
            temp1 += fock(i,(2*no_occ)+a) * t_ia(i,a);
            for(int j=0; j < 2*no_occ; j++)
            {
                for(int b=0; b < 2*(no_parts - no_occ); b++)
                {
                    temp2 += spin_eri(i,j,(2*no_occ)+a, (2*no_occ)+b) * t_ijab(i,j,a,b);
                    temp3 += spin_eri(i,j,(2*no_occ)+a, (2*no_occ)+b) * t_ia(i,a) * t_ia(j,b);
                }
            }
        }
    }
    e_cc = temp1 + 0.25 * temp2 + 0.5 * temp3;
    return e_cc;
}   


bool Hfock::converge(Eigen::Tensor<double,4> &ten1, Eigen::Tensor<double,4> &ten2)
{
    bool conv = false;
    Eigen::Tensor<double,4> diff2(no_occ, no_occ, no_vir, no_vir);
    diff2 = ten2 - ten1;
    diff2 = diff2.square();
    Eigen::Tensor<double,0> norm2 = diff2.mean(); //.sqrt();
    //cout<<"Here's norm2 "<<norm2<<std::endl;
    if(std::abs(cc_e_new - cc_e) < std::pow(10,-12) )
    {
        Matrix diff1;
        diff1 = t_ia_new - t_ia;
       // printf("t_ia difference : %8.10e \n", diff1.norm() );
       // printf("t_ijab difference : %8.10e \n", norm2(0) );
        if(std::abs(diff1.norm() ) < std::pow(10,-10) && std::abs(norm2(0) ) < std::pow(10,-10) )
        {
         conv = true;   
        }
    }
    return conv;
}

void Hfock::do_cc(int np,int no, string enucfile, string overlapfile, string kinfile,string nuc_attfile, string erifile)
{
  //do_hf reads in and sets: nuclear repulsion energy, two-electron integrals, no_parts(number of orbitals), no_occ(occupied orbitals), no_vir(virtual orbitals)
  do_hf(np, no, enucfile, overlapfile, kinfile, nuc_attfile, erifile);
  
//Initializing a tensor to store the spin orbital basis 2-electron integrals
  spin_eri = Eigen::Tensor<double,4> (14,14,14,14);
//Transforming the MO-basis 2-electron integrals
  transform_ints(eri_trans, spin_eri);

//Creating the spin orbital basis fock matrix using the transformed integrals
  fock = make_fock_new(spin_eri);

//Creating initial T1s and T2s
  t_ijab = initial_amps(fock, spin_eri);
  t_ia = Matrix::Zero(10, 4);
  //cout<<t_ia<<std::endl;

//Initializing a std::pair and storing the effective excitation operators tau and tau_t
  std::pair<Eigen::Tensor<double,4>, Eigen::Tensor<double,4> > taus = make_tau(t_ia, t_ijab);
  tau = taus.first;
  tau_t = taus.second;

//Updating the T1s and T2s using eqns 1-13 from Stanton, Bartlett et al, 1991
  t_ia_new = update_t_ia();
  t_ijab_new = update_t_ijab();

//Setting the old T1s and T2s to be the newly calculated T1s and T2s
  t_ia = t_ia_new;
  t_ijab = t_ijab_new;

//Calculating the first iteration CC energy
  cc_e = cc_energy(t_ia,t_ijab);
  printf("CC Energy iter 1 : %6.12f \n", cc_e);

//maximum number of iterations
  int maxiter = 100;

//Iterating until convergence is reached.
// 1. Update taus
// 2. Update T1s using the previous iteration's T1s and T2s
// 3. Update T2s using the previous iteration's T1s and T2s    
// 4. Calculate CC energy
// 5. Check for convergence; if not, set the old T1s and T2s to the calculated T1s and T2s and the old energy to the calculated energy
  for(int i=1; i < maxiter; i++) 
  {
      taus = make_tau(t_ia, t_ijab);
      tau = taus.first;
      tau_t = taus.second;
      t_ia_new = update_t_ia();
      t_ijab_new = update_t_ijab();
      cc_e_new = cc_energy(t_ia_new,t_ijab_new);
      printf("CC Energy iter %i : %6.12f \n", i+1, cc_e_new);
      if(converge(t_ijab,t_ijab_new) )
      {
          cout<<"Convergence reached. yay party."<<std::endl;
          printf("Total energy = %6.12f \n", newe+cc_e_new);
          break;
      }
      t_ia = t_ia_new;
      t_ijab = t_ijab_new;
      cc_e = cc_e_new;
  }
}
