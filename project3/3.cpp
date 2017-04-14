#include <iostream>
#include "hf.h"

int main()
{
  Hfock test;
  test.no_parts=7;
  test.no_occ=5;
  ifstream input("enuc.dat");
  input>>test.enuc;
  //cout<<test.enuc;

//reading in the integrals
  test.read_in("s.dat",test.overlap);
  test.read_in("t.dat",test.kinetic);
  test.read_in("v.dat",test.nuc_att);

//forming the core hamiltonian
  test.core_ham.resize(test.no_parts,test.no_parts);
  test.core_ham=test.kinetic + test.nuc_att;
  //cout<<'\n'<<test.core_ham<<'\n';

//reading in the two electron integrals
  test.read_2e("eri.dat",test.eri);
  //cout<<'\n'<<test.two_e<<'\n';

//building the orthogonalization matrix
  test.orthogonalize(test.overlap,test.ortho);
  //cout<<'\n'<<test.ortho<<'\n';

//transforming the fock matrix and getting the initial density matrix
  test.make_density(test.ortho,test.core_ham,test.fock,test.density);
  //cout<<'\n'<<test.fock<<'\n';
  //cout<<'\n'<<test.density<<'\n';

//computing the initial SCF energy
  test.newe=test.scf_energy(test.density,test.core_ham);
  cout<<"\nInitial SCF energy: "<<test.newe<<'\n';

//building the new fock matrix
  test.make_fock(test.density,test.fock);
  cout<<"\nInitial Fock matrix: \n"<<test.fock<<'\n';
  test.make_density(test.ortho,test.fock,test.newfock,test.density);
  test.olde=test.newe;
  test.newe=test.scf_energy(test.density,test.fock);
  cout<<"\nSecond SCF energy: "<<test.newe<<'\n';

//SCF procedure
  //Some notes:
  /*
    1. The fock matrix is built from the previous iteration's density and the core hamiltonian (make_fock)
    2. The fock matrix is converted into the orthonormal AO basis by similarity transformation (make_density)
    3. The converted fock matrix is diagonalized and its eigenvectors, transformed into the non-orthogonal AO basis, are used to form the density matrix (make_density)
    4. The SCF energy is then calculated using the density matrix, the fock matrix and the core hamiltonian, all in the non-orthogonal basis (scf_energy)
     */
  while(std::abs(test.olde-test.newe)>std::pow(10,-8) && test.converge_density())
  {    
    test.make_fock(test.density,test.fock);
    test.old_density=test.density;
    test.make_density(test.ortho,test.fock,test.newfock,test.density);
    test.olde=test.newe;
    test.newe=test.scf_energy(test.density,test.fock);
    cout<<"Total SCF energy: "<<test.newe<<'\n';
  }
  return 0;
}
