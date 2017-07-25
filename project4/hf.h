#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using std::cin;
using std::cout;
using std::vector;
using std::ifstream;
using std::string;

class Hfock
{
  public:
    Matrix overlap;
    Matrix kinetic;
    Matrix nuc_att;
    Matrix core_ham;
    Matrix ortho;
    Matrix fock;
    Matrix newfock;
    Matrix density;
    Matrix old_density;

    vector <double> eri;

    int no_parts;
    int no_occ;
    double enuc;
    double olde;
    double newe;

    void read_in(string ifile,Matrix &mat);
    void read_2e(string ifile,vector <double> &mat);
    void orthogonalize(Matrix &mat, Matrix &mat1);
    void make_density(Matrix &mat, Matrix &mat1, Matrix &mat2, Matrix &mat3);
    double scf_energy(Matrix &mat, Matrix &mat1);
    void make_fock(Matrix &mat, Matrix &mat1);
    bool converge_density();
    void mp2(Matrix &mat);

    //Hfock();
    //~Hfock();
};
