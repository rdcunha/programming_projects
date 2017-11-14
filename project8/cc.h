#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>

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
    Matrix t_ia;
    Matrix t_ia_new;
    Matrix err;
    Matrix B;

    Eigen::Tensor<double, 4> t_ijab;
    Eigen::Tensor<double, 4> t_ijab_new;
    Eigen::Tensor<double, 4> spin_eri;
    Eigen::Tensor<double, 4> tau_t;
    Eigen::Tensor<double, 4> tau;
    Eigen::Tensor<double, 3> fock_list;
    Eigen::Tensor<double, 3> err_list;

    vector <double> eri;
    vector <double> eri_trans;

    int no_parts;
    int no_occ;
    int no_vir;
    double enuc;
    double olde;
    double newe;
    double emp2;
    double cc_e;
    double cc_e_new;

    void read_in(string ifile,Matrix &mat);
    void read_2e(string ifile,vector <double> &mat);
    void orthogonalize(Matrix &mat, Matrix &mat1);
    void make_density(Matrix &mat, Matrix &mat1, Matrix &mat2, Matrix &mat3);
    double scf_energy(Matrix &mat, Matrix &mat1);
    void make_fock(Matrix &mat, Matrix &mat1);
    bool converge_density();
    void mp2(Matrix &mat);
    void do_hf(int np,int no, string enucfile, string overlapfile, string kinfile,string nuc_attfile, string erifile);
    void transform_ints(vector<double> vec, Eigen::Tensor<double, 4> &ten);
    Matrix make_fock_new(Eigen::Tensor<double, 4> &ten);
    Eigen::Tensor<double, 4> initial_amps(Matrix &mat, Eigen::Tensor<double, 4> &ten);
    std::pair<Eigen::Tensor<double, 4>, Eigen::Tensor<double, 4> > make_tau(Matrix &mat, Eigen::Tensor<double, 4> &ten);
    double make_F_ae(int a, int e);
    double make_F_mi(int m, int i);
    double make_F_me(int m, int e);
    double make_W_mnij(int m, int n, int i, int j);
    double make_W_abef(int a, int b, int e, int f);
    double make_W_mbej(int m, int b, int e, int j);
    Matrix update_t_ia();
    Eigen::Tensor<double, 4> update_t_ijab();
    double cc_energy(Matrix &t_ia, Eigen::Tensor<double,4> &t_ijab);
    bool converge(Eigen::Tensor<double,4> &ten1, Eigen::Tensor<double,4> &ten2);
    void do_cc(int np,int no, string enucfile, string overlapfile, string kinfile,string nuc_attfile, string erifile);
    void do_ccpt(int np,int no, string enucfile, string overlapfile, string kinfile,string nuc_attfile, string erifile);
    double make_tc_ijkabc(int i, int j, int k, int a, int b, int c);
    double make_td_ijkabc(int i, int j, int k, int a, int b, int c);
    double permute1(int i, int j, int k, int a, int b, int c, bool tc);
    double permute(int i, int j, int k, int a, int b, int c, bool tc);
    Matrix make_e(int i);
    int update_fock_list(int fock_count);
    bool converge_err(int fock_count);
    

    //Hfock();
    //~Hfock();
};
