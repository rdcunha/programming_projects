#include <string>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

class Molecule
{
  public:
    string iname;
    int no_atoms;
    int charge;
    vector <int> z_vals;
    Matrix geom;
    Matrix I;
    Matrix evals;
    string point_group;
    vector <double> cofm;

    void read_in();
    void print_geom();
    void rotate(double phi);
    void translate(double x, double y, double z);
    void com();
    double bond(int a1, int a2);
    double angle(int a1, int a2, int a3);
    double torsion(int a1, int a2, int a3, int a4);
    double oop(int a1, int a2, int a3, int a4);
    void inertia();
    int rotor_type();

    Molecule();
    ~Molecule();
};
