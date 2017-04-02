#include <iostream>
#include <fstream>
#include "molecule.h"

using namespace std;

int main()
{
  Molecule ho;
  ho.iname = "benzene_geom.txt";
  ho.read_in();
  ho.mass_weight("benzene_hess.dat");
  ho.freq();
  return 0;
}
