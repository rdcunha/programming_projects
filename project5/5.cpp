#include "cc.h"


int main()
{
    Hfock h2o;
    h2o.do_cc(7,5,"enuc.dat","s.dat","t.dat","v.dat","eri.dat");
    return 0;
}
