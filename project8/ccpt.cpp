#include "cc.h"

double Hfock::make_tc_ijkabc(int i, int j, int k, int a, int b, int c)
{
    double tc = 0;
    for(int e=0; e < 2*no_vir; e++)
    {
        tc += t_ijab(j, k, a, e) * spin_eri((2*no_occ)+e, i, (2*no_occ)+b, (2*no_occ)+c);
    }
    for(int m=0; m < 2*no_occ; m++)
    {
        tc -= t_ijab(i, m, b, c) * spin_eri(m, (2*no_occ)+a, j, k);
    }
    return tc;
}

double Hfock::make_td_ijkabc(int i, int j, int k, int a, int b, int c)
{
    double td = 0;
    td = t_ia(i, a) * spin_eri(j, k, (2*no_occ)+b, (2*no_occ)+c);
    return td;
}

double Hfock::permute1(int i, int j, int k, int a, int b, int c, bool tc)
{
    double val = 0;
    if(tc)
    {
    val = make_tc_ijkabc(i, j, k, a, b, c) - make_tc_ijkabc(j, i, k, a, b, c) - make_tc_ijkabc(k, j, i, a, b, c);
    }
    else
    {
    val = make_td_ijkabc(i, j, k, a, b, c) - make_td_ijkabc(j, i, k, a, b, c) - make_td_ijkabc(k, j, i, a, b, c);
    }
    return val;
}

double Hfock::permute(int i, int j, int k, int a, int b, int c, bool tc)
{
    double val = 0;
    val = permute1(i, j, k, a, b, c, tc) - permute1(i, j, k, b, a, c, tc) - permute1(i, j, k, c, b, a, tc);
    return val;
}


void Hfock::do_ccpt(int np,int no, string enucfile, string overlapfile, string kinfile,string nuc_attfile, string erifile)
{
    do_cc(np, no, enucfile, overlapfile, kinfile, nuc_attfile, erifile);
    double pt_energy = 0;
    /*
       Code here will calculate triples amplitudes and add contributions to the pt_energy correction
    */
    double tc_ijkabc = 0;
    double td_ijkabc = 0;
    double d_ijkabc = 0;
    for(int i=0; i < 2*no_occ; i++){
        for(int j=0; j < 2*no_occ; j++){
            for(int k=0; k < 2*no_occ; k++){
                for(int a=0; a < 2*no_vir; a++){
                    for(int b=0; b < 2*no_vir; b++){
                        for(int c=0; c < 2*no_vir; c++){
                            tc_ijkabc = 0;
                            td_ijkabc = 0;
                            d_ijkabc = 0;
                            //calculate tc_ijkabc 
                            tc_ijkabc = permute(i, j, k, a, b, c, true);
                            cout<<"connected t i: "<<i<<" j: "<<j<<" k: "<<k<<" a: "<<a<<" b: "<<b<<" c: "<<c<<" "<<tc_ijkabc<<std::endl;
                            td_ijkabc = permute(i, j, k, a, b, c, false);
                            //calculate denominator
                            d_ijkabc = fock(i,i) + fock(j,j) + fock(k,k) - fock((2*no_occ)+a,(2*no_occ)+a) - fock((2*no_occ)+b,(2*no_occ)+b) - fock((2*no_occ)+c,(2*no_occ)+c);
                            //divide by the denominator
                            tc_ijkabc /= d_ijkabc;
                            td_ijkabc /= d_ijkabc;

                            pt_energy += tc_ijkabc*d_ijkabc*(tc_ijkabc + td_ijkabc);
                        }
                    }
                }
            }
        }
    }
    pt_energy /= 36.0;
    printf("CCSD(T) correction to the energy: %6.12f \n", pt_energy);
    printf("Total energy: %6.12f \n", cc_e_new+pt_energy+newe);
}
