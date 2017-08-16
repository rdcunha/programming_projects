#include<omp.h>
#include<iostream>

using namespace std;

int main()
{
#pragma omp parallel
    {
    cout<<omp_get_max_threads()<<endl;
    }    
    return 0;
}
