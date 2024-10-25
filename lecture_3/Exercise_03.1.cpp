#include <iostream>
#include <cmath>
#include <iomanip>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
#include <armadillo>
#include <fstream>
#include <string>
#include "../classes/inc/block_average.h"
using namespace std;
using namespace arma;

int main(){
    // initialize variables
    int M=1e5;        // number of throws
    int N=800;         // number of blocks
    int N_step=1e4;    // number of step for computing disctretized call price
    int L=M/N;         // number of throws in each block
    double S0=100.;     // asset price at t=0
    double T=1.;        // delivery time
    double K=100.;      // strike price
    double r=0.1;      // risk-free interest rate
    double sigma=0.25; // volatility
    double mu=0.;       // mean value

    // Call price direct
    bool direct=true;
    call_price C_direct(S0,T,K,r,sigma,mu,direct);
    string ofile_Cdirect="Results/ex_03.1_call_option_price_direct.dat";
    C_direct.Average(N,L,ofile_Cdirect);
    // Call price discrete
    bool discrete=false;
    call_price C_discrete(S0,T,K,r,sigma,mu,discrete);
    string ofile_Cdiscrete="Results/ex_03.1_call_option_price_discrete.dat";
    C_discrete.Set_Nstep(N_step);
    C_discrete.Average(N,L,ofile_Cdiscrete);
    //put price direct
    put_price P_direct(S0,T,K,r,sigma,mu,direct);
    string ofile_Pdirect="Results/ex_03.1_put_option_price_direct.dat";
    P_direct.Average(N,L,ofile_Pdirect);
    //put price discrete
    put_price P_discrete(S0,T,K,r,sigma,mu,discrete);
    string ofile_Pdiscrete="Results/ex_03.1_put_option_price_discrete.dat";
    P_discrete.Set_Nstep(N_step);
    P_discrete.Average(N,L,ofile_Pdiscrete);
    return 0;
}
