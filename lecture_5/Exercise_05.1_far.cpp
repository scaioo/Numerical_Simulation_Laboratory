#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<armadillo>
#include <iomanip>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
#include "../classes/inc/func.h"
#include "../classes/inc/block_average.h"
#include "../classes/inc/func.h"
using namespace std;
using namespace arma;

void simulation_far(vec r_start, string filename, int i, bool gaussian,bool ground,string stringname, bool equilibation){
    Random rnd;
    int N=500; // number of blocks
    int L=10000; // number of throws in each block
    int dimension=3;
    int Nstep=10000; 
    double x_start=1.;
    double initial_increment=0.3;
    double target_ratio=0.5;
    double tolerance=0.001;
    double a=0;
    ofstream output;
    rnd.start();
    // create a function object
    Function* psi;
    if(ground) psi = new psi_100();
    else psi = new psi_210();
    norm2 r2;
    Metropolis metro(*psi,r_start,output,r2);
    metro.Set_spherical_coordinates(true);
    metro.Set_gaussian(gaussian);
    //set the length of the step for psi_100 & psi_210
    metro.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro.Set_string(stringname);
    if(equilibation) metro.Equilibrate(3e6,rnd);
    cout << "Metropolis algorithm for psi_100 uniform increment: a = " <<setprecision(4)<< metro.Get_a() << endl;
    metro.Average(N,L,filename);

}
int main(int argc, char* argv[]){

    // create a metropolis object
    vec r_start{100,100,100};
    vector<string> filename={"Results/ex_05.1_ground_uniform_far.dat","Results/ex_05.1_excited_uniform_far.dat","Results/ex_05.1_ground_gauss_far.dat","Results/ex_05.1_excited_gauss_far.dat"};
    vector<string> filename_eq={"Results/ex_05.1_ground_uniform_far_eq.dat","Results/ex_05.1_excited_uniform_far_eq.dat","Results/ex_05.1_ground_gauss_far_eq.dat","Results/ex_05.1_excited_gauss_far_eq.dat"};
    vector<string> stringname={"Results/RW_ground_uniform_far.dat","Results/RW_excited_uniform_far.dat","Results/RW_ground_gauss_far.dat","Results/RW_excited_gauss_far.dat"};
    vector<string> stringname_eq={"Results/RW_ground_uniform_far_eq.dat","Results/RW_excited_uniform_far_eq.dat","Results/RW_ground_gauss_far_eq.dat","Results/RW_excited_gauss_far_eq.dat"};
    vector<bool> gaussian={false,false,true,true};
    vector<bool> ground={true,false,true,false};
    for(int j=0;j<2;j++){
    for(int i=0;i<4;i++){
        //non equilibrate
        if(j==0)simulation_far(r_start,filename[i],i,gaussian[i],ground[i],stringname[i],false);
        //equilibrate
        else simulation_far(r_start,filename_eq[i],i,gaussian[i],ground[i],stringname_eq[i],true);
    }
    }
    return 0;
}
