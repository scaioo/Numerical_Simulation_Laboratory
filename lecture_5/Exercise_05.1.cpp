#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<armadillo>
#include <iomanip>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
#include "../classes/inc/block_average.h"
#include "../classes/inc/func.h"
using namespace std;
using namespace arma;
int main(int argc, char* argv[]){
    int N=500; // number of blocks
    int L=10000; // number of throws in each block
    int dimension=3;
    int Nstep=10000; 
    double x_start=1.;
    double initial_increment=0.3;
    double target_ratio=0.5;
    double tolerance=0.001;
    double a=0;
    // start random number generator
    Random rnd;
    rnd.start();
    // create a function object
    psi_100 psi; psi_210 psi2;
    norm2 r2;
    // create a metropolis object
    vec r_start{1,0,0};
    ofstream output;
    string filename_ground_unif="Results/ex_05.1_ground_uniform.dat",filename_excited_unif="Results/ex_05.1_excited_uniform.dat";
    string filename_ground_gauss="Results/ex_05.1_ground_gauss.dat",filename_excited_gauss="Results/ex_05.1_excited_gauss.dat";
    //------------------------------------------//
    //Create Metropolis for Uniform distribution
    //------------------------------------------//
    Metropolis metro_ground(psi,r_start,output,r2), metro_excited(psi2,r_start,output,r2);
    metro_ground.Set_gaussian(true);
    //set the spherical coordinates calculation
    metro_ground.Set_spherical_coordinates(true); //metro_excited_unif.Set_spherical_coordinates(true);
    
    //set the length of the step for psi_100 & psi_210
    metro_ground.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro_ground.Set_string("Results/RW_ground_gauss.dat");
    cout << "Metropolis algorithm for psi_100 uniform increment: a = " <<setprecision(4)<< metro_ground.Get_a() << endl;
    metro_excited.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro_excited.Set_string("Results/RW_excited_uniform.dat");
    cout << "Metropolis algorithm for psi_210 uniform increment: a = " <<setprecision(4)<< metro_excited.Get_a() << endl;
    metro_ground.Average(N,L,filename_ground_unif);
    metro_excited.Average(N,L,filename_excited_unif);

    //------------------------------------------//
    //   Metropolis for Gaussian distribution
    //------------------------------------------//
    //set the gaussian distribution for psi_100 & psi_210
    metro_ground.Set_gaussian(true); metro_excited.Set_gaussian(true);
    metro_ground.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro_ground.Set_string("Results/RW_ground_gauss.dat");
    cout << "Metropolis algorithm for psi_100 gaussian increment: a = " <<setprecision(3)<< metro_ground.Get_a() << endl; 
    metro_excited.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro_excited.Set_string("Results/RW_excited_gauss.dat");
    cout << "Metropolis algorithm for psi_210 gaussian increment: a = " <<setprecision(3)<< metro_excited.Get_a() << endl;
    
    // calculate Metropolis with a gaussian distribution for psi_100 & psi_210
    metro_ground.Set_r(r_start); metro_excited.Set_r(r_start);
    metro_ground.Average(N,L,filename_ground_gauss);
    metro_excited.Average(N,L,filename_excited_gauss);
    return 0;
}
