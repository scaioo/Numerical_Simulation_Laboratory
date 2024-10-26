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
    vec r_start{100,100,100};
    ofstream output;
    Metropolis metro_ground(psi,r_start,output,r2), metro_excited(psi2,r_start,output,r2);
    //set the spherical coordinates calculation
    metro_ground.Set_spherical_coordinates(true); metro_excited.Set_spherical_coordinates(true);
    string filename_ground_unif="Results/ex_05.1_ground_uniform_far_eq.dat",filename_excited_unif="Results/ex_05.1_excited_uniform_far_eq.dat";
    //set the length of the step for psi_100 & psi_210
    metro_ground.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    //Equilibration
    metro_ground.Equilibrate(3e6,rnd);
    metro_ground.Set_string("Results/RW_ground_uniform_far_eq.dat");
    cout << "Metropolis algorithm for psi_100 uniform increment: a = " <<setprecision(4)<< metro_ground.Get_a() << endl;
    metro_excited.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    //Equilibration
    metro_excited.Equilibrate(3e6,rnd);
    metro_excited.Set_string("Results/RW_excited_uniform_far_eq.dat");
    cout << "Metropolis algorithm for psi_210 uniform increment: a = " <<setprecision(4)<< metro_excited.Get_a() << endl;
    //calculate metropolis with a uniform distribution for psi_100 & psi_210
    metro_ground.Average(N,L,filename_ground_unif);
    metro_excited.Average(N,L,filename_excited_unif);
    //calculate metropolis with a gaussian distribution for psi_100 & psi_210
    string filename_ground_gauss="Results/ex_05.1_ground_gauss_far_eq.dat",filename_excited_gauss="Results/ex_05.1_excited_gauss_far_eq.dat";
    //set the gaussian distribution for psi_100 & psi_210
    metro_ground.Set_gaussian(true);metro_excited.Set_gaussian(true);
    metro_ground.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro_ground.Set_string("Results/RW_ground_gauss_far_eq.dat");
    metro_ground.Set_spherical_coordinates(true);
    cout << "Metropolis algorithm for psi_100 gaussian increment: a = " <<setprecision(3)<< metro_ground.Get_a() << endl;    
    metro_excited.Set_a(Nstep,x_start,initial_increment,target_ratio,tolerance,rnd,dimension);
    metro_excited.Set_string("Results/RW_excited_gauss_far_eq.dat");
    cout << "Metropolis algorithm for psi_210 gaussian increment: a = " <<setprecision(3)<< metro_excited.Get_a() << endl;
    metro_ground.Average(N,L,filename_ground_gauss);
    metro_excited.Average(N,L,filename_excited_gauss);
    return 0;
}