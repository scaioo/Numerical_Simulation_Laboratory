#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <armadillo>
#include "../functions/functions.hpp"
#include "../random_number/random.h"
using namespace std;
using namespace arma;

int main(){

    //initializing random number generator
    Random rnd;
    rnd.start();

    //declaration of variables
    int N_rw=1e5;                               //total number of ramdom walks
    int N_step=100;                              //number of step for each random walk
    double a=1;                                 // lattice step
    double discRW_r=0.0, contRW_r=0.0;          //distance from the origin
    vec p_discrete={0.0,0.0,0.0};    // coordinates
    vec p_continuous={0.0,0.0,0.0};  // coordinates

    // variables for the blocking method for discrete (discRW) and continuous (contRW) random walks
    double discRW_ave=0.0,discRW_av2=0.0, contRW_ave=0.0, contRW_av2=0.0;
    double discRW_sum_prog=0.0,discRW_su2_prog=0.0, contRW_sum_prog=0.0, contRW_su2_prog=0.0;
    double discRW_sum_prog_helper=0.0,discRW_su2_prog_helper=0.0, contRW_sum_prog_helper=0.0, contRW_su2_prog_helper=0.0;
    double discRW_err_prog=0.0, contRW_err_prog=0.0;

    //open files to save the results
    string filename1="Results/ex_02.2_discrete_random_walks.txt";
    ofstream ofile1;
    string filename2="Results/ex_02.2_continuous_random_walks.txt";
    ofstream ofile2;
    try{ofile1=openofile(filename1);}catch(invalid_argument &e){cerr << e.what() << endl; exit(-1);}
    try{ofile2=openofile(filename2);}catch(invalid_argument &e){cerr << e.what() << endl; exit(-1);}
    //writing the first line of the file
    ofile1 << "N_step; sqrt(<r^2>); error" << endl;
    ofile2 << "N_step; sqrt(<r^2>); error" << endl;

    double discRW_r_mean=0.0, contRW_r_mean=0.0;
    int N_block=100;
    int L=N_rw/N_block;
    for(int k=0; k<N_step;k++){
        // set to zero the variables to compute the block average for each step
        discRW_sum_prog_helper=0.0, discRW_su2_prog_helper=0.0;
        contRW_sum_prog_helper=0.0, contRW_su2_prog_helper=0.0;
        for (int m=0;m<N_block;m++){
            discRW_r_mean=0.0;contRW_r_mean=0.0;
            for(int l=0;l<L;l++){
                //running the Random walk
                p_discrete={0.0,0.0,0.0};p_continuous={0.0,0.0,0.0};
                for(int j=0;j<k+1;j++){
                    p_discrete=discrete_walk(rnd,p_discrete,a);
                    p_continuous=continuous_walk(rnd,p_continuous,a);
                }
                // distance of the RW
                discRW_r=distance2(p_discrete);
                contRW_r=distance2(p_continuous);
                // sum of the mean distace for the RW of a block
                discRW_r_mean+=discRW_r;
                contRW_r_mean+=contRW_r;
            }
            // computing the mean distance for the RW of a block
            discRW_r_mean/=L;
            contRW_r_mean/=L;
            // computing block averages
            //cout <<  discRW_ave << endl;
            discRW_ave=sqrt(discRW_r_mean);
            discRW_av2=discRW_ave*discRW_ave;
            discRW_sum_prog_helper+=discRW_ave;
            discRW_su2_prog_helper+=discRW_av2;
            discRW_sum_prog=discRW_sum_prog_helper/(m+1);
            discRW_su2_prog=discRW_su2_prog_helper/(m+1);
            discRW_err_prog=error(discRW_sum_prog,discRW_su2_prog,m);
            // computing the average distance from the origin for the continuous random walk
            contRW_ave=sqrt(contRW_r_mean);
            contRW_av2=contRW_ave*contRW_ave;
            contRW_sum_prog_helper+=contRW_ave;
            contRW_su2_prog_helper+=contRW_av2;
            contRW_sum_prog=contRW_sum_prog_helper/(m+1);
            contRW_su2_prog=contRW_su2_prog_helper/(m+1);
            contRW_err_prog=error(contRW_sum_prog,contRW_su2_prog,m);

        }
        // write results on files
        ofile1 << k+1 << ";\t" << discRW_sum_prog << ";\t" << discRW_err_prog << endl;
        ofile2 << k+1 << ";\t" << contRW_sum_prog << ";\t" << contRW_err_prog << endl;
    }
    //closing the files
    ofile1.close();
    ofile2.close();
    //saving the seed
    rnd.SaveSeed();
    return 0;
}