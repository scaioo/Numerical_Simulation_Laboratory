#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <vector>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
#include <bits/stdc++.h>
using namespace std;

int main() {
    // declaring variables
    vector<double> unif_sum(4,0);       // uniform N=1,2,10,100
    vector<double> exp_sum(4,0);        // exponential N=1,2,10,100
    vector<double> lorentz_sum(4,0);    // lorentzian N=1,2,10,100
    int realization=10000;
    valarray<int> N = {1,2,10,100}; // number of repetitions
    // initializing random number generator
    Random rnd;
    rnd.start(); 
    // initializing output file
    string ofilename_unif="Results/ex_01.2_uniform.txt";
    string ofilename_exp="Results/ex_01.2_exp.txt";
    string ofilename_lorentz="Results/ex_01.2_lorentz.txt";
    ofstream oFile_uniform=openofile(ofilename_unif);
    ofstream oFile_exp=openofile(ofilename_exp);
    ofstream oFile_lorentz=openofile(ofilename_lorentz);
    // initializing first line
    oFile_uniform << "Uniform_N=1;\tUniform_N=2;\tUniform_N=10;\tUniform_N=100" << endl;
    oFile_exp << "Exp_N=1;\tExp_N=2;\tExp_N=10;\tExp_N=100" << endl;
    oFile_lorentz << "Lorentz_N=1;\tLorentz_N=2;\tLorentz_N=10;\tLorentz_N=100" << endl;
    // looping over the realizations
    for(int j=0;j<realization;j++){ 
        // for N=1,2,10,100
        for(size_t k=0;k<N.size();k++){
            //cout << "Realization: " << j << " N: " << N[k] << endl;
            for(int i=0;i<N[k];i++){
                unif_sum[k]+=(rnd.Rannyu()/N[k]);
                exp_sum[k]+=(rnd.Exponential(1)/N[k]);
                lorentz_sum[k]+=(rnd.CauchyLorentz(0,1)/N[k]);
                //cout << "unif sum: " << unif_sum[k] << " Exp: " << exp_sum[k] << " Lorentz: " << lorentz_sum[k] << endl;
            }
        }
        // writing on file
        oFile_uniform << scientific << setprecision(5) << unif_sum[0] << ";\t" << unif_sum[1] << ";\t" << unif_sum[2] << ";\t" << unif_sum[3] << endl;
        oFile_exp << scientific << setprecision(5) << exp_sum[0] << ";\t" << exp_sum[1] << ";\t" << exp_sum[2] << ";\t" << exp_sum[3] << endl;
        oFile_lorentz << scientific << setprecision(5) << lorentz_sum[0] << ";\t" << lorentz_sum[1] << ";\t" << lorentz_sum[2] << ";\t" << lorentz_sum[3] << endl;
        // resetting the sum
        unif_sum.clear();
        exp_sum.clear();
        lorentz_sum.clear();
        // resizing the vector
        unif_sum.resize(4,0);
        exp_sum.resize(4,0);
        lorentz_sum.resize(4,0);
    }
    // closing the file
    oFile_uniform.close();
    oFile_exp.close();
    oFile_lorentz.close();

    return 0;
}