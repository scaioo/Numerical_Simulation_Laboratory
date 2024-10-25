    #include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <numeric>
#include <cmath>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
using namespace std;

int main(){
    //initialize the random number generator
    Random rnd;
    rnd.start();
    // saving results on output.txt
    string ofilename="Results/ex_02.1_integral.txt";
    ofstream outputFile;
    try{ outputFile=openofile(ofilename);} catch (const invalid_argument &e) {cerr << "Error: " << e.what() << endl;}
    //initialize first line 
    outputFile << "throws;\tuniform_Integral_value;\terror;\tTaylor_Integral_value;\terror" << endl;
    //initialize variables
    int N_throws=1000000, N_blocks=100;
    int L=N_throws/N_blocks;
    double rnd_number=0,Taylor_rnd_number_distribution=0;
    double function=0,Taylor_function=0;
    double sum_prog_helper=0, Taylor_sum_prog_helper=0, su2_prog_helper=0, Taylor_su2_prog_helper=0;
    double ave=0, av2=0, Taylor_ave=0, Taylor_av2=0;
    double sum_prog=0, su2_prog=0, Taylor_sum_prog=0, Taylor_su2_prog=0;
    double err_prog=0, Taylor_err_prog=0;
    double Px=0;
    //calculation of the integral
    for(int i=0;i<N_blocks;i++){
        function=0;
        Taylor_function=0;
        for(int j=0;j<L;j++){
            // calculation of the integral with uniform distribution
            rnd_number=rnd.Rannyu();
            function += M_PI/2*cos(M_PI*rnd_number/2);
            // calculation of the integral with Taylor distribution
            Taylor_rnd_number_distribution=1-sqrt(1-rnd.Rannyu());
            Px=2-2*Taylor_rnd_number_distribution;
            Taylor_function += M_PI/2*cos(M_PI*Taylor_rnd_number_distribution/2)/Px;
        }
            // uniform distribution block average
            ave=function/L;
            av2=ave*ave;
            sum_prog_helper+=ave;
            su2_prog_helper+=av2;
            sum_prog= sum_prog_helper/(i+1);
            su2_prog= su2_prog_helper/(i+1);
            err_prog=error(sum_prog,su2_prog,i);
            // Taylor distribution block average
            Taylor_ave=Taylor_function/L;
            Taylor_av2=Taylor_ave*Taylor_ave;
            Taylor_sum_prog_helper+=Taylor_ave;
            Taylor_su2_prog_helper+=Taylor_av2;
            Taylor_sum_prog= Taylor_sum_prog_helper/(i+1);
            Taylor_su2_prog= Taylor_su2_prog_helper/(i+1);
            Taylor_err_prog=error(Taylor_sum_prog,Taylor_su2_prog,i);

            // writing the results on the output file
            outputFile <<scientific<<setprecision(10)<< (i+1) << ";\t" << sum_prog << ";\t" << err_prog << ";\t" << Taylor_sum_prog << ";\t" << Taylor_err_prog << endl;
        
    }
    outputFile.close();
    rnd.SaveSeed();
    return 0;
}