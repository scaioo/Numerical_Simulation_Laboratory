#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <numeric>
#include <vector>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
using namespace std;

int main(){
    //Load Random number generator with default seeds
    Random rnd;
    rnd.start();
    // Parameters declaration
    int M = 10000000;             // Total number of throws
    int N = 1000;                // Number of blocks
    int L = M / N;              // Number of throws in each block
    // saving results on output.txt
    string ofilename="Risultati/ex_01.1_01.2_output.txt";
    ofstream outputFile=openofile(ofilename);
    // Initializing first line
    outputFile << "Number_of_throws;\tIntegral_value;\tIntegral_error;\tsigma_value;\tsigma_error;\tChi^2;\tprob_poisson" << endl;
    // Variables declaration 
    int throws=0;
    // EX 01.1 variables
    double ave=0;
    double av2=0;
    double randomNumbers=0;
    double sum_prog=0;
    double su2_prog=0;
    double sum_prog_helper=0;
    double su2_prog_helper=0;
    double err_prog=0;
    double sum=0;
    // EX 01.1.2 variables
    double sigma_sum=0;
    double sigma_ave=0;
    double sigma_av2=0;
    double sigma_sum_prog_helper=0;
    double sigma_su2_prog_helper=0;
    double sigma_sum_prog=0;
    double sigma_su2_prog=0;
    double sigma_err_prog=0;
    // EX 01.1.3 variables
    double Ei=0.5;
    double Chi2=0;
    double ave_Chi2_helper=0;
    double ave_Chi2=0;
    double prob_poisson=0;

    // ------------------------------------------------------------------------
    //                        EXCERCIZES 01.1.1 & 01.1.2 
    // ------------------------------------------------------------------------

    // Integral and sigma^2 caluclation
    for(int i=0; i<N;i++){
        //calculation of uncertainty
        sum=0;
        sigma_sum=0;
        for(int j=0;j<L;j++){
        //generate random number
        randomNumbers=rnd.Rannyu();
        //sum of random numbers and sigma^2 generated 
        sum+=randomNumbers;
        sigma_sum+=(randomNumbers-0.5)*(randomNumbers-0.5);
        }
        // Computing the Integral Value with uncertainty
        ave= sum/L;
        av2=ave*ave;
        sum_prog_helper+=ave; // Sum j{0..i} r_j
        su2_prog_helper+=av2; // Sum j{0..i} (r_j)^2
        sum_prog=sum_prog_helper/(i+1); // Cumulative average
        su2_prog=su2_prog_helper/(i+1); // Cumulative square average
        err_prog=error(sum_prog,su2_prog,i); // Statistical uncertainty
        // Computing sigma^2 values with uncertainty
        sigma_ave= sigma_sum/L;
        sigma_av2=sigma_ave*sigma_ave;
        sigma_sum_prog_helper+=sigma_ave; // Sum j{0..i} r_j
        sigma_su2_prog_helper+=sigma_av2; // Sum j{0..i} (r_j)^2
        sigma_sum_prog=sigma_sum_prog_helper/(i+1); // Cumulative average
        sigma_su2_prog=sigma_su2_prog_helper/(i+1); // Cumulative square average
        sigma_err_prog=error(sigma_sum_prog,sigma_su2_prog,i); // Statistical uncertainty
        // Computing Chi^2
        Chi2=(sum_prog-Ei)*(sum_prog-Ei)/Ei;
        ave_Chi2_helper+=Chi2;
        ave_Chi2=ave_Chi2_helper/(i+1);
        prob_poisson=ave_Chi2/N;
        // saving results on output.txt
        throws=(i+1)*L;
        outputFile<<scientific<<setprecision(5)<<throws<<";\t"<<sum_prog<<";\t"<<err_prog<<";\t"<<sigma_sum_prog<<";\t"<<sigma_err_prog<<";\t"<<ave_Chi2<<";\t"<<prob_poisson<<endl;
    }
    outputFile.close();

    // ------------------------------------------------------------------------
    //                        EXCERCIZES 01.1.3
    // ------------------------------------------------------------------------
    
    //setting variables
    Chi2=0;
    int iteration=100; // iteration of Chi^2's calculation  
    N=100; 
    L=10000;
    M=N*L;
    throws=0;
    int counter=0;
    double Chi2_helper=0;
    // saving results on output.txt
    string ofilename_chi2="Results/ex_01_output.2.txt";
    ofstream outputFile_chi2=openofile(ofilename_chi2);
    // Initializing first line
    outputFile_chi2 << "Number_of_throws;\tchi2" << endl;
    for(int k=0;k<iteration;k++){
        //computing Chi^2
        for(int i=0;i<N;i++){
            for(int j=0;j<L;j++){
                randomNumbers=rnd.Rannyu();
                if((static_cast<float>(i)/N) <= randomNumbers && randomNumbers < (static_cast<float>(i + 1)/N)){
                    counter++;
                }
            }
            Chi2_helper=static_cast<double>((counter-(L/N))*(counter-(L/N)))/(L/N);      
            Chi2+=Chi2_helper;
            counter=0;
        }
        throws=(k+1);
        outputFile_chi2<<scientific<<setprecision(10)<<throws<<";\t"<<Chi2<<endl;
        Chi2=0;
    }
    outputFile_chi2.close();
    string ofilename_chi2_2="Results/ex_01_output.3.txt";
    ofstream outputFile_chi2_2=openofile(ofilename_chi2_2);
    // Initializing first line
    outputFile_chi2_2 << "Number_of_throws;\tchi2" << endl;
    iteration=10000;
    for(int k=0;k<iteration;k++){
        //computing Chi^2
        for(int i=0;i<N;i++){
            for(int j=0;j<L;j++){
                randomNumbers=rnd.Rannyu();
                if((static_cast<float>(i)/N) <= randomNumbers && randomNumbers < (static_cast<float>(i + 1)/N)){
                    counter++;
                }
            }
            Chi2_helper=static_cast<double>((counter-(L/N))*(counter-(L/N)))/(L/N);      
            Chi2+=Chi2_helper;
            counter=0;
        }
        throws=(k+1);
        outputFile_chi2_2<<scientific<<setprecision(10)<<throws<<";\t"<<Chi2<<endl;
        Chi2=0;
    }
    outputFile_chi2_2.close();
    return 0;
}