#include "../inc/block_average.h"
#include "../../random_number/random.h"
#include "../../functions/functions.hpp"
#include "../inc/func.h"
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>
using namespace std;
using namespace arma;

//----------------------------------------------------------------
// CLASS TO COMPUTE BLOCK AVERAGE
//----------------------------------------------------------------
    //constructor
    BlockMeanCalculator::BlockMeanCalculator() {}
    // Metodo per calcolare la media a blocchi e la sua incertezza
    void BlockMeanCalculator::Average(const int N_Blocks,const int BlockSize,string filename) {
        Random rnd;
        rnd.start();
        // setting variables
        double ave=0; double av2=0;
        double sum_prog=0; double su2_prog=0;
        double sum_prog_helper=0; double su2_prog_helper=0;
        double err_prog=0; double sum=0.0; double value=0;
        // open file
        ofstream outputFile;
        try{ outputFile=openofile(filename);} catch (const invalid_argument &e) {cerr << "Error: " << e.what() << endl;}
        // write firs line of the file
        outputFile << "N;" << "\t mean value;" << "\t error" << endl;
        for (int i = 0; i < N_Blocks; ++i) {
            sum = 0.0;
            for (int j = 0; j < BlockSize; ++j) {
                sum += calculate(value,rnd);
            }
        ave= sum/BlockSize;
        av2=ave*ave;
        sum_prog_helper+=ave; // Sum j{0..i} r_j
        su2_prog_helper+=av2; // Sum j{0..i} (r_j)^2
        sum_prog=sum_prog_helper/(i+1); // Cumulative average
        su2_prog=su2_prog_helper/(i+1); // Cumulative square average
        err_prog=error(sum_prog,su2_prog,i); // Statistical uncertainty
        //if(err_prog==0 && i!=0) {throw invalid_argument("Errore nullo.");} // if the error is zero, the program stops
        outputFile << i+1 << " " << sum_prog << " " << err_prog << endl;       
        }
        _average={sum_prog,err_prog}; // save the average and its uncertainty
    }
    //default destructor
    BlockMeanCalculator::~BlockMeanCalculator() {}

//----------------------------------------------------------------    
// CALL OPTION PRICE EX 3.1
//----------------------------------------------------------------
    //constructor
    call_price::call_price(double S0, double T, double K, double r, double sigma,double mu,bool direct) : m_S0(S0), m_T(T), m_K(K), m_r(r), m_sigma(sigma), m_mu(mu), m_direct(direct) {}
    // method to set nstep
    void call_price::Set_Nstep(int step){m_Nstep=step;}
    // method to get nstep
    double call_price::Get_Nstep() {return m_Nstep;}
    // method to calculate the call option price
    double call_price::calculate(double value,Random &rnd) const{
        double St;
        if(m_direct){
            //here we use the direct method
            St=spot_price(m_S0,m_r,m_sigma,m_T,rnd.Gauss(m_mu,m_T));
        }else{
            if(m_Nstep==0){throw invalid_argument("calcolo discreto con N_step=0.");}
            //here we use the discretized method
            St=spot_price(m_S0,m_r,m_sigma,m_T,rnd,m_Nstep);
        }
        value=call_option_price(m_r,m_T,m_K,St);
        return value;
    } 
    //default destructor
    call_price::~call_price() {}

// -------------------------------------------
// PUT OPTION PRICE EX 3.1
// -------------------------------------------
    put_price::put_price(double S0, double T, double K, double r, double sigma,double mu,bool direct) : m_S0(S0), m_T(T), m_K(K), m_r(r), m_sigma(sigma), m_mu(mu), m_direct(direct) {}
    // method to set nstep
    void put_price::Set_Nstep(int step){m_Nstep=step;}
    // method to get nstep
    double put_price::Get_Nstep() {return m_Nstep;}
    // method to calculate the call option price
    double put_price::calculate(double value,Random &rnd) const{
        double St;
        if(m_direct){
            //here we use the direct method
            St=spot_price(m_S0,m_r,m_sigma,m_T,rnd.Gauss(m_mu,m_T));
        }else{
            if(m_Nstep==0){throw invalid_argument("calcolo discreto con N_step=0.");}
            //here we use the discretized method
            St=spot_price(m_S0,m_r,m_sigma,m_T,rnd,m_Nstep);
        }
        value=put_option_price(m_r,m_T,m_K,St);
        return value;
    }
    // default destructor
    put_price::~put_price() {}

// -------------------------------------------
// first Metropolis EX 5.1
// -------------------------------------------
    // Define the static member null_stream
    ofstream Metropolis::null_stream("/dev/null");
    // method to set nstep
    void Metropolis::Set_Nstep(int step){m_Nstep=step;}
    // method to get nstep
    double Metropolis::Get_Nstep() {return m_Nstep;}
    // set the value of the parameter a
    void Metropolis::Set_a(int Nstep,double x_start ,double initial_increment,double target_ratio,double tolerance,Random &rnd,int dimension){
        int Nhit=0;
        double a=x_start; double ratio=0; double increment=initial_increment;
        // Initialize vectors for the current position and new position
        vec x(dimension, arma::fill::zeros);
        vec xnew(dimension, arma::fill::zeros);
        vec rand_vec(dimension);
        // Spherical coordinate vectors
        vec coord(dimension);
        vec newcoord(dimension);
        for(;;){
            Nhit=0;
            for(int i=0;i<Nstep;i++){
                // Generate random displacements
                for(int d=0;d<dimension;d++){
                    if(m_gaussian){
                        rand_vec(d) = rnd.Gauss(0, a);
                    }else{
                        rand_vec(d) = rnd.Rannyu(-a, a);}
                }
                // Update the position
                xnew = x + rand_vec;
                // Convert to spherical coordinates if necessary
                if (m_spherical_coordinates) {
                    cartesianToSpherical(x, coord);
                    cartesianToSpherical(xnew, newcoord);
                } else {
                    coord = x;
                    newcoord = xnew;
                }
                // Calculate acceptance probability
                double p = pow(m_f(coord),2);
                double pnew = pow(m_f(newcoord),2);
                double A = min(1.,pnew/p);
                double r = rnd.Rannyu(0,1);
                // Accept or reject the new position
                if (r <= A) {Nhit++;x=xnew;}
                // Calculate the acceptance ratio
                ratio = static_cast<double>(Nhit)/Nstep;
            }
            // Check if the ratio is within the tolerance range
            if (fabs(ratio-target_ratio) < tolerance) {m_a=a;break;}
            // Adjust increment based on the difference
            if (ratio <= target_ratio) {a *= 0.9;}else{a *= 1.1;}
            if(a<0){a=increment;}
            }
    }
    // Get a
    double Metropolis::Get_a() {return m_a;}
    // method to calculate the function with metropolis
    double Metropolis::calculate(double value,Random &rnd) const{
            if(m_save_traj) _file.open(_string,std::ios::app);
            // Determine the dimensionality from m_r
            int dimension = m_r.n_elem;
            //cout << "dimension = " << dimension << endl;
            vec vec_rand(dimension);
            if (m_gaussian) {
                for (int i = 0; i < dimension; ++i) {
                    vec_rand(i) = rnd.Gauss(0, m_a);
                }
            } else {
                for (int i = 0; i < dimension; ++i) {
                    vec_rand(i) = rnd.Rannyu(-m_a, m_a);
                }
            }
            vec vec_rnew=m_r+vec_rand;
            // Handle spherical coordinates if necessary
            vec coord(dimension), newcoord(dimension);
            if(m_spherical_coordinates) {
                cartesianToSpherical(m_r, coord);
                cartesianToSpherical(vec_rnew, newcoord);
            } else {
                coord = m_r;
                newcoord = vec_rnew;
            }
            double p = pow(m_f(coord),2);
            double pnew = pow(m_f(newcoord),2);
            double A = min(1.,pnew/p);
            double accept = rnd.Rannyu(0,1);
            if (accept <= A) {
                    m_r = vec_rnew;
                    double value = m_g(m_r);
                    if (m_save_traj) {
                        for (int i = 0; i < dimension; ++i) {
                            _file << m_r(i) << (i == dimension - 1 ? "\n" : "\t");
                        }
                        _file.close();
                    }
                    return value;
                } else {
                    double value = m_g(m_r);
                    if (m_save_traj) _file.close();
                    return value;
                }
            }

    void Metropolis::Equilibrate(int Nstep,Random &rnd){
        // Determine the dimensionality from m_r
        for(int i=0;i<Nstep;i++){
            calculate(0,rnd);
        }
    }
    //default destructor
    Metropolis::~Metropolis() {;}

